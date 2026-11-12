#!/usr/bin/env python3
"""Post-process a finished SeqQuest BTO slab relaxation.

Run from the simulation folder:
    python3 postprocess_slab.py                (add --cleanup to tidy the folder)

This adapts the team's postprocess.sh, geom2xyz.sh and hist2xyz.py workflow.
Unlike postprocess.sh, it reads the scale, cell and run settings from THIS run's
lcao.in, so it can be used on any slab size.

Creates runN/ containing:
    lcao_runN.in/.out/.hist/.geom   copies of the main run files
    lcao_runN.xyz                   final structure in Angstroms
    lcao_runN_trajectory.xyz        all geometry steps in one file
    steps/stepN.xyz                 one file per step, for animate_xyz.py
    lcao_next.in                    this run's lcao.in with final coordinates
    force_convergence.csv/.png      force history (plot needs matplotlib)
    postprocess_report.txt          convergence and structure summary

Original files are not changed unless --cleanup is used. --cleanup moves files
to trashdir/runN/ (nothing is deleted) and replaces lcao.in with lcao_next.in.
Requires only standard Python 3.6+.
"""
import argparse
import csv
import glob
import math
import os
import re
import shutil
import subprocess
import sys

BOHR_TO_ANG = 0.529177
RY_PER_BOHR_TO_EV_PER_ANG = 13.60569 / BOHR_TO_ANG
RY_TO_EV = 13.60569
FORCE_RE = re.compile(
    r"rms,max force=\s*([-+0-9.DdEe]+)\s+([-+0-9.DdEe]+)\s+for component#\s*(\d+)"
)
ATOM_HEADER_RE = re.compile(r"atom,\s*type", re.I)


def fortran_float(text):
    return float(text.replace("D", "E").replace("d", "e"))


def is_number(text):
    try:
        fortran_float(text)
        return True
    except ValueError:
        return False


def element_name(label):
    match = re.match(r"[A-Za-z]+", label)
    return match.group(0).capitalize() if match else label


def read_lines(path):
    with open(path, "r", errors="replace") as handle:
        return handle.read().splitlines()


def next_data_line(lines, index):
    index += 1
    while index < len(lines) and not lines[index].strip():
        index += 1
    return index


def parse_lcao_in(path):
    lines = read_lines(path)
    settings = {
        "lines": lines,
        "scale": 1.0,
        "scale_found": False,
        "coordinates": "CARTESIAN",
        "functional": "not found",
        "dimension": "not found",
        "kgrid": "not found",
        "vectors": None,
        "natoms": None,
        "atom_header": None,
        "atom_start": None,
        "atom_labels": [],
        "atom_files": [],
        "geometry_block": [],
    }

    for index, raw in enumerate(lines):
        key = raw.strip().lower()
        value_index = next_data_line(lines, index)
        value = lines[value_index].strip() if value_index < len(lines) else ""
        if key == "scale":
            settings["scale"] = float(value.split()[0].replace("D", "E"))
            settings["scale_found"] = True
        elif key == "coordinates":
            settings["coordinates"] = value.upper()
        elif key == "functional":
            settings["functional"] = value
        elif key.startswith("dimension of system"):
            settings["dimension"] = value
        elif key == "kgrid":
            settings["kgrid"] = value
        elif key == "primitive lattice vectors":
            vectors = []
            cursor = index
            for _ in range(3):
                cursor = next_data_line(lines, cursor)
                vectors.append([float(x) for x in lines[cursor].split()[:3]])
            settings["vectors"] = vectors
        elif key == "number of atoms in unit cell":
            settings["natoms"] = int(value.split()[0])
            header = next_data_line(lines, value_index)
            settings["atom_header"] = header
            settings["atom_start"] = header + 1
            for line in lines[header + 1:header + 1 + settings["natoms"]]:
                settings["atom_labels"].append(line.split()[1])
        elif key == "atom file" and "=" in value:
            settings["atom_files"].append(value.split("=")[0].strip())
        elif key == "geometry relaxation":
            for line in lines[index:]:
                settings["geometry_block"].append(line.strip())
                if line.strip().lower() == "end geometry relaxation":
                    break

    if settings["natoms"] is None or settings["vectors"] is None:
        sys.exit("Could not find the atom list or lattice vectors in lcao.in.")
    return settings


def input_units_to_ang(settings):
    return BOHR_TO_ANG * settings["scale"]


def cell_lengths_ang(settings):
    factor = input_units_to_ang(settings)
    vectors = settings["vectors"]
    return [vectors[0][0] * factor, vectors[1][1] * factor, vectors[2][2] * factor]


def is_orthogonal(settings):
    vectors = settings["vectors"]
    off_diagonal = [vectors[i][j] for i in range(3) for j in range(3) if i != j]
    return all(abs(value) < 1e-8 for value in off_diagonal)


def resolve_label(label, settings, index):
    if label.isdigit():
        type_index = int(label) - 1
        if 0 <= type_index < len(settings["atom_files"]):
            return element_name(settings["atom_files"][type_index])
    if index < len(settings["atom_labels"]):
        fallback = settings["atom_labels"][index]
        if label.isdigit() and not fallback.isdigit():
            return element_name(fallback)
    return element_name(label)


def parse_geom(path, settings):
    """Read the last atom block in lcao.geom, in lcao.in coordinate units."""
    blocks = [[]]
    for line in read_lines(path):
        if ATOM_HEADER_RE.search(line):
            blocks.append([])
            continue
        parts = line.split()
        if len(parts) == 5 and parts[0].isdigit() and all(is_number(x) for x in parts[2:]):
            label, coords = parts[1], parts[2:]
        elif len(parts) >= 4 and all(is_number(x) for x in parts[1:4]):
            label, coords = parts[0], parts[1:4]
        else:
            continue
        blocks[-1].append((label, [fortran_float(x) for x in coords]))

    atoms = next((block for block in reversed(blocks) if block), [])
    if len(atoms) != settings["natoms"]:
        sys.exit(
            "lcao.geom has {} atoms but lcao.in has {}.".format(
                len(atoms), settings["natoms"]
            )
        )
    labels = [resolve_label(label, settings, i) for i, (label, _) in enumerate(atoms)]
    return labels, [coords for _, coords in atoms]


def parse_hist(path):
    lines = read_lines(path)
    header = {"cell": None, "types": [], "natoms": None, "type_indices": []}
    frames = []
    index = 0
    current = None

    while index < len(lines):
        line = lines[index].strip()
        upper = line.upper()
        if upper.startswith("@CELL VECTORS"):
            rows = [lines[index + offset].split() for offset in (1, 2, 3)]
            header["cell"] = [
                fortran_float(rows[0][0]),
                fortran_float(rows[1][1]),
                fortran_float(rows[2][2]),
            ]
            index += 4
            continue
        if upper.startswith("@TYPES"):
            parts = line.split()
            if len(parts) > 1 and parts[1].isdigit():
                count = int(parts[1])
            else:
                index = next_data_line(lines, index)
                count = int(lines[index].strip())
            for _ in range(count):
                index = next_data_line(lines, index)
                header["types"].append(lines[index].split()[0])
            index += 1
            continue
        if upper.startswith("@NUMBER OF ATOMS"):
            index = next_data_line(lines, index)
            header["natoms"] = int(lines[index].split()[0])
            while len(header["type_indices"]) < header["natoms"]:
                index = next_data_line(lines, index)
                header["type_indices"].extend(int(x) for x in lines[index].split())
            index += 1
            continue
        if upper.startswith("@GSTEP"):
            index = next_data_line(lines, index)
            current = {"step": lines[index].strip(), "coords": [], "energy": None}
            frames.append(current)
            index += 1
            continue
        if upper.startswith("@COORDINATES") and current is not None:
            for _ in range(header["natoms"]):
                index = next_data_line(lines, index)
                current["coords"].append(
                    [fortran_float(x) for x in lines[index].split()[:3]]
                )
            index += 1
            continue
        if "ENERGY" in upper and upper.startswith("@") and current is not None:
            value_index = next_data_line(lines, index)
            if value_index < len(lines):
                parts = lines[value_index].split()
                if parts and is_number(parts[0]):
                    current["energy"] = fortran_float(parts[0])
        index += 1

    frames = [frame for frame in frames if len(frame["coords"]) == header["natoms"]]
    return header, frames


def parse_forces(lines):
    history = []
    for line in lines:
        match = FORCE_RE.search(line)
        if match:
            rms = fortran_float(match.group(1))
            maximum = fortran_float(match.group(2))
            history.append(
                {
                    "rms_ry": rms,
                    "max_ry": maximum,
                    "rms_ev": rms * RY_PER_BOHR_TO_EV_PER_ANG,
                    "max_ev": maximum * RY_PER_BOHR_TO_EV_PER_ANG,
                    "component": int(match.group(3)),
                }
            )
    return history


def find_final_energy(lines):
    preferred = re.compile(r"(final|relaxed).*energy|energy.*(final|relaxed)", re.I)
    fallback = re.compile(r"total energy", re.I)
    for pattern in (preferred, fallback):
        matches = [line.strip() for line in lines if pattern.search(line)]
        matches = [line for line in matches if re.search(r"-?\d+\.\d+", line)]
        if matches:
            line = matches[-1]
            value = float(re.findall(r"-?\d+\.\d+(?:[DdEe][-+]?\d+)?", line)[-1]
                          .replace("D", "E").replace("d", "e"))
            return line, value
    return None, None


def minimum_image(delta, length):
    return delta - length * round(delta / length)


def xy_distance(first, second, cell):
    dx = minimum_image(first[0] - second[0], cell[0])
    dy = minimum_image(first[1] - second[1], cell[1])
    dz = first[2] - second[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def find_layers(labels, coords, gap):
    order = sorted(range(len(coords)), key=lambda i: coords[i][2])
    layers = [[order[0]]]
    for previous, current in zip(order, order[1:]):
        if coords[current][2] - coords[previous][2] > gap:
            layers.append([])
        layers[-1].append(current)
    return layers


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def layer_summary(labels, coords, layers):
    rows = []
    for number, layer in enumerate(layers, start=1):
        counts = {}
        for index in layer:
            counts[labels[index]] = counts.get(labels[index], 0) + 1
        cation_z = [coords[i][2] for i in layer if labels[i] != "O"]
        oxygen_z = [coords[i][2] for i in layer if labels[i] == "O"]
        all_z = [coords[i][2] for i in layer]
        rows.append(
            {
                "number": number,
                "composition": " ".join(
                    "{}{}".format(element, counts[element]) for element in sorted(counts)
                ),
                "z": mean(all_z),
                "rumpling": mean(cation_z) - mean(oxygen_z)
                if cation_z and oxygen_z
                else float("nan"),
            }
        )
    return rows


def ti_o_bonds(labels, coords, cell, cutoff):
    bonds = []
    oxygen = [i for i, label in enumerate(labels) if label == "O"]
    for titanium in [i for i, label in enumerate(labels) if label == "Ti"]:
        lengths = sorted(
            distance
            for distance in (xy_distance(coords[titanium], coords[o], cell) for o in oxygen)
            if distance <= cutoff
        )
        if lengths:
            bonds.append((titanium + 1, len(lengths), lengths[0], lengths[-1]))
    return bonds


def closest_pair(coords, cell):
    best = (float("inf"), None, None)
    for first in range(len(coords)):
        for second in range(first + 1, len(coords)):
            distance = xy_distance(coords[first], coords[second], cell)
            if distance < best[0]:
                best = (distance, first, second)
    return best


def write_xyz(path, labels, coords, comment):
    with open(path, "w") as handle:
        handle.write("{}\n{}\n".format(len(coords), comment))
        for label, (x, y, z) in zip(labels, coords):
            handle.write("{} {:.10f} {:.10f} {:.10f}\n".format(label, x, y, z))


def write_next_input(path, settings, coords):
    lines = settings["lines"]
    start = settings["atom_start"]
    stop = start + settings["natoms"]
    atom_lines = []
    # Keep the labels exactly as lcao.in wrote them (element names or type numbers).
    for number, (label, position) in enumerate(zip(settings["atom_labels"], coords), start=1):
        atom_lines.append(
            " {} {} {:14.10f} {:14.10f} {:14.10f}".format(number, label, *position)
        )
    with open(path, "w") as handle:
        handle.write("\n".join(lines[:start] + atom_lines + lines[stop:]) + "\n")


def next_run_number():
    numbers = [
        int(path[3:])
        for path in glob.glob("run*")
        if os.path.isdir(path) and path[3:].isdigit()
    ]
    return max(numbers + [0]) + 1


def latest_job_id():
    ids = []
    for path in glob.glob("slurm-*.out"):
        match = re.fullmatch(r"slurm-(\d+)\.out", path)
        if match:
            ids.append(int(match.group(1)))
    return str(max(ids)) if ids else None


def job_status(job_id):
    if not job_id or not shutil.which("sacct"):
        return "sacct unavailable"
    command = ["sacct", "-X", "-n", "-j", job_id, "--format=State,ExitCode,Elapsed"]
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            universal_newlines=True,
            timeout=20,
        )
        return " ".join(result.stdout.split()) or "no sacct record"
    except (OSError, subprocess.SubprocessError):
        return "sacct unavailable"


def plot_forces(path, history, target):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    steps = list(range(1, len(history) + 1))
    figure, axis = plt.subplots(figsize=(7, 4.5))
    axis.semilogy(steps, [row["max_ev"] for row in history], marker="o", label="Maximum force")
    axis.semilogy(steps, [row["rms_ev"] for row in history], marker="s", label="RMS force")
    axis.axhline(target, color="black", linestyle="--", label="Target")
    axis.set_xlabel("Geometry step")
    axis.set_ylabel("Force (eV/Å)")
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return True


def cleanup(run_number):
    trash = os.path.join("trashdir", "run{}".format(run_number))
    os.makedirs(trash, exist_ok=True)
    keep_suffixes = (".atm", ".sh", ".job", ".py")
    this_script = os.path.basename(__file__)
    for item in os.listdir("."):
        if item in ("trashdir", "lcao.in", this_script):
            continue
        if os.path.isdir(item) and re.fullmatch(r"run\d+", item):
            continue
        if item.endswith(keep_suffixes):
            continue
        shutil.move(item, os.path.join(trash, item))


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--job-id", help="Slurm job ID (default: newest slurm-*.out)")
    parser.add_argument("--ftol", type=float, default=0.05,
                        help="Maximum force target in eV/Å (default: 0.05)")
    parser.add_argument("--layer-gap", type=float, default=1.0,
                        help="z gap in Å that separates atomic layers (default: 1.0)")
    parser.add_argument("--ti-o-cutoff", type=float, default=2.6,
                        help="Ti-O bond cutoff in Å (default: 2.6)")
    parser.add_argument("--cleanup", action="store_true",
                        help="Move extra files to trashdir and install lcao_next.in as lcao.in")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    for required in ("lcao.in", "lcao.out", "lcao.geom"):
        if not os.path.isfile(required):
            sys.exit("Missing {}. Run this from the simulation folder.".format(required))

    settings = parse_lcao_in("lcao.in")
    if settings["coordinates"] != "CARTESIAN":
        sys.exit("Only CARTESIAN coordinates are supported, found: {}".format(settings["coordinates"]))
    if not is_orthogonal(settings):
        sys.exit("This script only supports orthogonal cells.")

    output_lines = read_lines("lcao.out")
    forces = parse_forces(output_lines)
    labels, geom_coords = parse_geom("lcao.geom", settings)
    factor = input_units_to_ang(settings)
    cell = cell_lengths_ang(settings)
    final_coords = [[value * factor for value in row] for row in geom_coords]

    warnings = []
    hist_header, frames = ({}, [])
    if os.path.isfile("lcao.hist"):
        try:
            hist_header, frames = parse_hist("lcao.hist")
        except (ValueError, IndexError):
            warnings.append("lcao.hist has an unexpected format, so no trajectory was made.")
    hist_factor = None
    if frames and hist_header.get("cell"):
        hist_factor = cell[0] / hist_header["cell"][0]
        if hist_header["natoms"] != settings["natoms"]:
            warnings.append("lcao.hist and lcao.in list different numbers of atoms.")
            frames = []
    elif os.path.isfile("lcao.hist") and not warnings:
        warnings.append("Could not read coordinates from lcao.hist.")

    job_id = arguments.job_id or latest_job_id()
    run_number = next_run_number()
    run_dir = "run{}".format(run_number)
    prefix = os.path.join(run_dir, "lcao_run{}".format(run_number))
    os.makedirs(run_dir)

    for extension in ("in", "out", "hist", "geom"):
        source = "lcao.{}".format(extension)
        if os.path.isfile(source):
            shutil.copy2(source, "{}.{}".format(prefix, extension))
    for path in glob.glob("*.job") + glob.glob("slurm-{}.*".format(job_id)):
        shutil.copy2(path, run_dir)

    write_xyz(prefix + ".xyz", labels, final_coords, "Final geometry from lcao.geom (Angstrom)")
    write_next_input(os.path.join(run_dir, "lcao_next.in"), settings, geom_coords)

    initial_coords = None
    if frames:
        steps_dir = os.path.join(run_dir, "steps")
        os.makedirs(steps_dir)
        hist_labels = [
            element_name(hist_header["types"][index - 1])
            for index in hist_header["type_indices"]
        ]
        half_cell = [length / 2.0 for length in hist_header["cell"]]
        with open(prefix + "_trajectory.xyz", "w") as trajectory:
            for number, frame in enumerate(frames, start=1):
                converted = []
                for x, y, z in frame["coords"]:
                    converted.append(
                        [
                            (x - hist_header["cell"][0] if x > half_cell[0] else x) * hist_factor,
                            (y - hist_header["cell"][1] if y > half_cell[1] else y) * hist_factor,
                            (z - half_cell[2]) * hist_factor,
                        ]
                    )
                comment = "GSTEP {} (Angstrom)".format(frame["step"])
                write_xyz(os.path.join(steps_dir, "step{}.xyz".format(number)),
                          hist_labels, converted, comment)
                trajectory.write("{}\n{}\n".format(len(converted), comment))
                for label, (x, y, z) in zip(hist_labels, converted):
                    trajectory.write("{} {:.10f} {:.10f} {:.10f}\n".format(label, x, y, z))
        initial_coords = [
            [value * hist_factor for value in row] for row in frames[0]["coords"]
        ]
    if initial_coords is None:
        initial_coords = [
            [float(x) * factor for x in line.split()[2:5]]
            for line in settings["lines"][
                settings["atom_start"]:settings["atom_start"] + settings["natoms"]
            ]
        ]

    with open(os.path.join(run_dir, "force_convergence.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "rms_Ry_per_bohr", "max_Ry_per_bohr",
                         "rms_eV_per_A", "max_eV_per_A", "max_component"])
        for step, row in enumerate(forces, start=1):
            writer.writerow([step, row["rms_ry"], row["max_ry"], row["rms_ev"],
                             row["max_ev"], row["component"]])
    plotted = bool(forces) and plot_forces(
        os.path.join(run_dir, "force_convergence.png"), forces, arguments.ftol
    )

    report = []
    add = report.append
    add("SEQQUEST SLAB POST-PROCESSING REPORT")
    add("Run folder: {}   Slurm job: {}   sacct: {}".format(run_dir, job_id, job_status(job_id)))
    error_file = "slurm-{}.err".format(job_id)
    if os.path.isfile(error_file):
        size = os.path.getsize(error_file)
        add("{}: {}".format(error_file, "empty" if size == 0 else "{} bytes, read it".format(size)))

    add("")
    add("INPUT SETTINGS FROM lcao.in")
    add("Functional: {}   Dimension: {}   k-grid: {}".format(
        settings["functional"], settings["dimension"], settings["kgrid"]))
    add("Scale: {} ({})   Cell: {:.4f} x {:.4f} x {:.4f} Å".format(
        settings["scale"],
        "found" if settings["scale_found"] else "NOT found, using 1.0",
        cell[0], cell[1], cell[2]))
    for line in settings["geometry_block"]:
        add("  " + line)

    add("")
    add("CONVERGENCE")
    if forces:
        last = forces[-1]
        add("Geometry steps with forces: {}".format(len(forces)))
        add("Force steps in lcao.out / frames in lcao.hist: {} / {}".format(len(forces), len(frames)))
        add("Start max force: {:.4f} eV/Å".format(forces[0]["max_ev"]))
        add("Final max force: {:.4f} eV/Å ({:.3E} Ry/bohr)".format(last["max_ev"], last["max_ry"]))
        add("Final RMS force: {:.4f} eV/Å".format(last["rms_ev"]))
        best = min(range(len(forces)), key=lambda i: forces[i]["max_ev"]) + 1
        add("Lowest max force: {:.4f} eV/Å at step {}".format(forces[best - 1]["max_ev"], best))
        met = last["max_ev"] <= arguments.ftol
        add("Target {:.3f} eV/Å: {}".format(arguments.ftol, "MET" if met else "NOT MET"))
        if not met:
            warnings.append("Final maximum force is above the target.")
    else:
        warnings.append("No 'rms,max force' lines were found in lcao.out.")
        add("No force history found.")

    energy_line, energy_value = find_final_energy(output_lines)
    add("")
    add("ENERGY")
    if energy_line:
        add("Line used: " + energy_line)
        add("Value: {} (as printed)".format(energy_value))
        if "Ry" in energy_line and not re.search(r"\beV\b", energy_line):
            add("Equivalent: {:.6f} eV".format(energy_value * RY_TO_EV))
        add("Confirm this is the final relaxed energy before using it as E_slab.")
    else:
        warnings.append("No final energy line was found in lcao.out.")
        add("Not found. See the end of lcao.out below.")
    frame_energies = [frame["energy"] for frame in frames if frame["energy"] is not None]
    if len(frame_energies) >= 2:
        add("Energy change in lcao.hist, first to last step: {:.6f}".format(
            frame_energies[-1] - frame_energies[0]))

    counts = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    add("")
    add("STRUCTURE")
    add("Atoms: {}   Composition: {}".format(
        len(labels), " ".join("{}{}".format(k, counts[k]) for k in sorted(counts))))
    distance, first, second = closest_pair(final_coords, cell)
    add("Closest pair: {:.3f} Å ({}{} - {}{})".format(
        distance, labels[first], first + 1, labels[second], second + 1))
    if distance < 1.4:
        warnings.append("Closest atoms are under 1.4 Å. Check the units and geometry.")

    displacements = []
    for atom in range(len(labels)):
        displacements.append(xy_distance(final_coords[atom], initial_coords[atom], cell))
    shift = [mean([final_coords[i][axis] - initial_coords[i][axis]
                   for i in range(len(labels))]) for axis in (2,)]
    if abs(shift[0]) > 1.0:
        adjusted = [[row[0], row[1], row[2] + shift[0]] for row in initial_coords]
        displacements = [xy_distance(final_coords[i], adjusted[i], cell) for i in range(len(labels))]
        add("Removed a uniform z offset of {:.3f} Å between lcao.geom and lcao.hist.".format(shift[0]))
    largest = max(range(len(labels)), key=lambda i: displacements[i])
    add("Displacement from first step: RMS {:.3f} Å, maximum {:.3f} Å ({}{})".format(
        math.sqrt(mean([d * d for d in displacements])), displacements[largest],
        labels[largest], largest + 1))
    still = [i + 1 for i, d in enumerate(displacements) if d < 1e-4]
    add("Atoms that did not move: {}".format(len(still)))

    layers = find_layers(labels, final_coords, arguments.layer_gap)
    final_rows = layer_summary(labels, final_coords, layers)
    initial_rows = layer_summary(labels, initial_coords, layers)
    add("")
    add("LAYERS (bottom to top; cation minus oxygen z is the rumpling)")
    add("{:>5} {:<16} {:>10} {:>13} {:>13} {:>12}".format(
        "Layer", "Composition", "z (Å)", "Rumple (Å)", "Start (Å)", "Spacing (Å)"))
    for number, row in enumerate(final_rows):
        spacing = ""
        if number:
            spacing = "{:.3f}".format(row["z"] - final_rows[number - 1]["z"])
        add("{:>5} {:<16} {:>10.3f} {:>13.3f} {:>13.3f} {:>12}".format(
            row["number"], row["composition"], row["z"], row["rumpling"],
            initial_rows[number]["rumpling"], spacing))
    add("Bottom termination: {}   Top termination: {}".format(
        final_rows[0]["composition"], final_rows[-1]["composition"]))

    bonds = ti_o_bonds(labels, final_coords, cell, arguments.ti_o_cutoff)
    add("")
    add("Ti-O BONDS (cutoff {:.2f} Å)".format(arguments.ti_o_cutoff))
    if bonds:
        add("Shortest: {:.3f} Å   Longest: {:.3f} Å".format(
            min(b[2] for b in bonds), max(b[3] for b in bonds)))
        add("{:>6} {:>12} {:>12} {:>12}".format("Ti", "Neighbors", "Min (Å)", "Max (Å)"))
        for atom, neighbors, shortest, longest in bonds:
            add("{:>6} {:>12} {:>12.3f} {:>12.3f}".format(atom, neighbors, shortest, longest))

    add("")
    add("WARNINGS")
    report.extend(["- " + warning for warning in warnings] or ["None"])
    add("")
    add("LAST 15 LINES OF lcao.out")
    report.extend(output_lines[-15:])
    add("")
    add("Force plot: {}".format("force_convergence.png" if plotted else "not made (needs matplotlib)"))

    report_path = os.path.join(run_dir, "postprocess_report.txt")
    with open(report_path, "w") as handle:
        handle.write("\n".join(report) + "\n")
    print("\n".join(report))
    print("\nSaved everything in {}/".format(run_dir))

    if arguments.cleanup:
        shutil.copy2(os.path.join(run_dir, "lcao_next.in"), "lcao.in.new")
        cleanup(run_number)
        os.replace(os.path.join("trashdir", "run{}".format(run_number), "lcao.in.new"), "lcao.in")
        print("Moved extra files to trashdir/run{}/ and updated lcao.in.".format(run_number))


if __name__ == "__main__":
    main()