#!/usr/bin/env python3
"""Make a smooth Jmol movie from the slab trajectory.

Run from the simulation folder:
    python3 fix_trajectory.py

Reads  run1/lcao_run1_trajectory.xyz  and  run1/lcao_run1.xyz
Writes run1/lcao_run1_movie.xyz

Atoms that sit near half the cell can jump to the other side of the slab between
frames. This puts every atom in every frame on the same periodic image as the
final structure and removes the constant z offset. It does not change any
distances between atoms, apart from choosing the periodic image. The input files
are not changed.
"""
import argparse


def read_frames(path):
    with open(path) as handle:
        lines = handle.read().splitlines()
    frames = []
    index = 0
    while index < len(lines) and lines[index].strip():
        natoms = int(lines[index])
        comment = lines[index + 1]
        atoms = []
        for line in lines[index + 2:index + 2 + natoms]:
            parts = line.split()
            atoms.append((parts[0], [float(value) for value in parts[1:4]]))
        if len(atoms) != natoms:
            raise SystemExit("{} ends partway through a frame.".format(path))
        frames.append((comment, atoms))
        index += natoms + 2
    return frames


def nearest_image(value, reference, length):
    return value - length * round((value - reference) / length)


def main():
    parser = argparse.ArgumentParser(description="Make a smooth Jmol movie.")
    parser.add_argument("--run", default="run1", help="Run folder (default: run1)")
    parser.add_argument("--cell", type=float, default=8.0654,
                        help="x and y cell length in Angstroms (default: 8.0654)")
    arguments = parser.parse_args()

    run = arguments.run
    prefix = "{}/lcao_{}".format(run, run)
    frames = read_frames(prefix + "_trajectory.xyz")
    final = read_frames(prefix + ".xyz")[-1][1]

    last_atoms = frames[-1][1]
    if len(last_atoms) != len(final):
        raise SystemExit("The trajectory and final structure have different atom counts.")
    for (first_label, _), (second_label, _) in zip(last_atoms, final):
        if first_label != second_label:
            raise SystemExit("The trajectory and final structure list atoms in a different order.")

    z_offset = sum(a[1][2] - b[1][2] for a, b in zip(final, last_atoms)) / len(final)
    largest_shift = 0.0
    output = prefix + "_movie.xyz"
    with open(output, "w") as handle:
        for comment, atoms in frames:
            handle.write("{}\n{}\n".format(len(atoms), comment))
            for (label, (x, y, z)), (_, reference) in zip(atoms, final):
                new_x = nearest_image(x, reference[0], arguments.cell)
                new_y = nearest_image(y, reference[1], arguments.cell)
                new_z = z + z_offset
                largest_shift = max(
                    largest_shift,
                    abs(new_x - reference[0]),
                    abs(new_y - reference[1]),
                    abs(new_z - reference[2]),
                )
                handle.write("{} {:.10f} {:.10f} {:.10f}\n".format(label, new_x, new_y, new_z))

    print("Wrote {} with {} frames.".format(output, len(frames)))
    print("Largest distance from the final position along x, y or z: {:.3f} A".format(largest_shift))
    if largest_shift > 1.0:
        print("WARNING: that is larger than expected. Check the movie and tell Claude.")


if __name__ == "__main__":
    main()