#!/usr/bin/env python3
import os
import glob
import re

def alphanum_key(s):
    """Splits a string into numeric and non-numeric parts for natural sorting"""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def read_hist_file(filename):
    """
    Reads a multi-step .hist file.
    Returns a dictionary with common header information and a list of steps.
    Each step is a dictionary with:
      - step: the simulation step number (from @GSTEP)
      - coordinates: list of (x, y, z) floats for that step
    The header includes:
      - cell: (cell_x, cell_y, cell_z)
      - types: list of atom type strings (index 0 corresponds to type index 1)
      - num_atoms: total number of atoms
      - atom_type_indices: list of integers indicating the type for each atom
    """
    header = {}
    steps = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    i = 0
    # Read header info until we hit the first @GSTEP
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith('@CELL VECTORS'):
            # Next three lines are cell vectors
            cell_lines = [lines[i+1].strip(), lines[i+2].strip(), lines[i+3].strip()]
            cell_x = float(cell_lines[0].split()[0])
            cell_y = float(cell_lines[1].split()[1])
            cell_z = float(cell_lines[2].split()[2])
            header["cell"] = (cell_x, cell_y, cell_z)
            i += 4
            continue
        if line.startswith('@TYPES'):
            parts = line.split()
            try:
                num_types = int(parts[1])
            except:
                num_types = int(lines[i+1].strip())
                i += 1
            types = []
            for j in range(num_types):
                types.append(lines[i+1+j].strip())
            header["types"] = types
            i += 1 + num_types
            continue
        if line.startswith('@NUMBER OF ATOMS'):
            i += 1
            # Skip any empty lines until we get the number of atoms
            while i < len(lines) and not lines[i].strip():
                i += 1
            header["num_atoms"] = int(lines[i].strip())
            i += 1
            # Now read the atom type indices
            type_indices = []
            while len(type_indices) < header["num_atoms"]:
                parts = lines[i].strip().split()
                type_indices.extend([int(x) for x in parts])
                i += 1
            header["atom_type_indices"] = type_indices
            continue
        # Once we see the first @GSTEP, we assume header is done
        if line.startswith('@GSTEP'):
            break
        i += 1

    # Now process each simulation step
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('@GSTEP'):
            # Next non-empty line is the step number
            i += 1
            while i < len(lines) and not lines[i].strip():
                i += 1
            step_val = lines[i].strip()
            try:
                step_val = int(step_val)
            except:
                pass
            i += 1
            # Advance until we see @COORDINATES
            while i < len(lines) and not lines[i].strip().startswith('@COORDINATES'):
                i += 1
            if i >= len(lines):
                break
            i += 1  # skip the @COORDINATES line
            coords = []
            # Read as many lines as there are atoms
            for _ in range(header["num_atoms"]):
                # Skip any empty lines
                while i < len(lines) and not lines[i].strip():
                    i += 1
                if i >= len(lines):
                    break
                parts = lines[i].strip().split()
                if len(parts) < 3:
                    i += 1
                    continue
                x, y, z = map(float, parts[:3])
                coords.append((x, y, z))
                i += 1
            steps.append({"step": step_val, "coordinates": coords})
        else:
            i += 1

    return header, steps

def convert_coords(hist_coords, cell, scale_factor=4.0327222580):
    """
    Unwraps the coordinates from the .hist file
    For x and y: if coordinate > (cell_length / 2), subtract cell_length
    For z: subtract half the cell_z
    Returns a list of converted (x, y, z)
    """
    cell_x, cell_y, cell_z = cell
    half_cell_x = cell_x / 2.0
    half_cell_y = cell_y / 2.0
    half_cell_z = cell_z / 2.0
    conv = []
    for (x, y, z) in hist_coords:
        new_x = x - cell_x if x > half_cell_x else x
        new_y = y - cell_y if y > half_cell_y else y
        new_z = z - half_cell_z
        conv.append((new_x, new_y, new_z))
    return conv

def write_xyz(filename, conv_coords, atom_type_indices, types, step_label):
    """
    Writes a single .xyz file
    First line: number of atoms
    Second line: a comment line including the step label
    Then, one line per atom: atom label and converted coordinates
    """
    with open(filename, 'w') as f:
        f.write(f"{len(conv_coords)}\n")
        f.write(f"{step_label}\n")
        for idx, (x, y, z) in enumerate(conv_coords):
            # Atom type indices are 1-indexed
            atom_label = types[atom_type_indices[idx]-1] if 0 <= atom_type_indices[idx]-1 < len(types) else "X"
            f.write(f"{atom_label} {x: .8f} {y: .8f} {z: .8f}\n")

def process_all_hist_files():
    """
    Reads all .hist files in the current folder, extracts all GSTEP blocks,
    and writes them as consecutive .xyz files named step1.xyz, step2.xyz, ... stepN.xyz
    If run again, existing step*.xyz files will be removed so they are overwritten
    """
    # Remove any existing step*.xyz files
    old_files = glob.glob("step*.xyz")
    for file in old_files:
        os.remove(file)
        print(f"Removed old file: {file}")
    
    # Get a list of all .hist files
    hist_files = glob.glob("*.hist")
    if not hist_files:
        print("No .hist files found.")
        return
    
    # Use natural sorting to ensure proper numerical order
    hist_files = sorted(hist_files, key=alphanum_key)
    
    global_counter = 1
    # Loop over each .hist file
    for hist_filename in hist_files:
        print(f"Processing {hist_filename} ...")
        header, steps = read_hist_file(hist_filename)
        cell = header["cell"]
        types = header["types"]
        atom_type_indices = header["atom_type_indices"]
        
        for step_data in steps:
            step_val = step_data["step"]
            hist_coords = step_data["coordinates"]
            conv_coords = convert_coords(hist_coords, cell)
            # Create output filename using a global counter
            out_filename = f"step{global_counter}.xyz"
            # Optionally include the original file and step info in the comment
            step_label = f"File: {hist_filename}, GSTEP: {step_val}"
            write_xyz(out_filename, conv_coords, atom_type_indices, types, step_label)
            print(f"Wrote {out_filename}")
            global_counter += 1

if __name__ == '__main__':
    process_all_hist_files()

