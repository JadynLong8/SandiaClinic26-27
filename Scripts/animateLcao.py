#!/usr/bin/env python3
"""
animate_lcao.py
Extracts atomic coordinates and types from a SeqQuest lcao.hist file, 
converts from Ångström to Bohr, and creates an animated 3D scatter plot 
with bond lines, saved as a GIF.

NOTE: The Matplotlib backend is set to 'Agg' for safe execution on HPC systems 
(like Bridges-2) where no graphical display is available.
"""

import re
import sys
import numpy as np
import matplotlib

# CRUCIAL: Set the non-GUI backend for HPC execution
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation, PillowWriter

# Conversion constant
ANGSTROM_TO_BOHR = 1.0 / 0.529177  # ≈ 1.8897261

# Define a simple color map for common elements (CPK coloring).
# You can add more elements as needed.
ATOM_COLOR_MAP = {
    'H': 'white',      # Hydrogen
    'C': 'black',      # Carbon
    'N': 'blue',       # Nitrogen
    'O': 'red',        # Oxygen
    'P': 'orange',     # Phosphorus
    'Si': 'gray',      # Silicon
    'Pt': 'darkgray',  # Platinum
    'DEFAULT': 'purple'  # Color for any other atom
}

# Define a global bond threshold in Ångströms.
# This is a simple approximation. Any two atoms closer than this
# distance will be drawn with a bond.
# *** You may need to adjust this value for your system! ***
BOND_THRESHOLD_ANGSTROM = 2.0
BOND_THRESHOLD_BOHR = BOND_THRESHOLD_ANGSTROM * ANGSTROM_TO_BOHR


def parse_lcao_hist(filename):
    """
    Extract coordinates and atom types for each @GSTEP block in lcao.hist.
    This version reads atom types from the file's header.

    Returns a tuple: (all_coords, atom_types)
    - all_coords: A list of NumPy arrays, one for each step.
    - atom_types: A single list of atom type strings (e.g., ['Pt', 'C', 'H', ...])
    """
    try:
        with open(filename, "r") as f:
            content = f.read()
    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        sys.exit(1)

    all_coords = []
    atom_types = []

    # --- Step A - Parse Atom Types from header ---
    try:
        # 1. Get the type names (e.g., ['O', 'H', 'C', ...])
        # Looks for the block between @TYPES and @CELL VECTORS
        types_match = re.search(r"@TYPES\s+\d+\s+([\s\S]*?)@CELL VECTORS", content)
        if not types_match:
            print("Error: Could not find '@TYPES ... @CELL VECTORS' block.")
            sys.exit(1)
            
        type_names = types_match.group(1).strip().splitlines()
        # -> type_names = ['O', 'H', 'C', 'Pt', 'P', 'Si']

        # 2. Get the index mapping (e.g., [1, 6, 6, 3, ...])
        # Looks for the block between @NUMBER OF ATOMS and @TYPES
        atoms_match = re.search(r"@NUMBER OF ATOMS\s+\d+\s+([\s\S]*?)@TYPES", content)
        if not atoms_match:
            print("Error: Could not find '@NUMBER OF ATOMS ... @TYPES' block.")
            sys.exit(1)
            
        # Flatten all numbers into one list
        type_indices_str = atoms_match.group(1).strip().split()
        type_indices = [int(i) for i in type_indices_str]
        # -> type_indices = [1, 6, 6, 3, 3, ..., 4]

        # 3. Create the final atom_types list by mapping indices to names
        # We use (i - 1) because the file is 1-indexed but Python lists are 0-indexed
        atom_types = [type_names[i - 1] for i in type_indices]
    
    except Exception as e:
        print(f"Error parsing file header: {e}")
        print("Please check your file's @TYPES and @NUMBER OF ATOMS format.")
        sys.exit(1)

    # --- Step B - Parse Coordinates for each GSTEP ---
    gsteps = content.split("@GSTEP")[1:]
    if not gsteps:
        print("Error: No @GSTEP blocks found in file.")
        sys.exit(1)

    for block_num, block in enumerate(gsteps):
        match = re.search(r"@COORDINATES\s+([\s\S]*?)@ESCF", block)
        if not match:
            print(f"Warning: No @COORDINATES ... @ESCF block found in GSTEP {block_num + 1}")
            continue

        coords_text = match.group(1).strip()
        coords = []
        for line in coords_text.splitlines():
            parts = line.split()
            
            # Now we only expect 3 parts: X, Y, Z
            if len(parts) == 3: 
                try:
                    x, y, z = map(float, parts)
                    coords.append([x * ANGSTROM_TO_BOHR,
                                   y * ANGSTROM_TO_BOHR,
                                   z * ANGSTROM_TO_BOHR])
                except ValueError:
                    # Skip lines that don't contain three valid floats
                    continue
        
        if coords:
            # Sanity check: make sure number of coords matches number of atoms
            if len(coords) == len(atom_types):
                all_coords.append(np.array(coords))
            else:
                print(f"Warning: GSTEP {block_num + 1} has {len(coords)} coordinates, "
                      f"but header defined {len(atom_types)} atoms. Skipping this step.")

    return all_coords, atom_types


def find_bonds(coords, threshold):
    """
    Finds pairs of atoms (by index) that are within the threshold distance.
    `coords` is a single snapshot (N_atoms, 3) array.
    `threshold` is the distance in the same units as coords (Bohr).
    Returns a list of tuples: [(0, 1), (0, 2), (1, 5), ...]
    """
    bonds = []
    n_atoms = len(coords)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):  # Start from i+1 to avoid double-counting
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance < threshold:
                bonds.append((i, j))
    print(f"🧬 Found {len(bonds)} bonds based on {BOND_THRESHOLD_ANGSTROM} Å threshold.")
    return bonds


def animate_coordinates(all_coords, atom_types, output_gif="lcao_animation.gif"):
    """
    Create and save a 3D scatter animation of atomic positions with bonds.
    Each atom TYPE will have a unique, consistent color.
    """
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    
    # Set explicit axis labels
    ax.set_xlabel("X (Bohr)")
    ax.set_ylabel("Y (Bohr)")
    ax.set_zlabel("Z (Bohr)")
    
    # --- Setup for type-based colors and bonds ---
    
    # 1. Get initial data and find bonds
    initial_data = all_coords[0]
    bonds = find_bonds(initial_data, BOND_THRESHOLD_BOHR)

    # 2. Map atom types to colors using the dictionary
    colors = [ATOM_COLOR_MAP.get(atom_type, ATOM_COLOR_MAP['DEFAULT']) 
              for atom_type in atom_types]

    # 3. Draw initial scatter plot (atoms)
    #    We add edgecolors='black' so white 'H' atoms are visible
    scat = ax.scatter(initial_data[:, 0], initial_data[:, 1], initial_data[:, 2], 
                      s=60, c=colors, edgecolors='black', lw=0.5)

    # 4. Draw initial bond lines
    bond_lines = []
    for i, j in bonds:
        atom1 = initial_data[i]
        atom2 = initial_data[j]
        # Using '-' for solid line and c='gray' for bond color (fixes warning)
        line, = ax.plot([atom1[0], atom2[0]], 
                        [atom1[1], atom2[1]], 
                        [atom1[2], atom2[2]], '-', lw=1.0, c='gray') 
        bond_lines.append(line)

    # --- End new setup ---

    # Fix plot limits based on overall min/max
    all_points = np.concatenate(all_coords)
    min_vals = all_points.min(axis=0)
    max_vals = all_points.max(axis=0)
    
    # Set buffer for prettier plot
    buffer = 2.0 
    ax.set_xlim(min_vals[0] - buffer, max_vals[0] + buffer)
    ax.set_ylim(min_vals[1] - buffer, max_vals[1] + buffer)
    ax.set_zlim(min_vals[2] - buffer, max_vals[2] + buffer)
    
    # ax.set_aspect('equal', 'box') is NOT supported in 3D plots

    def update(frame):
        data = all_coords[frame]
        
        # Update the 3D scatter plot data (atoms)
        scat._offsets3d = (data[:, 0], data[:, 1], data[:, 2])
        ax.set_title(f"GSTEP Frame {frame + 1}/{len(all_coords)} (Bohr units)")
        
        # Update all the bond lines
        for k, (i, j) in enumerate(bonds):
            atom1 = data[i]
            atom2 = data[j]
            bond_lines[k].set_data_3d([atom1[0], atom2[0]], 
                                      [atom1[1], atom2[1]], 
                                      [atom1[2], atom2[2]])
            
        # Return all modified artists (atoms + bonds)
        return [scat] + bond_lines

    # Create the animation object
    ani = FuncAnimation(fig, update, frames=len(all_coords),
                        interval=200, blit=False, repeat=True)
    
    # Save the animation using PillowWriter (requires pillow to be installed)
    ani.save(output_gif, writer=PillowWriter(fps=5))
    print(f"✅ Animation saved as {output_gif}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python animate_lcao.py <lcao.hist>")
        sys.exit(1)
    
    filename = sys.argv[1]
    print(f"📂 Reading file: {filename}")

    # Get both coords and atom types
    all_coords, atom_types = parse_lcao_hist(filename)
    
    print(f"🔢 Found {len(all_coords)} geometry steps")
    print(f"⚛️ Found {len(atom_types)} atoms in the system")

    if not all_coords or not atom_types:
        print("⚠️ No valid coordinate blocks found in file. Exiting.")
        sys.exit(1)

    # Generate output file name based on input
    output_gif = filename.replace(".hist", "_animation.gif") if ".hist" in filename else f"{filename}_animation.gif"
    
    # Pass atom_types to the animation function
    animate_coordinates(all_coords, atom_types, output_gif)

if __name__ == "__main__":
    main()