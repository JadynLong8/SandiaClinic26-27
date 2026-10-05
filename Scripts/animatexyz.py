import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def alphanum_key(s):
    """Splits a string into numeric and non-numeric parts for natural sorting"""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def load_xyz(filename):
    """
    Loads an .xyz file and returns a NumPy array of coordinates and a list of element symbols.
    Assumes the first line is the number of atoms and the second is a comment
    """
    coords = []
    elements = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    try:
        natoms = int(lines[0].strip())
    except (ValueError, IndexError):
        # Fallback if the first line isn't a valid number or file is too short
        print(f"Warning: Could not parse atom count from {filename}. Attempting to read all lines.")
        natoms = len(lines) - 2 # fallback
        if natoms < 0:
            return np.empty((0, 3)), [] # Return empty if file is corrupt
            
    # Read exactly natoms lines starting from line 3 (index 2)
    for line in lines[2:2+natoms]:
        parts = line.split()
        if len(parts) >= 4:
            element = parts[0]
            try:
                x, y, z = map(float, parts[1:4])
            except ValueError:
                print(f"Warning: Skipping malformed line in {filename}: {line.strip()}")
                continue
            coords.append([x, y, z])
            elements.append(element)
            
    # Ensure coords is a 2D array even if empty
    if coords:
        coords_array = np.array(coords)
    else:
        coords_array = np.empty((0, 3))
    return coords_array, elements

def compute_bonds(coords, cutoff_rest=2.1, cutoff_last=1.8, cutoff_cross=None):
    """
    Compute bonds for a set of atom coordinates
    
    For atoms, the bond cutoff is chosen based on the group:
      - Atoms in the "last 27" of the array use cutoff_last
      - Atoms not in the last 27 use cutoff_rest
      - For bonds between atoms from different groups, cutoff_cross is used
    
    If cutoff_cross is not provided, it defaults to the smaller of cutoff_rest and cutoff_last
    
    Returns a list of bond pairs, where each pair is a tuple of two coordinate arrays
    """
    if cutoff_cross is None:
        cutoff_cross = min(cutoff_rest, cutoff_last)
    
    bonds = []
    n_atoms = coords.shape[0]
    # Define indices for the last 27 atoms
    last_group_indices = set(range(max(0, n_atoms - 27), n_atoms))
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            distance = np.linalg.norm(coords[i] - coords[j])
            # Determine which cutoff to use
            if i in last_group_indices and j in last_group_indices:
                cutoff = cutoff_last
            elif i not in last_group_indices and j not in last_group_indices:
                cutoff = cutoff_rest
            else:
                cutoff = cutoff_cross
            if distance < cutoff:
                bonds.append((coords[i], coords[j]))
    return bonds

def main():
    # Set your bond cutoff values
    cutoff_rest = 2.2   # For bonds among atoms not in the last 27
    cutoff_last = 2     # For bonds among atoms in the last 27
    cutoff_cross = min(cutoff_rest, cutoff_last)  # For bonds between the two groups
    
    # Get list of all .xyz files in the current directory
    xyz_files = glob.glob("*.xyz")
    if not xyz_files:
        print("No .xyz files found!")
        return

    # Sort the files in natural human order
    xyz_files = sorted(xyz_files, key=alphanum_key)
    print("Files found (in order):")
    for f in xyz_files:
        print(f)

    # Define the color mapping for each element
    color_map = {
        "Ti": "blue",
        "O": "red",
        "Ba": "green",
        "C": "brown",
        "H": "orange",
        "Si": "pink",
        "Pt": "purple"
    }

    trajectories = []    # List to store coordinate arrays for each frame
    file_names = []      # List to store the corresponding file names
    element_list = None  # To store the element symbols (assumed consistent across frames)

    for filename in xyz_files:
        coords, elems = load_xyz(filename)
        if coords.shape[0] == 0:
            print(f"Skipping empty or invalid file: {filename}")
            continue
            
        trajectories.append(coords)
        file_names.append(filename)
        if element_list is None:
            element_list = elems

    if not trajectories:
        print("No valid trajectory data was loaded. Exiting.")
        return

    # Build a list of colors for each atom using the element symbols.
    # If an element is not in the mapping, default to gray
    atom_colors = [color_map.get(elem, "gray") for elem in element_list]

    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set axis limits based on all coordinates from all frames
    all_coords = np.concatenate(trajectories)
    ax.set_xlim(np.min(all_coords[:, 0]), np.max(all_coords[:, 0]))
    ax.set_ylim(np.min(all_coords[:, 1]), np.max(all_coords[:, 1]))
    ax.set_zlim(np.min(all_coords[:, 2]), np.max(all_coords[:, 2]))

    # Create a scatter plot for the first frame with per-atom colors
    scat = ax.scatter(trajectories[0][:, 0], trajectories[0][:, 1], trajectories[0][:, 2],
                      s=50, c=atom_colors)

    # Global list to store bond line objects so we can remove them in each update
    bond_lines = []

    def update(frame):
        nonlocal bond_lines
        coords = trajectories[frame]
        # Update scatter plot positions
        scat._offsets3d = (coords[:, 0], coords[:, 1], coords[:, 2])
        # Remove previous bond lines
        for line in bond_lines:
            line.remove()
        bond_lines.clear()

        # Compute bonds for current coordinates using the adjusted cutoffs
        bonds = compute_bonds(coords, cutoff_rest, cutoff_last, cutoff_cross)
        # Draw each bond as a line
        for (atom1, atom2) in bonds:
            line, = ax.plot([atom1[0], atom2[0]],
                            [atom1[1], atom2[1]],
                            [atom1[2], atom2[2]], color='k', lw=1)
            bond_lines.append(line)

        # Set the title to show current frame and file name
        ax.set_title(f"Frame {frame+1}: {file_names[frame]}")
        return [scat] + bond_lines

    # Increase animation speed by reducing the interval (e.g., 100ms)
    ani = FuncAnimation(fig, update, frames=len(trajectories), interval=100, blit=False)

    # --- THIS IS THE CORRECTED PART ---
    # Instead of plt.show(), save the animation to a file.
    # This is required on a headless server like Bridges-2.
    print("Saving animation to 'molecule_animation.mp4'...")
    
    # You may need to run 'module load ffmpeg' on Bridges-2 before this script
    ani.save('molecule_animation.mp4', writer='ffmpeg', dpi=150)
    
    print("...Animation saved successfully.")
    # plt.show() # We comment this out as it won't work

if __name__ == "__main__":
    main()