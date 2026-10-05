# Roman De Santos (rdesantos@hmc.edu)
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025-2026


# how to run: python3 <path to this script> <path to target lcao.neb_geom> <file prefix>
# This program converts the .neb_geom into seperate files labled <file prefix>n.xyz where n is the image number in the .neb_geom file
# Note a user must update the atom_types manualy at the moment for each different NEB work, so it is recomended to make distinct copies of this file 
# for distinct NEB jobs.


import sys
import os

def convert_bohr_to_angstrom(input_path, output_prefix):
    # Standard conversion factor
    BOHR_TO_ANG = 0.52917721067
    
    # Atom sequence (Must Match NEB Atom order)

    atom_types = ["C", "H","H","H","H"]

    # atom_types = [
    #     "O", "Si", "Si", "C", "C", "C", "C", "C", "C", "C", "C", 
    #     "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", 
    #     "H", "H", "H", "H", "H", "H", "H", "Si", "O", "O", "Si", 
    #     "Si", "C", "C", "C", "H", "C", "C", "C", "C", "H", "H", 
    #     "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", 
    #     "H", "H", "H", "H", "H", "H", "H", "H"    
    # ]

    # atom_types = [
    #    "O", "Si", "Si", "C", "C", "H", "H", "H", "Si", "O", "O", "Si", "Si", "H", "H", "H", "H",
    #     "H", "H", "H", "H", "H", "H", "H", "H", "H"  
    # ]
    num_atoms = len(atom_types)
    
    if not os.path.exists(input_path):
        print(f"Error: Input file {input_path} not found.")
        return

    current_image_coords = []
    block_count = 0

    with open(input_path, 'r') as f:
        for line in f:
            parts = line.split()
            
            # Process lines that contain exactly 3 numerical coordinates
            if len(parts) == 3:
                try:
                    x = float(parts[0]) * BOHR_TO_ANG
                    y = float(parts[1]) * BOHR_TO_ANG
                    z = float(parts[2]) * BOHR_TO_ANG
                    
                    # Match coordinate to atom type index
                    atom_symbol = atom_types[len(current_image_coords)]
                    current_image_coords.append(f"{atom_symbol:<2} {x:15.10f} {y:15.10f} {z:15.10f}")
                    
                    # Once we have a full block of 63 atoms, write it to a file
                    if len(current_image_coords) == num_atoms:
                        # Determine tag based on your specific order: einitial, efinal, then images
                        if block_count == 0:
                            tag = "einitial"
                        elif block_count == 1:
                            tag = "efinal"
                        else:
                            tag = f"image{block_count - 1}"
                        
                        # Create filename (e.g., ConvertedImages_image1.xyz)
                        filename = f"{output_prefix}_{tag}.xyz"
                        
                        with open(filename, 'w') as out_f:
                            out_f.write(f"{num_atoms}\n")
                            out_f.write(f"{tag}\n")
                            for coord in current_image_coords:
                                out_f.write(f"{coord}\n")
                        
                        print(f"Created: {filename}")
                        
                        # Reset for next block
                        current_image_coords = []
                        block_count += 1
                        
                except ValueError:
                    # Skip lines that aren't coordinate data (headers/energy strings)
                    continue

    print(f"Done. Processed {block_count} total blocks.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 NEBimages2xyz.py <input_file> <output_prefix>")
    else:
        # sys.argv[2] now acts as a prefix for the individual filenames
        convert_bohr_to_angstrom(sys.argv[1], sys.argv[2])