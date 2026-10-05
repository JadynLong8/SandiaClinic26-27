#!/bin/bash
#SBATCH --job-name=postprocess
#SBATCH --output=postprocess.out
#SBATCH --error=postprocess.err
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Luis Lorenzana (HMC 25') [Initially for Slabs]
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025
# Email: llorenzana@g.hmc.edu

# Roman De Santos (HMC 26') [Adapted script for tBuPa clusters (dimension 0)]
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025-26
# Email: rdesantos@hmc.edu

# Postprocess script:
# Inputs: lcao.in and *.atm
# Outputs:
# Top Level Directory: 
# - lcao.in (with updated coordinate positions)
# - Will keep *.atm, *.sh, *.job, everything else will go in trashdir 
# - runN Directory: Everytime you run postprocess.sh; N+1
# -- lcao_runN.in (Original lcao.in file)
# -- lcao_runN.out (Output file for N run)
# -- lcao_runN.hist (History file with coordinates and forces at every G step)
# -- lcao_runN.geom (Last updated coordinate postion of structure)
# -- lcao_runN.xyz (Visualization file to see strecture in VESTA or Jmol)

###############################################################################
########################## Determine new run number ###########################
###############################################################################
max=0
for d in run*; do
    if [ -d "$d" ]; then
        num="${d#run}"
        # Only consider numeric directories
        if [[ "$num" =~ ^[0-9]+$ ]]; then
            if [ "$num" -gt "$max" ]; then
                max="$num"
            fi
        fi
    fi
done
run_num=$((max+1))
run_dir="run${run_num}"
mkdir -p "$run_dir"
echo "Using run number: ${run_num} (directory: ${run_dir})"

###############################################################################
######### Archive old lcao.in, lcao.out, and lcao.xyz (if present) ############
###############################################################################
archive_file () {
    local fname="$1" # original file name
    local ext="$2"   # extension (e.g., in, out, xyz)
    if [ -f "$fname" ]; then
        newname="lcao_run${run_num}.${ext}"
        mv "$fname" "$newname"
        mv "$newname" "$run_dir/"
        echo "Archived $fname as $newname into ${run_dir}/"
    fi
}

# lcao.geom and lcao.hist will be archived later (needed for processing)
archive_file "lcao.in" "in"
archive_file "lcao.out" "out"
archive_file "lcao.xyz" "xyz"

###############################################################################
##################### Find the .geom file to process ##########################
###############################################################################
geom_file_count=$(ls -1 *.geom 2>/dev/null | wc -l)
if [ "$geom_file_count" -eq 0 ]; then
    echo "No .geom file found in the current directory. Exiting."
    exit 1
elif [ "$geom_file_count" -gt 1 ]; then
    echo "Multiple .geom files found. Using the first one found."
fi
INPUT=$(ls *.geom | head -n 1)
echo "Processing input file: $INPUT"

###############################################################################
########## Generate new lcao.in file (remains in current directory) ###########
###############################################################################

############### YOU MAY NEED TO CHANGE THIS TO FIT YOUR SYSTEM ################


# Count the number of atom lines in the .geom file
ATOM_COUNT=$(awk '($0 !~ /atom,/) && NF>=4 { count++ } END { print count }' "$INPUT")

# Write the fixed header portion to new lcao.in
cat << 'EOF' > "lcao.in"
do setup
do iters
do force
do relax
do post
setup data
title
title2
functional
 PBE-SP
spin polarization
 1.0
dimension of system (0=cluster ... 3=bulk)
 0
coordinates
CARTESIAN
scale
 1
primitive lattice vectors
  60.0000000000     0.00000000000     0.00000000000
  0.00000000000     60.0000000000     0.00000000000
  0.00000000000     0.00000000000     60.0000000000
grid dimensions
       96  96  96
atom types
 4
atom file
 O = O.atm
atom file
 H = H.atm
atom file
 C = C.atm
atom file
 P = P.atm
number of atoms in unit cell
EOF

###############################################################################

# Append the computed atom count.
printf "  %d\n" "$ATOM_COUNT" >> "lcao.in"

# Append atom block header.
echo "atom, type, position;  step#     0" >> "lcao.in"

# Process the .geom file and append the atom data.
awk -v fmt="   %-5s  %-2s   %14.10f   %14.10f   %14.10f\n" '
/atom, type, position;/ { next }
NF == 5 {
    printf fmt, $1, $2, $3, $4, $5;
}
' "$INPUT" >> "lcao.in"

# Append the fixed footer.
cat << 'EOF' >> "lcao.in"
origin offset
 0 0 0
end setup phase data
run phase input data
history
 8
no ges
states
 200
temperature
 0.01
iterations
 300
blend ratio
0.300
convergence criterion
 0.02
geometry relaxation
gmethod
asd
gsteps
 60
end geometry relaxation
end run phase data
EOF

echo "New lcao.in file generated in the current directory."

###############################################################################
########################## Generate new lcao.xyz file #########################
###############################################################################

# Conversion factor: bohr to angstrom times scale (0.529177)
FACTOR=0.529177

# Count non-empty lines (excluding header lines that contain "atom, type, position;")
NUM_ATOMS=$(grep -v 'atom, type, position;' "$INPUT" | grep -cv '^\s*$')

# Create temporary lcao.xyz file.
tmp_xyz="tmp_lcao.xyz"
{
    echo "$NUM_ATOMS"
    echo ".geom to .xyz and scaled by factor of $FACTOR"
    awk -v factor="$FACTOR" '
    {
        if ($0 ~ /atom, type, position;/) next;
        if (NF==5 && $1 ~ /^[0-9]+$/) {
            symbol = $2;
            x = $3 * factor;
            y = $4 * factor;
            z = $5 * factor;
        } else if (NF>=4) {
        symbol = $2;
        x = $3 * factor;
        y = $4 * factor;
        z = $5 * factor;
        } else {
            next;
        }
        printf "%s %.10f %.10f %.10f\n", symbol, x, y, z;
    }' "$INPUT"
} > "$tmp_xyz"

# Rename new lcao.xyz as lcao_runN.xyz and move it into runN.
new_xyz="lcao_run${run_num}.xyz"
mv "$tmp_xyz" "$new_xyz"
mv "$new_xyz" "$run_dir/"
echo "New lcao.xyz file generated as ${run_dir}/${new_xyz}."

# Archive lcao.geom and lcao.hist (if present) into runN.
if [ -f "lcao.geom" ]; then
    newname="lcao_run${run_num}.geom"
    mv "lcao.geom" "$newname"
    mv "$newname" "$run_dir/"
    echo "Archived lcao.geom as $newname into ${run_dir}/"
fi

if [ -f "lcao.hist" ]; then
    newname="lcao_run${run_num}.hist"
    mv "lcao.hist" "$newname"
    mv "$newname" "$run_dir/"
    echo "Archived lcao.hist as $newname into ${run_dir}/"
fi

###############################################################################
############################### Cleanup Section ###############################
###############################################################################
echo "Starting cleanup..."

# In the cleanup, we want to leave in the current directory:
# - All *.atm files
# - All *.sh files
# - All *.job files
# - The new lcao.in file
# - All runN directories (which contain archived old files)
# Everything else will be moved into trashdir.

mkdir -p trashdir

for item in *; do
    # Skip trashdir itself
    if [ "$item" == "trashdir" ]; then
        continue
    fi
    # For directories: keep those starting with "run"
    if [ -d "$item" ]; then
        if [[ "$item" == run* ]]; then
            continue
        else
            mv "$item" trashdir/
        fi
    else
        # For files, preserve *.atm, *.sh, *.job, and the new lcao.in.
        case "$item" in
            *.atm|*.sh|*.job|lcao.in)
                continue
                ;;
            *)
                mv "$item" trashdir/
                ;;
        esac
    fi
done
