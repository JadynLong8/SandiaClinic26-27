#!/bin/bash
#SBATCH --job-name=postprocess
#SBATCH --output=postprocess.out
#SBATCH --error=postprocess.err
#SBATCH --time=00:05:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks=1

# Roman De Santos (HMC 26') [Refactored for NEB process and readability]
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025-26
# Email: rdesantos@hmc.edu

# Postprocess script:
# Inputs: lcao.in and lcao.neb_geom
# Outputs:
# Top Level Directory: 
# - lcao.in (with updated coordinate positions)
# - Will keep *.atm, *.sh, *.job, everything else will go in trashdir 
# - runN Directory: Everytime you run postprocess.sh; N+1
# -- lcao.in (Original lcao.in file)
# -- lcao.out (Output file for N run)
# -- lcao.hist (History file with coordinates and forces at every G step)
# -- lcao.geom (Last updated coordinate postion of structure)
# -- lcao.xyz (Visualization file to see strecture in VESTA or Jmol or IQmol)


# --- 1. IDENTIFY RUN DIRECTORY ---
max_val=-1
for d in run*/ ; do
    if [[ $d =~ run([0-9]+)/ ]]; then
        num="${BASH_REMATCH[1]}"
        if (( num > max_val )); then
            max_val=$num
        fi
    fi
done

next_val=$((max_val + 1))
RUN_DIR="run${next_val}"
BASE_PATH=$(pwd)
FULL_RUN_PATH="${BASE_PATH}/${RUN_DIR}"

echo "Creating archive: $RUN_DIR"
mkdir -p "$RUN_DIR"
mkdir -p "trashdir"

# --- 2. MOVE AND PRESERVE RESULTS ---
# We move the current files into the archive folder first.
# This preserves the "old" lcao.in inside the runN directory.
mv lcao.out lcao.neb_geom lcao.neb_post lcao.stat lcao.hist lcao.in lcao.geom[0-9]* "$RUN_DIR" 2>/dev/null

# Clean up other auxiliary files
shopt -s extglob
mv !(*.sh|*.job|*.jobfile|*.atm|run*|trashdir) trashdir/ 2>/dev/null

# --- 3. POST-PROCESSING (Inside Run Directory) ---
echo "Processing XYZ images in $RUN_DIR..."
module purge
module load anaconda3/2022.10
module load ffmpeg

cd "$FULL_RUN_PATH"
python3 /ocean/projects/mat250017p/desantos/scripts/NEBimages2xyz.py "${FULL_RUN_PATH}/lcao.neb_geom" "out"

# --- 4. GENERATE NEW LCAO.IN AT BASE_PATH ---
echo "Splicing coordinates for the next iteration..."

# We reference the files we just moved into the run directory
line_in=$(grep -n "image number in NEB" lcao.in | head -n 1 | cut -d: -f1)
line_geom=$(grep -n "image number in NEB" lcao.neb_geom | head -n 1 | cut -d: -f1)

if [[ -n "$line_in" && -n "$line_geom" ]]; then
    # Create the NEW file back at the project root (BASE_PATH)
    # 1. Take the parameters/header from the archived lcao.in
    head -n $((line_in - 1)) lcao.in > "${BASE_PATH}/lcao.in"
    # 2. Append the updated coordinates from the archived lcao.neb_geom
    tail -n +$line_geom lcao.neb_geom >> "${BASE_PATH}/lcao.in"
    
    echo "SUCCESS: New lcao.in generated at ${BASE_PATH}"
    echo "REFERENCE: Original input preserved at ${FULL_RUN_PATH}/lcao.in"
else
    echo "ERROR: 'image number in NEB' keyword missing. Could not generate new lcao.in."
fi

cd "$BASE_PATH"
echo "Workflow complete. Ready for next submission."