#!/bin/bash
#SBATCH --job-name=geom_xyz
#SBATCH --output=geom_xyz.out
#SBATCH --error=geom_xyz.err
#SBATCH --time=00:05:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks=1

# Luis Lorenzana
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025

# Refactored by Roman De Santos (rdesantos@hmc.edu)
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025-2026


# how to run: sbatch convert_geom2xyz.sh input.geom [output.xyz]
# if no output file is specified, the script will create one using the input filename with a .xyz extension

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 input.geom [output.xyz]"
    exit 1
fi

INPUT="$1"
if [ "$#" -eq 2 ]; then
    OUTPUT="$2"
else
    OUTPUT="${INPUT%.*}.xyz"
fi

FACTOR=0.529177  # conversion from bohr to angstrom
# count non-empty lines excluding header lines that contain "atom, type, position;"
NUM_ATOMS=$(grep -v 'atom, type, position;' "$INPUT" | grep -cv '^\s*$')

{
    # write the header for the .xyz file.
    echo "$NUM_ATOMS"
    echo ".geom to .xyz and scaled by factor of $FACTOR"

    # process the .geom file:
    # - Skip any line containing the header text
    # - If the line has 5 fields (index, symbol, x, y, z), drop the index
    # - Otherwise, assume the line contains at least 4 fields: symbol, x, y, z
    # - Multiply the coordinates by the factor
awk -v factor="$FACTOR" '
    /^[[:space:]]*$/ { next }                       # Skip empty lines
    /atom, type, position;/ { next }                # Skip header lines
    {
        # Find last 4 fields (symbol, x, y, z) even if prefix like AT1 exists
        n = NF
        symbol = $(n-3)
        x = $(n-2) * factor
        y = $(n-1) * factor
        z = $(n) * factor
        printf "%s %.10f %.10f %.10f\n", symbol, x, y, z
    }
' "$INPUT"

} > "$OUTPUT" 

echo "Conversion complete. Output written to $OUTPUT"