#!/bin/bash
#SBATCH --job-name=geom_xyz
#SBATCH --output=geom_xyz.out
#SBATCH --error=geom_xyz.err
#SBATCH --time=00:05:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks=10

# Luis Lorenzana
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025

# Refactored by Roman De Santos (rdesantos@hmc.edu)
# Harvey Mudd College
# Sandia National Laboratories Clinic 2025-2026

# how to run: sbatch convert_geom2xyz.sh input.geom [output.xyz]
# if no output file is specified, the script will create one using the input filename with a .xyz extension

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 input.xyz [output.geom]"
    exit 1
fi

INPUT="$1"
if [ "$#" -eq 2 ]; then
    OUTPUT="$2"
else
    OUTPUT="${INPUT%.*}.geom"
fi

FACTOR=1.8897259886  # angstrom → bohr

{
    # Optional header for .geom format
    echo "atom, type, position; step#     0"

    # Process .xyz file
    awk -v factor="$FACTOR" '
        BEGIN { OFMT="%.10f"; n=1 }

        # Skip the first two lines (atom count and comment)
        NR <= 2 { next }

        # Each remaining line: Symbol X Y Z
        NF >= 4 {
            symbol = $1
            x = $2 * factor
            y = $3 * factor
            z = $4 * factor
            printf "AT%-3d %-3s %14.10f %14.10f %14.10f\n", n, symbol, x, y, z
            n++
        }
    ' "$INPUT"
} > "$OUTPUT"

echo "✅ Conversion complete. Output written to $OUTPUT"
