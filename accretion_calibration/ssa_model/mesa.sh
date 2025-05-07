#!/bin/bash

# Step 1: Set threading
export OMP_NUM_THREADS=20

# Step 2: Default parameters
Y=0.2703           # initial helium mass fraction
Z=0.018            # initial metallicity
ALPHA=2.22         # mixing-length parameter
LOGS="LOGS"        # output directory
FAST=0             # fast-flag
FUTURE=0           # future-flag
MAX_AGE=4.572      # target age in Gyr

# Step 3: Parse command-line flags
while [ "$#" -gt 0 ]; do
  case "$1" in
    -Y) Y="$2"; shift 2 ;;                   # override helium
    -Z) Z="$2"; shift 2 ;;                   # override metallicity
    -a) ALPHA="$2"; shift 2 ;;               # override αMLT
    -f) FAST=1;     shift 1 ;;               # fast mode
    -F) FUTURE=1;   shift 1 ;;               # include future evolution
    -t) MAX_AGE="$2"; shift 2 ;;             # override max age
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Step 4: Ensure log directory exists
mkdir -p "$LOGS"

# Step 5: Loop over inlist_main and apply controls
for INLIST in inlist_main; do
  echo "=== Calibrating using $INLIST ==="

  # 5.1 Set mixing‐length parameter
  shmesa change "$INLIST" controls mixing_length_alpha "$ALPHA"

  # 5.2 Set initial helium via x_ctrl(7)
  shmesa change "$INLIST" controls "x_ctrl(7)"            "$Y"    # initial Helium Abundance

  # 5.3 Set initial metals
  shmesa change "$INLIST" controls initial_z            "$Z"

  # 5.4 Archive this inlist
  cp "$INLIST" "$LOGS/"

  # 5.5 Run MESA
  #./star "$INLIST"
  ./rn_accrete
done

