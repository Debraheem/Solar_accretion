#!/bin/bash

export OMP_NUM_THREADS=10

Y=0.26
Z=0.018
ALPHA=2.22
LOGS="LOGS"
FAST=0
FUTURE=0
MAX_AGE=4.568

while [ "$#" -gt 0 ]; do
  case "$1" in
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) ALPHA="$2"; shift 2;;
    -f) FAST=1; shift 1;;
    -F) FUTURE=1; shift 1;;
    -t) MAX_AGE="$2"; shift 2;;
    *) echo "unknown option: $1" >&2; exit 1;;
  esac
done

for INLIST in inlist_calibrate; do 
    shmesa change "$INLIST" mixing_length_alpha $ALPHA
    shmesa change "$INLIST" new_Y $Y
    shmesa change "$INLIST" new_Z $Z
    cp "$INLIST" "$LOGS"
    #./star inlist_calibrate
    ./rn
done 
