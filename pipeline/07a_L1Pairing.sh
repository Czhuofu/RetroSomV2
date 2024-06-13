#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i m r s n] -- Pairing the supporting reads for the same MEI

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID of the donor
    -m  masterpath
    -r  RetroSom version control
    -s  strandness (0, -strand; 1, +strand)
    -n  maximum number of supporting reads to be considered as a putative soamtic insertion (default 100)"

maxreads=100
while getopts ":ho:i:m:r:s:n:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    s) strand="$OPTARG"
       ;;
    n) maxreads="$OPTARG"
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   esac
done
$masterpath/LINE/05_Pairing/01_Stacking.pl $outpath $sub $ver $strand $maxreads
$masterpath/LINE/05_Pairing/02_pairs.pl $outpath $sub $ver $strand
