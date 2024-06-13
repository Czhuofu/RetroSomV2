#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i m r g s] -- Create feature matrix for the L1 SR supporting reads. Requires one core, and 8gb memory
Requires bedtools, blast and exonerate

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath (default ~/masterpath)
    -r  RetroSom version control (default 1)
    -g  hg38 or hg19
    -s  strand (1, +strand; 0, -strand)"

ver=1
masterpath=~/masterpath
while getopts ":ho:i:m:r:g:s:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    s) strand="$OPTARG"
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

$masterpath/LINE/03_SR_matrix/01_L1SR_matrix.pl $outpath $sub $ver $strand $hg $masterpath
