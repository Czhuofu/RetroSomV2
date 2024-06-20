#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i m r s n g p e] -- Pairing the supporting reads for the same MEI

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID of the donor
    -m  masterpath
    -r  RetroSom version control
    -s  strandness (0, -strand; 1, +strand)
    -n  maximum number of supporting reads to be considered as a putative soamtic insertion (default 100)
    -g  reference genome (hg38, hg19, b37)
    -p  p_value cutoff (default p<0.1)
    -e  control/normal tissue"

maxreads=100
pcutoff=0.1
while getopts ":ho:i:m:r:s:n:g:p:e:" opt; do
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
    g) hg="$OPTARG"
       ;;
    p) pcutoff="$OPTARG"
       ;;
    e) cont="$OPTARG"
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
Rscript $masterpath/ALU/06_Predict/01.Model.Stacking.r $outpath $ver $sub $masterpath
$masterpath/ALU/06_Predict/02_post_filtering.pl $outpath $ver $sub $cont $masterpath $hg $pcutoff
$masterpath/ALU/06_Predict/04_plot_Alu.pl $outpath $sub $ver 0.1 $masterpath $hg
