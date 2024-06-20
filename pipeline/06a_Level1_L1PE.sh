#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i m r] -- Predict the confidence of each supporting read to be real MEI using RF, NB and LR

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID for the donor
    -m  masterpath
    -r  RetroSom version control"

while getopts ":ho:i:m:r:" opt; do
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

set -x

source $masterpath/pipeline/remove_R_lib_path.sh
remove_home_R_library_path

Rscript $masterpath/LINE/02_PE_level1/01.RF.LogR.NB.r $outpath $ver 1 $sub $masterpath
$masterpath/LINE/02_PE_level1/02_combine_RF.pl $outpath $ver 1 $sub
$masterpath/LINE/02_PE_level1/03_3probs.pl $outpath $ver 1 $sub

Rscript $masterpath/LINE/02_PE_level1/01.RF.LogR.NB.r $outpath $ver 0 $sub $masterpath
$masterpath/LINE/02_PE_level1/02_combine_RF.pl $outpath $ver 0 $sub
$masterpath/LINE/02_PE_level1/03_3probs.pl $outpath $ver 0 $sub
