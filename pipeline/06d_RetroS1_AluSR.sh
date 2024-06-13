#!/bin/bash

usage="$(basename "$0") [-h] [-o i m r g s] -- Predict the confidence of each supporting read to be real MEI using RF, NB and LR

where:
    -h  show this help text
    -o  output folder path
    -m  masterpath
    -r  RetroSom version control"

while getopts ":ho:m:r:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
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
Rscript $masterpath/ALU/04_SR_level0/01.RF.r $outpath $ver $masterpath
$masterpath/ALU/04_SR_level0/02_combine_all.sh $outpath $ver $masterpath
