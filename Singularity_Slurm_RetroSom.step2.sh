#! /bin/bash

usage="$(basename "$0") [-h] [-o i m r g n p e l s u f] -- Discovering somatic MEI insertions supported with >=2 reads.

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath
    -r  version control for RetroSom (default 1)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)
    -n  maximum number of supporting reads to be considered as a putative soamtic insertion (default 100)
    -p  p_value cutoff (default p<0.1)
    -e  control/normal tissue
    -l  number of sequencing datasets
    -s  singularity image file path including image file name
    -u  slurm short partition options (default : '-A aeurban -p batch --mem=50gb --time=100:00:00' )
    -f  delete intermediate file (default=FALSE) "

ver=1
hg=hg38
maxread=100
pcutoff=0.1
slurmshortags="\"-A aeurban -p batch --mem=50gb --time=100:00:00\""
slurmlongags="\"-A aeurban -p batch --time=160:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=25gb\""
deleteIntermediateFiles=false

while getopts ":ho:i:m:r:g:n:p:e:l:s:u:f:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    i) tag="$OPTARG"
       ;;
    m) masterpath=$(readlink -f $OPTARG)
       ;;
    r) ver="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    n) maxreads="$OPTARG"
       ;;
    p) pcutoff="$OPTARG"
       ;;
    e) cont="$OPTARG"
       ;;
    l) libs="$OPTARG"
       ;;
    s) sifimagepath="$OPTARG"
       ;;
    u) slurmshortags="$OPTARG"
       ;;
    f) deleteIntermediateFiles="$OPTARG"
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

######################################
### Step0: Creating output folders ###
######################################
sub=$tag'_NoModel'
mkdir -p $outpath/$sub
mkdir -p $outpath/$sub/reads
mkdir -p $outpath/$sub/align
mkdir -p $outpath/$sub/QC
mkdir -p $outpath/$sub/temp
mkdir -p $outpath/$sub/script
mkdir -p $outpath/$sub/retro_v$ver
cd $outpath/$sub/script
cp $masterpath/pipeline/*sh $outpath/$sub/script
export TMPDIR=$outpath/$sub/script

slurm_sc="-o $outpath/$sub/logs/%x.%A.output -e $outpath/$sub/logs/%x.%A.err $slurmshortags" 

###########################################
### Step7: Pairing the supporting reads ###
###########################################
jid7_0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07_combine.lib $sifimagepath -o $outpath -i $tag -m $masterpath -r $ver -s 0 -l $libs -g $hg" | sbatch -J Combine0 $slurm_sc | awk '{print $4}')
jid7_1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07_combine.lib $sifimagepath -o $outpath -i $tag -m $masterpath -r $ver -s 1 -l $libs -g $hg" | sbatch -J Combine1 $slurm_sc | awk '{print $4}')

jid7a0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07a_L1Pairing $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" | sbatch -J L1Pair0 $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')
jid7a1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07a_L1Pairing $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" | sbatch -J L1Pair1 $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')
jid7b0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07b_AluPairing $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" | sbatch -J AluPair0 $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')
jid7b1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07b_AluPairing $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" | sbatch -J AluPair1 $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')

########################################
### Step8: Prediction for read pairs ###
########################################
jid8a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 08a_L1Predict $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont" | sbatch -J L1Predict  $slurm_sc --dependency=afterok:$jid7a0:$jid7a1 | awk '{print $4}')
jid8b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 08b_AluPredict $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont" | sbatch -J AluPredict  $slurm_sc --dependency=afterok:$jid7b0:$jid7b1 | awk '{print $4}')

########################################
### Step9: delete intermediate files ###
########################################
if deleteIntermediateFiles; then 
   jid9=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 09_delete_non_essentialFiles $sifimagepath -t $outpath" | sbatch -J DeleteFiles $slurm_sc --dependency=afterok:$jid8a:$jid8b | awk '{print $4}')
else
   echo "intermediate file will not be deleted"
fi