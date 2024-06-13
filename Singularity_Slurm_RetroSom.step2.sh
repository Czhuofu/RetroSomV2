#! /bin/bash

usage="$(basename "$0") [-h] [-o i m r g t a b c n p e l] -- Discovering somatic MEI insertions supported with >=2 reads.

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -m  masterpath
    -r  version control for RetroSom (default 1)
    -g  reference genome (default hg38, supporting hg38, hg19 and b37)
    -t  type of input (1=raw_sequencing_reads;2=BAM_to_be_realigned; 3=BAM_to_be_cleaned; 4=cleaned_BAM)
    -a  seqeuncing read1 (required if input are raw sequencing reads, input==1)
    -b  sequencing read2 (required if input are raw sequencing reads, input==1)
    -c  input BAM file (input==2 or 3)
    -n  maximum number of supporting reads to be considered as a putative soamtic insertion (default 100)
    -p  p_value cutoff (default p<0.1)
    -e  control/normal tissue
    -l  number of sequencing datasets"

ver=1
hg=hg38
maxread=100
pcutoff=0.1
while getopts ":ho:i:m:r:t:g:a:b:c:n:p:e:l:" opt; do
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
    t) datatype="$OPTARG"
       ;;
    g) hg="$OPTARG"
       ;;
    a) read1=$(readlink -f $OPTARG)
       ;;
    b) read2=$(readlink -f $OPTARG)
       ;;
    c) bam=$(readlink -f $OPTARG)
       ;;
    n) maxreads="$OPTARG"
       ;;
    p) pcutoff="$OPTARG"
       ;;
    e) cont="$OPTARG"
       ;;
    l) libs="$OPTARG"
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
mkdir $outpath/$sub
mkdir $outpath/$sub/reads
mkdir $outpath/$sub/align
mkdir $outpath/$sub/QC
mkdir $outpath/$sub/temp
mkdir $outpath/$sub/script
mkdir $outpath/$sub/retro_v$ver
cd $outpath/$sub/script
cp $masterpath/pipeline/*sh $outpath/$sub/script

### updating the location of the reference sequences ###
echo -e "Alu\t$masterpath/refTE/sequence/ALU.fa" > $masterpath/refTE/TE_ALHS.bed
echo -e "L1\t$masterpath/refTE/sequence/L1HS.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "HERVK\t$masterpath/refTE/sequence/HERVK.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "HERVH\t$masterpath/refTE/sequence/HERVH.fa" >> $masterpath/refTE/TE_ALHS.bed
echo -e "SVA\t$masterpath/refTE/sequence/SVA.fa" >> $masterpath/refTE/TE_ALHS.bed

#################################
### Step1: Sequence Alignment ###
################################# 
# Input Type 1: raw sequencing reads ### 
# Input Type 2: aligned reads to be realigned ### 
# Input Type 3: aligned reads (no realignment) ### 
slurm_sc="-o %x.%A.output -e %x.%A.output -A aeurban -p batch --mem=50gb --time=160:00:00"
slurm_mc="-o %x.%A.output -e %x.%A.output -A aeurban -p batch --time=160:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=25gb"

export TMPDIR=$outpath/$sub/script
###########################################
### Step7: Pairing the supporting reads ###
###########################################
jid7_0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07_combine.lib $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $tag -m $masterpath -r $ver -s 0 -l $libs -g $hg" | sbatch -J Combine $slurm_sc | awk '{print $4}')
jid7_1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07_combine.lib $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $tag -m $masterpath -r $ver -s 1 -l $libs -g $hg" | sbatch -J Combine $slurm_sc | awk '{print $4}')

jid7a0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07a_L1Pairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" | sbatch -J L1Pair $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')
jid7a1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07a_L1Pairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" | sbatch -J L1Pair $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')
jid7b0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07b_AluPairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" | sbatch -J AluPair $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')
jid7b1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 07b_AluPairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" | sbatch -J AluPair $slurm_sc --dependency=afterok:$jid7_0:$jid7_1 | awk '{print $4}')

########################################
### Step8: Prediction for read pairs ###
########################################
jid8a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 08a_L1Predict $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont" | sbatch -J L1Predict  $slurm_sc --dependency=afterok:$jid7a0:$jid7a1 | awk '{print $4}')
jid8b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$masterpath --app 08b_AluPredict $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont" | sbatch -J AluPredict  $slurm_sc --dependency=afterok:$jid7b0:$jid7b1 | awk '{print $4}')
