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
    -l  number of sequencing datasets
    -s  singularity image file path including image file name"

ver=1
hg=hg38
maxread=100
pcutoff=0.1
while getopts ":ho:i:m:r:t:g:a:b:c:n:p:e:l:s:" opt; do
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
    s) sifimagepath="$OPTARG"
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
### Step0.0: Renaming directories if it is not available ###
######################################



######################################
### Step0.1: Creating output folders ###
######################################
create_directories() {
   local outpath="$1"
   local sub="$2"
   local ver="$3"
   local masterpath="$4"

   # Create directories
   mkdir -p "$outpath/$sub/reads"
   mkdir -p "$outpath/$sub/align"
   mkdir -p "$outpath/$sub/QC"
   mkdir -p "$outpath/$sub/temp"
   mkdir -p "$outpath/$sub/logs"
   mkdir -p "$outpath/$sub/script"
   mkdir -p "$outpath/$sub/retro_v$ver"

   # Copy script files
   cd "$outpath/$sub/script"
   cp "$masterpath/pipeline/"*sh .
   
}

cont_sub=$cont'_NoModel'
create_directories $outpath $cont_sub $ver $masterpath

sub=$tag'_NoModel'
create_directories $outpath $sub $ver $masterpath

# sub=$tag'_NoModel'
# mkdir $outpath/$sub
# mkdir $outpath/$sub/reads
# mkdir $outpath/$sub/align
# mkdir $outpath/$sub/QC
# mkdir $outpath/$sub/temp
# mkdir $outpath/$sub/logs
# mkdir $outpath/$sub/script
# mkdir $outpath/$sub/retro_v$ver
# cd $outpath/$sub/script
# cp $masterpath/pipeline/*sh $outpath/$sub/script # (Temp comment)

# cont_sub=$cont'_NoModel'
# mkdir $outpath/$cont_dir
# mkdir $outpath/$cont_dir/reads
# mkdir $outpath/$cont_dir/align
# mkdir $outpath/$cont_dir/QC
# mkdir $outpath/$cont_dir/temp
# mkdir $outpath/$cont_dir/script
# mkdir $outpath/$cont_dir/retro_v$ver


## Change it to cotain all the samples as input

#echo $sub > $outpath/$sub/list.txt

### updating the location of the reference sequences ###
## Temp comment -S
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
# slurm_sc="-o %x.%A.output -e %x.%A.output -A aeurban -p batch --mem=50gb --time=160:00:00"
# slurm_mc="-o %x.%A.output -e %x.%A.output -A aeurban -p batch --time=160:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=25gb"

slurm_sc="-o $outpath/$sub/logs/%x.%A.output -e $outpath/$sub/logs/%x.%A.err -p medium --mem=50gb --time=3-23:00"
slurm_lc="-o $outpath/$sub/logs/%x.%A.output -e $outpath/$sub/logs/%x.%A.err -p long --mem=50gb --time=15-23:00"
slurm_mc="-o $outpath/$sub/logs/%x.%A.output -e $outpath/$sub/logs/%x.%A.err -p medium --time=3-23:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=25gb"
tmppath=/n/scratch/users/j/jp394


export TMPDIR=$outpath/$sub/script
###########################################
### Step7: Pairing the supporting reads ###
###########################################

array=($tag $cont)

for ((i=0; i<"${#array[@]}"; i++)); do
   echo "${array[i]}"
   for ((j=0; j<2; j++)); do
      singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07_combine.lib $sifimagepath  -o $outpath -i ${array[i]} -m $masterpath -r $ver -s $j -l $libs -g $hg     
   done
done 

# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07_combine-$tag.lib $sifimagepath  -o $outpath -i $tag -m $masterpath -r $ver -s 0 -l $libs -g $hg
# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07_combine-$tag.lib $sifimagepath  -o $outpath -i $tag -m $masterpath -r $ver -s 1 -l $libs -g $hg

# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07_combine-$cont.lib $sifimagepath  -o $outpath -i $cont -m $masterpath -r $ver -s 0 -l $libs -g $hg
# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07_combine-$cont.lib $sifimagepath  -o $outpath -i $cont -m $masterpath -r $ver -s 1 -l $libs -g $hg

array=($sub $cont_sub)
step7_apps=("07a_L1Pairing" "07b_AluPairing")
for ((i=0; i<"${#step7_apps[@]}"; i++)); do
   for ((j=0; j<"${#array[@]}"; j++)); do
      for ((k=0; k<2; k++)); do
         singularity run -B $outpath,$tmppath,$masterpath,$bam --app ${step7_apps[i]} $sifimagepath  -o $outpath -i ${array[i]} -m $masterpath -r $ver -s $k -n $maxreads
      done
   done
done


# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07a_L1Pairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads
# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07a_L1Pairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads

# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07b_AluPairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads
# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07b_AluPairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads


# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07a_L1Pairing $sifimagepath  -o $outpath -i $cont_sub -m $masterpath -r $ver -s 0 -n $maxreads
# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07a_L1Pairing $sifimagepath  -o $outpath -i $cont_sub -m $masterpath -r $ver -s 1 -n $maxreads

# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07b_AluPairing $sifimagepath  -o $outpath -i $cont_sub -m $masterpath -r $ver -s 0 -n $maxreads
# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07b_AluPairing $sifimagepath  -o $outpath -i $cont_sub -m $masterpath -r $ver -s 1 -n $maxreads


# jid7a0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07a_L1Pairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" )
# jid7a1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07a_L1Pairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" )
# jid7b0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07b_AluPairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" )
# jid7b1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 07b_AluPairing $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" )

# ########################################
# ### Step8: Prediction for read pairs ###
# ########################################

step8_apps=("08a_L1Predict" "08b_AluPredict")
for ((i=0; i<"${#step8_apps[@]}"; i++)); do
   singularity run -B $outpath,$tmppath,$masterpath,$bam --app ${step8_apps[i]} $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont_sub
done

# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 08a_L1Predict $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont_sub
# singularity run -B $outpath,$tmppath,$masterpath,$bam --app 08b_AluPredict $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont_sub

# jid8a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 08a_L1Predict $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont" )
# jid8b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 08b_AluPredict $sifimagepath  -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont")