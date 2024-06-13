#! /bin/bash

### rebuilt after SCG's update on singularity ###
### add -B /local/scratch/xwzhu/ for exonerate and Sort, as the temporary folder ###

usage="$(basename "$0") [-h] [-o i m r g t a b c n p e] -- Discovering somatic MEI insertions supported with >=2 reads.

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
    -e  control/normal tissue"

ver=1
hg=hg38
maxread=100
pcutoff=0.1
while getopts ":ho:i:m:r:t:g:a:b:c:n:p:e:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath=$(readlink -f $OPTARG)
       ;;
    i) sub="$OPTARG"
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
slurm_sc="-o %x.%A.output -e %x.%A.output -A aeurban -p batch --mem=50gb --time=100:00:00"
slurm_mc="-o %x.%A.output -e %x.%A.output -A aeurban -p batch --time=160:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=25gb"
tmppath=/local/scratch/xwzhu/
if [ "$datatype" == 1 ]
then
    jid1a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$read1,$read2 --app 01a_align_reads $masterpath/pipeline/RetroSomV2.6.sif -1 $read1 -2 $read2 -o $outpath -i $sub" | sbatch -J AlignReads $slurm_mc | awk '{print $4}')
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01c_clean_alignment $masterpath/pipeline/RetroSomV2.6.sif -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J CleanBam $slurm_sc --dependency=afterok:$jid1a | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01d_qc_alignment $masterpath/pipeline/RetroSomV2.6.sif -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J QCBam $slurm_sc --dependency=afterok:$jid1c | awk '{print $4}')
elif [ "$datatype" == 2 ]
then
    jid1b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01b_realign_BAM $masterpath/pipeline/RetroSomV2.6.sif -a $bam -o $outpath -i $sub" | sbatch -J RealignBAM $slurm_sc | awk '{print $4}')
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01c_clean_alignment $masterpath/pipeline/RetroSomV2.6.sif -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J CleanBAM $slurm_sc --dependency=afterok:$jid1b | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01d_qc_alignment $masterpath/pipeline/RetroSomV2.6.sif -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J QCBams $slurm_sc --dependency=afterok:$jid1c | awk '{print $4}')
elif [ "$datatype" == 3 ]
then
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01c_clean_alignment $masterpath/pipeline/RetroSomV2.6.sif -a $bam -o $outpath -i $sub" | sbatch -J CleanBam $slurm_sc | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01d_qc_alignment $masterpath/pipeline/RetroSomV2.6.sif -a $bam -o $outpath -i $sub" | sbatch -J QCBam $slurm_sc --dependency=afterok:$jid1c | awk '{print $4}')
elif [ "$datatype" == 4 ]
then
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01e_linkBAM $masterpath/pipeline/RetroSomV2.6.sif -a $bam -o $outpath -i $sub" | sbatch -J CleanBam $slurm_sc | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01d_qc_alignment $masterpath/pipeline/RetroSomV2.6.sif -a $bam -o $outpath -i $sub" | sbatch -J QCBam $slurm_sc --dependency=afterok:$jid1c | awk '{print $4}')
fi

##################################################
### Step2: Discover candidate supporting reads ###
##################################################
retro=retro_v$ver
if [ "$datatype" == 0 ]
then
jid2=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 02_RetroSom $masterpath/pipeline/RetroSomV2.6.sif $sub $outpath $retro $masterpath" | sbatch -J RetroDiscover $slurm_mc | awk '{print $4}')
else
jid2=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 02_RetroSom $masterpath/pipeline/RetroSomV2.6.sif $sub $outpath $retro $masterpath" | sbatch -J RetroDiscover $slurm_mc --dependency=afterok:$jid1c | awk '{print $4}')
fi

##############################################
### Step3: Putative MEIs without filtering ###
##############################################
jid3a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 03a_preprocess $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg" | sbatch -J PreProc $slurm_sc --dependency=afterok:$jid2 | awk '{print $4}')
# insertions in -strand #
jid3b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 03b_call_putative_MEI $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 -f 1" | sbatch -J PutMEI $slurm_sc --dependency=afterok:$jid2:$jid3a | awk '{print $4}')
# insertions in +strand #
jid3c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 03b_call_putative_MEI $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 -f 1" | sbatch -J PutMEI $slurm_sc --dependency=afterok:$jid2:$jid3a | awk '{print $4}')
######################################################
### Step4: Remapping L1HS or AluY specific Alleles ###
######################################################
# Realign L1 supporting reads #
jid4a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 04a_remap_L1 $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg" | sbatch -J L1MAP $slurm_sc --dependency=afterok:$jid3a | awk '{print $4}')
# Realign Alu supporting reads #
jid4b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 04b_remap_Alu $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg" | sbatch -J AluMAP $slurm_sc --dependency=afterok:$jid3a | awk '{print $4}')

####################################################
### Step5: Matricies for L1/Alu supporting reads ###
####################################################
# -strand #
jid5a0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05a_matrix_gen_L1PE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J L1PEMat $slurm_sc --dependency=afterok:$jid3b:$jid4a | awk '{print $4}')
jid5b0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05b_matrix_gen_L1SR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J L1SRMat $slurm_sc --dependency=afterok:$jid3b:$jid4a | awk '{print $4}')
jid5c0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05c_matrix_gen_AluPE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J AluPEMat $slurm_sc --dependency=afterok:$jid3b:$jid4b | awk '{print $4}')
jid5d0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05d_matrix_gen_AluSR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J AluSRMat $slurm_sc --dependency=afterok:$jid3b:$jid4b | awk '{print $4}')
# +strand #
jid5a1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05a_matrix_gen_L1PE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J L1PEMat $slurm_sc --dependency=afterok:$jid3c:$jid4a | awk '{print $4}')
jid5b1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05b_matrix_gen_L1SR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J L1SRMat $slurm_sc --dependency=afterok:$jid3c:$jid4a | awk '{print $4}')
jid5c1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05c_matrix_gen_AluPE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J AluPEMat $slurm_sc --dependency=afterok:$jid3c:$jid4b | awk '{print $4}')
jid5d1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05d_matrix_gen_AluSR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J AluSRMat $slurm_sc --dependency=afterok:$jid3c:$jid4b | awk '{print $4}')

###################################################
### Step6: Level1 prediction with RF, NB and LR ###
###################################################
jid6a0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06a_RetroS1_L1PE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J L1PEl1 $slurm_sc --dependency=afterok:$jid5a0:$jid5a1 | awk '{print $4}')
jid6b0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06b_RetroS1_L1SR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J L1SRl1 $slurm_sc --dependency=afterok:$jid5b0:$jid5b1 | awk '{print $4}')
jid6c0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06c_RetroS1_AluPE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J AluPEl1 $slurm_sc --dependency=afterok:$jid5c0:$jid5c1 | awk '{print $4}')
jid6d0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06d_RetroS1_AluSR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J AluSRl1 $slurm_sc --dependency=afterok:$jid5d0:$jid5d1 | awk '{print $4}')

jid6a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06a_Level1_L1PE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J L1PEl1 $slurm_sc --dependency=afterok:$jid5a0:$jid5a1 | awk '{print $4}')
jid6b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06b_Level1_L1SR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J L1SRl1 $slurm_sc --dependency=afterok:$jid5b0:$jid5b1 | awk '{print $4}')
jid6c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06c_Level1_AluPE $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J AluPEl1 $slurm_sc --dependency=afterok:$jid5c0:$jid5c1 | awk '{print $4}')
jid6d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06d_Level1_AluSR $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J AluSRl1 $slurm_sc --dependency=afterok:$jid5d0:$jid5d1 | awk '{print $4}')

###########################################
### Step7: Pairing the supporting reads ###
###########################################
#jid7a0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 07a_L1Pairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" | sbatch -J L1Pair $slurm_sc --dependency=afterok:$jid1d:$jid6a:$jid6b | awk '{print $4}')
#jid7a1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 07a_L1Pairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" | sbatch -J L1Pair $slurm_sc --dependency=afterok:$jid1d:$jid6a:$jid6b | awk '{print $4}')
#jid7b0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 07b_AluPairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 0 -n $maxreads" | sbatch -J AluPair $slurm_sc --dependency=afterok:$jid1d:$jid6c:$jid6d | awk '{print $4}')
#jid7b1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 07b_AluPairing $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -s 1 -n $maxreads" | sbatch -J AluPair $slurm_sc --dependency=afterok:$jid1d:$jid6c:$jid6d | awk '{print $4}')

########################################
### Step8: Prediction for read pairs ###
########################################
#jid8a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 08a_L1Predict $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont" | sbatch -J L1Predict  $slurm_sc --dependency=afterok:$jid7a0:$jid7a1 | awk '{print $4}')
#jid8b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 08b_AluPredict $masterpath/pipeline/RetroSomV2.6.sif -o $outpath -i $sub -m $masterpath -r $ver -g $hg -p $pcutoff -e $cont" | sbatch -J AluPredict  $slurm_sc --dependency=afterok:$jid7b0:$jid7b1 | awk '{print $4}')
