#! /bin/bash

### rebuilt after SCG's update on singularity ###
### add -B /local/scratch/xwzhu/ for exonerate and Sort, as the temporary folder ###

### Change script (07/24/2024, Junseok Park, junseok.park@childrens.harvard.edu)
### Running script is changed for O2 execution
### Change Singularity Image Path
## Original: /pipeline/RetroSomV2.5.sif
## O2:/n/app/singularity/containers/jp394/RetroSomV2.5.sif
## Original code: $sifimagepath 
## Changed code: $sifimagepath
## added option: s (singularity), f (temp file location)


usage="$(basename "$0") [-h] [-o i m r g t a b c n p e s f u l] -- Discovering somatic MEI insertions supported with >=2 reads.

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
    -s  singularity image file path including image file name
    -f  temporary directory for exonerate and Sort, as the temporary folder (default: /local/scratch/xwzhu/)
    -u  slurm short partition options (default : '-A aeurban -p batch --mem=50gb --time=100:00:00' )
    -l  slurm medium partition options (default: '-A aeurban -p batch --time=160:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=25gb' )"

# Initialize variables with default values
ver=1
hg=hg38
maxread=100
pcutoff=0.1
slurmshortags="\"-A aeurban -p batch --mem=50gb --time=100:00:00\""
slurmlongags="\"-A aeurban -p batch --time=160:00:00 --ntasks=1 --cpus-per-task=10 --mem-per-cpu=25gb\""
tmppath=/local/scratch/xwzhu/

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
    s) sifimagepath="$OPTARG"
       ;;
    f) tmppath="$OPTARG"
       ;;
    u) slurmshortags="$OPTARG"
       ;;
    l) slurmlongags="$OPTARG"
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
mkdir -p $outpath/$sub
mkdir -p $outpath/$sub/reads
mkdir -p $outpath/$sub/align
mkdir -p $outpath/$sub/QC
mkdir -p $outpath/$sub/temp
mkdir -p $outpath/$sub/script
mkdir -p $outpath/$sub/retro_v$ver
mkdir -p $outpath/$sub/logs
cd $outpath/$sub/script
cp $masterpath/pipeline/*sh $outpath/$sub/script
export TMPDIR=$outpath/$sub/script

# ### updating the location of the reference sequences ###
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
slurm_sc="-o $outpath/$sub/logs/%x.%A.output -e $outpath/$sub/logs/%x.%A.err $slurmshortags" 
slurm_mc="-o $outpath/$sub/logs/%x.%A.output -e $outpath/$sub/logs/%x.%A.err $slurmlongags" 

if [ "$datatype" == 1 ]
then
    jid1a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$read1,$read2 --app 01a_align_reads $sifimagepath -1 $read1 -2 $read2 -o $outpath -i $sub" | sbatch -J 1.1.1.AlignReads $slurm_lc | awk '{print $4}')
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01c_clean_alignment $sifimagepath -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J 1.1.2.CleanBam $slurm_lc --dependency=afterok:$jid1a | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01d_qc_alignment $sifimagepath -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J 1.1.3.QCBam $slurm_lc --dependency=afterok:$jid1c | awk '{print $4}')
elif [ "$datatype" == 2 ]
then
    jid1b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01b_realign_BAM $sifimagepath -a $bam -o $outpath -i $sub" | sbatch -J 1.2.1.RealignBAM $slurm_lc | awk '{print $4}')
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01c_clean_alignment $sifimagepath -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J 1.2.2.CleanBAM $slurm_lc --dependency=afterok:$jid1b | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 01d_qc_alignment $sifimagepath -a $outpath/$sub/align/$sub.sort.bam -o $outpath -i $sub" | sbatch -J 1.2.3.QCBams $slurm_lc --dependency=afterok:$jid1c | awk '{print $4}')
elif [ "$datatype" == 3 ]
then
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01c_clean_alignment $sifimagepath -a $bam -o $outpath -i $sub" | sbatch -J 1.3.1.CleanBam $slurm_lc | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01d_qc_alignment $sifimagepath -a $bam -o $outpath -i $sub" | sbatch -J 1.3.2.QCBam $slurm_lc --dependency=afterok:$jid1c | awk '{print $4}')
elif [ "$datatype" == 4 ]
then
    jid1c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01e_linkBAM $sifimagepath -a $bam -o $outpath -i $sub" | sbatch -J 1.4.1.CleanBam $slurm_lc | awk '{print $4}')
    jid1d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 01d_qc_alignment $sifimagepath -a $bam -o $outpath -i $sub" | sbatch -J 1.4.2.QCBam $slurm_lc --dependency=afterok:$jid1c | awk '{print $4}')
fi

##################################################
### Step2: Discover candidate supporting reads ###
##################################################
retro=retro_v$ver
if [ "$datatype" == 0 ]
then
   jid2=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 02_RetroSom $sifimagepath $sub $outpath $retro $masterpath" | sbatch -J 2.RetroDiscover $slurm_mc | awk '{print $4}')
else
   jid2=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 02_RetroSom $sifimagepath $sub $outpath $retro $masterpath" | sbatch -J 2.RetroDiscover $slurm_mc --dependency=afterok:$jid1c | awk '{print $4}')
fi

##############################################
### Step3: Putative MEIs without filtering ###
##############################################
jid3a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath,$bam --app 03a_preprocess $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg" | sbatch -J 3.1.PreProc $slurm_sc --dependency=afterok:$jid2 | awk '{print $4}')
# insertions in -strand #
jid3b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 03b_call_putative_MEI $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0 -f 1" | sbatch -J 3.2.PutMEI $slurm_sc --dependency=afterok:$jid2:$jid3a | awk '{print $4}')
# insertions in +strand #
jid3c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 03b_call_putative_MEI $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1 -f 1" | sbatch -J 3.3.PutMEI $slurm_sc --dependency=afterok:$jid2:$jid3a | awk '{print $4}')
######################################################
### Step4: Remapping L1HS or AluY specific Alleles ###
######################################################
# Realign L1 supporting reads #
jid4a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 04a_remap_L1 $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg" | sbatch -J 4.L1MAP $slurm_sc --dependency=afterok:$jid3a | awk '{print $4}')
# Realign Alu supporting reads #
jid4b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 04b_remap_Alu $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg" | sbatch -J 4.AluMAP $slurm_sc --dependency=afterok:$jid3a | awk '{print $4}')

####################################################
### Step5: Matricies for L1/Alu supporting reads ###
####################################################
# -strand #
jid5a0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05a_matrix_gen_L1PE $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J 5.1.L1PEMat $slurm_sc --dependency=afterok:$jid3b:$jid4a | awk '{print $4}')
jid5b0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05b_matrix_gen_L1SR $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J 5.1.L1SRMat $slurm_sc --dependency=afterok:$jid3b:$jid4a | awk '{print $4}')
jid5c0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05c_matrix_gen_AluPE $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J 5.1.AluPEMat $slurm_sc --dependency=afterok:$jid3b:$jid4b | awk '{print $4}')
jid5d0=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05d_matrix_gen_AluSR $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 0" | sbatch -J 5.1.AluSRMat $slurm_sc --dependency=afterok:$jid3b:$jid4b | awk '{print $4}')
# +strand #
jid5a1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05a_matrix_gen_L1PE $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J 5.2.L1PEMat $slurm_sc --dependency=afterok:$jid3c:$jid4a | awk '{print $4}')
jid5b1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05b_matrix_gen_L1SR $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J 5.2.L1SRMat $slurm_sc --dependency=afterok:$jid3c:$jid4a | awk '{print $4}')
jid5c1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05c_matrix_gen_AluPE $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J 5.2.AluPEMat $slurm_sc --dependency=afterok:$jid3c:$jid4b | awk '{print $4}')
jid5d1=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 05d_matrix_gen_AluSR $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver -g $hg -s 1" | sbatch -J 5.2.AluSRMat $slurm_sc --dependency=afterok:$jid3c:$jid4b | awk '{print $4}')


###################################################
### Step6: Level1 prediction with RF, NB and LR ###
###################################################
jid6a=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06a_Level1_L1PE $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J 6.2.L1PEl1 $slurm_sc --dependency=afterok:$jid5a0:$jid5a1 | awk '{print $4}')
jid6b=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06b_Level1_L1SR $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J 6.2.L1SRl1 $slurm_sc --dependency=afterok:$jid5b0:$jid5b1  | awk '{print $4}')
jid6c=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06c_Level1_AluPE $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J 6.2.AluPEl1 $slurm_sc --dependency=afterok:$jid5c0:$jid5c1 | awk '{print $4}')
jid6d=$(echo -e '#!/bin/sh\n'"singularity run -B $outpath,$tmppath,$masterpath --app 06d_Level1_AluSR $sifimagepath -o $outpath -i $sub -m $masterpath -r $ver" | sbatch -J 6.2.AluSRl1 $slurm_sc --dependency=afterok:$jid5d0:$jid5d1  | awk '{print $4}')


