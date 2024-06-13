#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i g a] -- BWA re-alignment and post-processing from aligned BAM file
Requires six cores, and 4-10G memory per core. Requires BWA (tested for v0.7.12), java, samtools and p
icard.

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID
    -g  path to reference genome hg38 version GRCh38DH
    -a  source BAM file

Please note that post-processing is slow, and it may take ~100hours to align sequencing reads from 30X W
GS data"

while getopts ":ho:i:g:a:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) sub="$OPTARG"
       ;;
    g) refG="$OPTARG"
       ;;
    a) bam="$OPTARG"
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

Itag="@RG\tID:$sub\tSM:sample1\tLB:library1\tPL:Illumina"
mkdir $outpath/$sub/align

### Alignment and post processing ###
### using BWA version 0.7.xx ###
### using human reference genome GRCh38DH ###
samtools sort -@ 6 -T $sub -n \
   $bam | \
   samtools fastq - | \
   bwa mem -t 6 -B 4 -O 6 -E 1 -M -p \
   -R $Itag \
   $refG/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   - \
   | /tools/bwakit/0.7.15/k8 /tools/bwakit/0.7.15/bwa-postalt.js -p \
   $refG/GRCh38_full_analysis_set_plus_decoy_hla-extra.fa \
   $refG/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt \
   | samtools view -bSh - \
   > $outpath/$sub/align/$sub.bam

### sort the aligned file by coordinates ###
samtools sort -@ 6 -T $sub \
   $outpath/$sub/align/$sub.bam \
   -o $outpath/$sub/align/$sub.sort.bam

