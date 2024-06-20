#!/bin/bash -l

usage="$(basename "$0") [-h] [-o i m r s l g] -- Pairing the supporting reads for the same MEI

where:
    -h  show this help text
    -o  output folder path
    -i  subject ID of the donor
    -m  masterpath
    -r  RetroSom version control
    -s  strandness (0, -strand; 1, +strand)
    -l  number of the subset sequencing datasets
    -g  version of the human reference genome"

while getopts ":ho:i:m:r:s:l:g:" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    o) outpath="$OPTARG"
       ;;
    i) subject="$OPTARG"
       ;;
    m) masterpath="$OPTARG"
       ;;
    r) ver="$OPTARG"
       ;;
    s) strand="$OPTARG"
       ;;
    l) readgroup="$OPTARG"
       ;;
    g) hg="$OPTARG"
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

#echo $1  ### subject ID ###
#echo $2  ### path ###
#echo retro_v$ver\_$strand ### output folder ###
#echo $4  ### reference genome hg19 G37 hg38... ###
#echo $5 ##strandness ###
sub=$subject\_NoModel
sub2=$subject\-

TE=LINE
tmpfolder=$outpath/$sub

#mkdir $outpath/$sub
mkdir $outpath/$sub/retro_v$ver\_$strand
mkdir $outpath/$sub/retro_v$ver\_$strand/$TE

### combine 6 invidual libraries ###
### combine SR support reads ###
rm -f $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover


for (( c=1; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$5="lib"lib"_"$5' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.tabe.discover >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover
     done

### combine the PE reads ###
rm -f $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
for (( c=1; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$5="lib"lib"_"$5' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.$TE.nodup.sites >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
     done

### make calls ###
awk '{if ($4 ~ /L1/) print}' $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.SR.calls

### PE calling ###
### L1 PE calls ###
$masterpath/utls/06_ana_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE 1000
mv $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.calls $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.nodup.calls
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub retro_v$ver\_$strand $TE $masterpath

windowBed \
  -w 100 -v \
  -a $outpath/$sub/retro_v$ver\_$strand/LINE/$sub.LINE.SR.PE.calls \
  -b $masterpath/refTE/position/$hg\.fa_LINE1\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/retro_v$ver\_$strand/LINE/$sub.LINE.noref.calls

$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub retro_v$ver\_$strand LINE

### combine the prob3 ###
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$TE/$sub2\1.pe.prob3.txt > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.pe.prob3.txt
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$TE/$sub2\1.sr.prob3.txt > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.sr.prob3.txt
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.pe.LINE.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.pe.LINE.matrix
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.sr.LINE.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.sr.LINE.matrix
if [ "$strand" == '1' ]
   then
      rm -f $outpath/$sub/insert.$sub.txt
   fi

for (( c=1; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.pe.prob3.txt | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.pe.prob3.txt
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.sr.prob3.txt | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.sr.prob3.txt
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.pe.LINE.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.pe.LINE.matrix
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.LINE.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.LINE.matrix
       if [ "$strand" == '1' ]
         then
           awk -v lib=$c '$3="lib"lib' OFS='\t' $outpath/$sub2$c/insert.$sub2$c.txt >>  $outpath/$sub/insert.$sub.txt
         fi
     done

TE=ALU
mkdir $outpath/$sub/retro_v$ver\_$strand/$TE

### combine the PE reads ###
rm -f $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
for (( c=1; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$5="lib"lib"_"$5' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.$TE.nodup.sites >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.novel.sites
     done

### make calls ###
awk '{if ($4 ~ /Alu/) print}' $outpath/$sub/retro_v$ver\_$strand/$sub.sr.tabe.discover | grep 'OK' | \
   awk '{print $1"\t"$2"\t"$3"\t"$5}' | sort -T $tmpfolder -u | \
   sort -T $tmpfolder -k1,1 -k2,3n | \
   mergeBed -d 40 -c 4 -o distinct -delim ";" \
   -i stdin | \
   awk '{split ($4, num, ";"); print $1"\t"$2"\t"$3"\tsr,"length(num)","$4}' \
   > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.SR.calls

### PE calling ###
### L1 PE calls ###
$masterpath/utls/06_ana_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE 1000
mv $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.calls $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.$TE.PE.nodup.calls
$masterpath/utls/08_refine_depth.pl $sub $outpath/$sub retro_v$ver\_$strand $TE

### combine SR and PE ###
date '+%m/%d/%y %H:%M:%S'
echo "Calling PE phase ... Combine SR and PE"
$masterpath/utls/09_merge.SR.PE.support.sh $sub $outpath/$sub retro_v$ver\_$strand $TE $masterpath

windowBed \
  -w 100 -v \
  -a $outpath/$sub/retro_v$ver\_$strand/ALU/$sub.ALU.SR.PE.calls \
  -b $masterpath/refTE/position/$hg\.fa_ALU\_$strand.bed \
  | windowBed -w 10 -v \
      -a stdin \
      -b $masterpath/refTE/position/$hg.mask.bed \
      > $outpath/$sub/retro_v$ver\_$strand/ALU/$sub.ALU.noref.calls

$masterpath/utls/16_prefilter_r1r2.pl $sub $outpath/$sub retro_v$ver\_$strand ALU

### combine the prob3 ###
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$TE/$sub2\1.pe.prob3.txt > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.pe.prob3.txt
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$TE/$sub2\1.sr.prob3.txt > $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.sr.prob3.txt
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.pe.ALU.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.pe.ALU.matrix
head -n 1 $outpath/$sub2\1/retro_v$ver\_$strand/$sub2\1.sr.ALU.matrix > $outpath/$sub/retro_v$ver\_$strand/$sub.sr.ALU.matrix

for (( c=1; c<=$readgroup; c++ ))
     do
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.pe.prob3.txt | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.pe.prob3.txt
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$TE/$sub2$c.sr.prob3.txt | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$TE/$sub.sr.prob3.txt
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.pe.ALU.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.pe.ALU.matrix
       awk -v lib=$c '$4="lib"lib"_"$4' OFS='\t' $outpath/$sub2$c/retro_v$ver\_$strand/$sub2$c.sr.ALU.matrix | awk 'FNR > 1' >> $outpath/$sub/retro_v$ver\_$strand/$sub.sr.ALU.matrix
     done

