#!/bin/bash -l

sub=$1
outpath=$2
Retro=$3
masterpath=$4

date '+%m/%d/%y %H:%M:%S'
echo "Discovery phase ... begins"
$masterpath/RetroSeq/retroseq.prun.pl -discover \
  -align \
  -srmode \
  -minclip 20 \
  -len 26 \
  -srcands $outpath/$sub/$Retro/$sub.sr.discover \
  -bam $outpath/$sub/align/$sub.final.bam \
  -eref $masterpath/refTE/TE_ALHS.bed \
  -output $outpath/$sub/$Retro/$sub.discover

### seperate the supporting reads from either strand ###
$masterpath/utls/21_par_direction.pl $sub $outpath $Retro
date '+%m/%d/%y %H:%M:%S'
echo "Discovery phase ... finish"

