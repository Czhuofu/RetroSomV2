#!/bin/bash

# Original script: Xiaowei Zhu
# Modified: Junseok Park

# Changed the ln -s to ln -s (soft link)
# Use array and loop

workpath=$1
oldname=$2
newname=$3


if [ ! -d "$workpath" ]; then
    echo "$workpath does not exist"
    exit 0 
fi

cd $workpath

strands=("0" "1")
TEs=("LINE" "ALU")

if [ ! -d "./$oldname" ]; then
    echo "./$oldname does not exist"
    exit 0 
fi

if [ -d "./$newname" ]; then
    echo "./$newname does exist"
    #exit 0 
else
    ln -s ./$oldname ./$newname
fi


if [ -e "./$newname/insert.$oldname.txt" ]; then

    ln -srf ./$newname/insert.$oldname.txt ./$newname/insert.$newname.txt

else 
    echo " ./$newname/insert.$oldname.txt is not existing"
fi


for strand in "${strands[@]}"; do
    if [ -e "./$newname/retro_v1_$strand/$oldname.sr.tabe.discover" ]; then
        ln -srf ./$newname/retro_v1_$strand/$oldname.sr.tabe.discover ./$newname/retro_v1_$strand/$newname.sr.tabe.discover
    fi

    for TE in "${TEs[@]}"; do
        if [ -e "./$newname/retro_v1_$strand/$TE/$oldname.$TE.nodup.sites" ]; then
            ln -srf ./$newname/retro_v1_$strand/$TE/$oldname.$TE.nodup.sites ./$newname/retro_v1_$strand/$TE/$newname.$TE.nodup.sites
        fi
        if [ -e "./$newname/retro_v1_$strand/$TE/$oldname.pe.prob3.txt" ]; then
            ln -srf ./$newname/retro_v1_$strand/$TE/$oldname.pe.prob3.txt ./$newname/retro_v1_$strand/$TE/$newname.pe.prob3.txt
        fi
        if [ -e "./$newname/retro_v1_$strand/$TE/$oldname.sr.prob3.txt" ]; then
            ln -srf ./$newname/retro_v1_$strand/$TE/$oldname.sr.prob3.txt ./$newname/retro_v1_$strand/$TE/$newname.sr.prob3.txt
        fi
        if [ -e "./$newname/retro_v1_$strand/$oldname.pe.$TE.matrix" ]; then
            ln -srf ./$newname/retro_v1_$strand/$oldname.pe.$TE.matrix ./$newname/retro_v1_$strand/$newname.pe.$TE.matrix
        fi
        if [ -e "./$newname/retro_v1_$strand/$oldname.sr.$TE.matrix" ]; then
            ln -srf ./$newname/retro_v1_$strand/$oldname.sr.$TE.matrix ./$newname/retro_v1_$strand/$newname.sr.$TE.matrix
        fi
    done
done