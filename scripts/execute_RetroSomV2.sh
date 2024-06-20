#!/bin/bash

LANG=en_US.UTF-8
LANGUAGE=en_US
LC_CTYPE="en_US.UTF-8"
LC_NUMERIC=ro_RO.UTF-8
LC_TIME=ro_RO.UTF-8
LC_COLLATE="en_US.UTF-8"
LC_MONETARY=ro_RO.UTF-8
LC_MESSAGES="en_US.UTF-8"
LC_PAPER=ro_RO.UTF-8
LC_NAME=ro_RO.UTF-8
LC_ADDRESS=ro_RO.UTF-8
LC_TELEPHONE=ro_RO.UTF-8
LC_MEASUREMENT=ro_RO.UTF-8
LC_IDENTIFICATION=ro_RO.UTF-8
LC_ALL=

R_VERSION=3.5.1

module load gcc/6.2.0 R/$R_VERSION

echo 'R_LIBS_USER="~/R-'${R_VERSION}'/library"' > $HOME/.Renviron
export R_LIBS_USER="~/R-"${R_VERSION}

module load bedtools/2.27.1 picard/2.27.5 samtools/1.15.1

RESULT_PATH=""

RETROPATH=""
SAMPLE_DIR=""
SINGULARITY_IMAGE=""
SLURM_PARTITION=medium
SAMPLE_PREFIX=.bam

# Create result directory if it is not existing
if [ ! -d "${RESULT_PATH}" ]; then
    mkdir -p ${RESULT_PATH}
fi

# Print the array elements
while IFS= read -r line; do

    IFS=' ' read -ra array <<< "$line"
    SAMPLE_ID=${array[0]}
    CONT_ID=${array[1]}

    echo "$SAMPLE_ID is loaded"
    # Create sample directory if it is not exists
    if [ -d "${RESULT_PATH}/${SAMPLE_ID}" ]; then
        rm -r ${RESULT_PATH}/${SAMPLE_ID}
        echo "Directory ${RESULT_PATH}/${SAMPLE_ID} has been removed."
    fi

    BAMFILE_PATH=$SAMPLE_DIR/${SAMPLE_ID}${SAMPLE_PREFIX}

    if [ -e "$BAMFILE_PATH" ]; then
        command="${RETROPATH}/Singularity_Slurm_RetroSomV2.5_test.sh -o $RESULT_PATH \
            -i $SAMPLE_ID \
            -e $CONT_ID \
            -m $RETROPATH \
            -r 1
            -g hg38 \
            -t 3 \
            -c $BAMFILE_PATH \
            -s $SINGULARITY_IMAGE 
            "

        echo $command

        eval $command
    else
        echo "$BAMFILE_PATH is not existing"
    fi

done < "sampleIds.txt"


