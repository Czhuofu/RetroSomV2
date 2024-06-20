#!/bin/bash

# check wether you are in interactive session or not
HOST_NAME=$(hostname)

if [[ $HOST_NAME != "compute-"* ]]; then
    echo "run interactive session"
    srun --pty -p interactive -t 0-12:0:0 --mem 16GB -n 4 /bin/bash
fi

BASE_DIRECTORY=/n/data1/bch/genetics/lee
REF_DIRECTORY=${BASE_DIRECTORY}/projects/SMaHT/refs
WORK_DIRECTORY=${BASE_DIRECTORY}/projects/SMaHT/results/RetroSom
CONFIG_DIRECTORY=${BASE_DIRECTORY}/projects/SMaHT/benchmark/toolTest_pipeline/config

## omitted --fakeroot
singularity shell \
    -B /n/data1/bch/genetics/lee/:/n/data1/bch/genetics/lee/ \
    -B /n/scratch/users/j/jp394/:/n/scratch/users/j/jp394/ \
    -B /n/no_backup2/bch/lee/data/:/n/no_backup2/bch/lee/data/ \
    /n/app/singularity/containers/jp394/RetroSomV2.5.sif