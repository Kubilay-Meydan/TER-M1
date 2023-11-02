#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y

WORKFLOW="workflows/${1}.smk"

if [ "$2" != "" ]; then
    RUN_NAME="$1"_"$2"
else
    RUN_NAME=$1
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}
cp config/test_pe_config.yaml ${FOLDER}/${RUN_NAME}_config.yaml


snakemake -s ${WORKFLOW} \
--use-conda \
--configfile config/test_pe_config.yaml \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100
