#!/bin/bash
#sbatch -N1 -c2 -t20-22 -o /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/26_plus_genomes/log/out.log -e /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/26_plus_genomes/log/er.log --mem=128gb -p celltypes ~/Matthew/code/SpeciesReferenceProcessing/annotations/cactus_26p.sh
#to setup followed install instructions, created XDG_RUNTIME_DIR, provided slurm config

eval "$(conda shell.bash hook)"
conda activate cactus
export PYTHONPATH=/home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/lib:$PYTHONPATH
export TOIL_SLURM_ARGS="--nice=1000"
export PATH=$PATH:/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/ucsc_tools/hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64
export PATH=$PATH:/home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/bin
export XDG_RUNTIME_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/XDG_RUNTIME_DIR


out_dir=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/26_plus_genomes
cd $out_dir

cp /home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/src/cactus/cactus_progressive_config.xml ${out_dir}/config-slurm.xml
sed -i ${out_dir}/config-slurm.xml -e 's/blast chunkSize="30000000"/blast chunkSize="90000000"/g'
sed -i ${out_dir}/config-slurm.xml -e 's/dechunkBatchSize="1000"/dechunkBatchSize="200"/g'

singularity version
#singularity pull quay.io/comparative-genomics-toolkit/cactus:v2.6.13
#echo "pulled"
#singularity shell cactus:v2.6.13.sif

#cactus-prepare ${out_dir}/out \
#--outDir ${out_dir}/ \
#--outSeqFile ${out_dir}/26pg.txt \
#--outHal ${out_dir}/26pg.hal \
#--jobStore /home/matthew.schmitz/Matthew/js \
#--configFile ${out_dir}/config-slurm.xml \
###--batchSystem slurm --batchLogsDir /home/matthew.schmitz/Matthew/log

cactus ${out_dir}/out \
${out_dir}/26pg.txt \
26pg.hal \
--realTimeLogging \
--batchSystem slurm \
--configFile ${out_dir}/config-slurm.xml \
--batchLogsDir /home/matthew.schmitz/Matthew/logs \
--binariesMode local \
--consMemory 128gb \
--consCores 4 

