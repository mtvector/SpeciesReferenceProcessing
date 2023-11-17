#!/bin/bash
#sbatch -N1 -c2 -t2-5 -o /home/matthew.schmitz/Matthew/genome/hal/test/log/out.log -e /home/matthew.schmitz/Matthew/genome/hal/test/log/er.log --mem=128gb -p celltypes ~/Matthew/code/SpeciesReferenceProcessing/annotations/cactus_test.sh
#to setup followed install instructions, created XDG_RUNTIME_DIR, provided slurm config

eval "$(conda shell.bash hook)"
conda activate cactus
export PYTHONPATH=/home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/lib:$PYTHONPATH
export TOIL_SLURM_ARGS="--nice=1000"
export PATH=$PATH:/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/ucsc_tools/hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64
export PATH=$PATH:/home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/bin
export XDG_RUNTIME_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/XDG_RUNTIME_DIR


out_dir=/home/matthew.schmitz/Matthew/genome/hal/test
cd $out_dir

cp /home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/src/cactus/cactus_progressive_config.xml ${out_dir}/config-slurm.xml
sed -i ${out_dir}/config-slurm.xml -e 's/blast chunkSize="30000000"/blast chunkSize="90000000"/g'
sed -i ${out_dir}/config-slurm.xml -e 's/dechunkBatchSize="1000"/dechunkBatchSize="200"/g'

cp /home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/examples/evolverMammals.txt ${out_dir}/

singularity version
#singularity pull quay.io/comparative-genomics-toolkit/cactus:v2.6.13
#echo "pulled"
#singularity shell cactus:v2.6.13.sif

#cactus-prepare /home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/examples/evolverMammals.txt \
#--outDir ${out_dir}/ \
#--outSeqFile ${out_dir}/evolverMammals.txt \
#--outHal ${out_dir}/evolverMammals.hal \
#--jobStore /home/matthew.schmitz/Matthew/js \
#--configFile ${out_dir}/config-slurm.xml \
###--batchSystem slurm --batchLogsDir /home/matthew.schmitz/Matthew/log

cactus ${out_dir}/out \
/home/matthew.schmitz/Matthew/utils/cactus-bin-v2.6.13/examples/evolverMammals.txt \
evolverMammals.hal \
--realTimeLogging \
--batchSystem slurm \
--configFile ${out_dir}/config-slurm.xml \
--batchLogsDir /home/matthew.schmitz/Matthew/logs \
--binariesMode local \
--consMemory 128gb \
--consCores 4 

