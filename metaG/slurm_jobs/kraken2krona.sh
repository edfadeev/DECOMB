#!/bin/bash
#
#SBATCH --job-name=Kraken2krona
#SBATCH --cpus-per-task=2
#SBATCH --mem=3GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load kraken/2.0.8-beta
module load krona/2.7

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#define kraken tools directory
KRAK_TOOLS_DIR=/proj/DECOMB/source/KrakenTools

#run the program
python $KRAK_TOOLS_DIR/kreport2krona.py -r ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.report -o ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.krona

ktImportText ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.krona -o ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.krona.html
