#!/bin/bash
#
#SBATCH --job-name=Kraken2_classificaion
#SBATCH --cpus-per-task=12
#SBATCH --mem=48GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load kraken/2.0.8-beta

#Set up the path to the working directory
WORKDIR=/proj/DECOMB/analysis/process_metaG/
cd $WORKDIR

#run the program
kraken2 --paired --use-names --db /mirror/kraken2/kraken_standard_db --threads 12 --output ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.output --report ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.report --classified-out ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_classified_#.fastq.gz ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R1.fastq.gz ./TRIM/CDT3KANXX_${SLURM_ARRAY_TASK_ID}_QC_R2.fastq.gz

#produce report in krona
module load krona/2.7

#define kraken tools directory
KRAK_TOOLS_DIR=/proj/DECOMB/source/KrakenTools

#run the program
python $KRAK_TOOLS_DIR/kreport2krona.py -r ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.report -o ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.krona

ktImportText ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.krona -o ./kraken2/CDT3KANXX_${SLURM_ARRAY_TASK_ID}.krona.html
