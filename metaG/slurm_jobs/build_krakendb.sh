#!/bin/bash
#
#SBATCH --job-name=krakendb-build
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
module load kraken/2.0.7

#build database 
# Initialize an empty kraken database in the fast, local temporary directory if this job
kraken-build --download-taxonomy --db $TMPDIR/krakendb
# Add genomes
kraken-build --download-library bacteria --db $TMPDIR/krakendb
# Build the database index
kraken-build --build --db $TMPDIR/krakendb
# Remove all temporary files, except the final kraken database
kraken-build --clean --db $TMPDIR/krakendb
# Move the final kraken database into a network folder, from where it can be used for further kraken jobs
mv $TMPDIR/krakendb /proj/DECOMB/analysis/process_metaG/kraken2
