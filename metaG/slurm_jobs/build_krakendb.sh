#!/bin/bash
#
#SBATCH --job-name=krakendb-build
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.out
#SBATCH --error=/proj/DECOMB/analysis/process_metaG/Log/%x-%j.err

#load module
# Build the database index
module load kraken/2.0.8-beta

mkdir $TMPDIR/krakendb

tar xvfz ./kraken2/kraken_ref/k2_standard_20201202.tar.gz --directory $TMPDIR/krakendb

kraken2-build --build --threads 12 --db $TMPDIR/krakendb
# Remove all temporary files, except the final kraken database
kraken2-build --clean --threads 12 --db $TMPDIR/krakendb
# Move the final kraken database into a network folder, from where it can be used for further kraken jobs
mv $TMPDIR/krakendb /proj/DECOMB/analysis/process_metaG/kraken2
