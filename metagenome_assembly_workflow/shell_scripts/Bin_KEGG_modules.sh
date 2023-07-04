#!/bin/bash
#
#SBATCH --job-name=bin_KOfam_assignment
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --output=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.out
#SBATCH --error=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio/Log/%x-%j.err

#Set up the path to the working directory
WORKDIR=/scratch/oceanography/efadeev/DECOMB/analysis/metaG_anvio
cd $WORKDIR

readarray -t BINS < $WORKDIR/06_BINS/Refined_bins.txt

for bin in ${BINS[@]}; do

anvi-run-kegg-kofams -c $WORKDIR/06_BINS/REFINED/$bin/CONTIGS.db \
--kegg-data-dir /proj/DECOMB/source/KOfam -T 20 \
--just-do-it
done

#produce metabolic metrices for visualization
anvi-estimate-metabolism -e $WORKDIR/07_METABOLISM/selected-bins-collection.txt \
                         -O $WORKDIR/07_METABOLISM/Bins \
                         --matrix-format --kegg-data-dir /proj/DECOMB/source/KOfam

#produce interactive visualization
anvi-matrix-to-newick $WORKDIR/07_METABOLISM/Bins-completeness-MATRIX.txt

# dry run to get the profile db:
anvi-interactive -d $WORKDIR/07_METABOLISM/Bins-completeness-MATRIX.txt \
                 -p $WORKDIR/07_METABOLISM/Bins_metabolism_PROFILE.db \
                 --manual-mode \
                 --dry-run

# import the state file (changes from the Anvio tutorial)
anvi-import-state -s $WORKDIR/07_METABOLISM/bins-heatmap-layout.json \
                  -p $WORKDIR/07_METABOLISM/Bins_metabolism_PROFILE.db \
                  -n default

#make a list of the relevant metabolic pathways

# learn where the MODULES.db is:
export ANVIO_MODULES_DB=`python -c "import anvio; import os; print(os.path.join(os.path.dirname(anvio.__file__), '/proj/DECOMB/source/KOfam/MODULES.db'))"`

# start an empty file:
echo -e "module\tclass\tcategory\tsubcategory\tname" > $WORKDIR/07_METABOLISM/modules_info.txt

# get module classes:
sqlite3 $ANVIO_MODULES_DB "select module, data_value from kegg_modules where data_name='CLASS'" | \
    sed 's/; /|/g' | \
    tr '|' '\t' >> $WORKDIR/07_METABOLISM/module_class.txt

# get module names:
sqlite3 $ANVIO_MODULES_DB "select module,data_value from kegg_modules where data_name='NAME'" | \
    tr '|' '\t' > $WORKDIR/07_METABOLISM/module_names.txt

# join everything
paste $WORKDIR/07_METABOLISM/module_class.txt <(cut -f 2 $WORKDIR/07_METABOLISM/module_names.txt ) >> $WORKDIR/07_METABOLISM/modules_info.txt

# empty the trash bin:
rm $WORKDIR/07_METABOLISM/module_names.txt $WORKDIR/07_METABOLISM/module_class.txt

#import modules to the profile
anvi-import-misc-data $WORKDIR/07_METABOLISM/modules_info.txt -p $WORKDIR/07_METABOLISM/Bins_metabolism_PROFILE.db -t items

#produce tables with KO hits and modules
anvi-estimate-metabolism -e $WORKDIR/07_METABOLISM/selected-bins-collection.txt\
                         -O $WORKDIR/07_METABOLISM/Bins \
                         --kegg-data-dir /proj/DECOMB/source/KOfam \
                         --kegg-output-modes kofam_hits,modules
