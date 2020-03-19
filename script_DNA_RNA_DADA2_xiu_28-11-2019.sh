# This script is for the 16S amplicon sequence analysis (515F-806R) of DNA/cDNA analysis
# Sequeincing on Illumina Miseq platform (150*2), Argonne (DNA and cDNA samples run seperately)
# RNA-run 60 samples; DNA-run 104 samples, ratio=1.7; where after demux, sequence counts in total: DNA 7,798,981; RNA 9,852,975, ratio=1.44
# Updated on 30-11-2019 by Jia Xiu

# if submit job to Peregrine HPC as .sh file, i.e. import_multiplex_seq.sh, add following header to .sh file
#!/bin/bash
#SBATCH --job-name=
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=12GB
#SBATCH -o import-%j.out
#SBATCH -e import-%j.error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=you-email@rug.nl

# change directory
cd $PWD


# using keemei to validate sample-metadata.tsv files


# load QIIME2
module load QIIME2/2019.10


# demutiplex raw sequences get from argonne

# Import files interactively
# Import the multiplexed sequences

cd raw_seq/dna-demux-sequences/

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path MANIFEST-dna \
  --output-path paired-end-demux-dna.qza \
  --input-format PairedEndFastqManifestPhred33

cd /raw_seq/cdna-demux-sequences

sed 's/\t/,/g' MANIFEST-cdna > MANIFEST-cdna.txt

mv MANIFEST-cdna.txt MANIFEST-cdna 

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path MANIFEST-cdna \
  --output-path paired-end-demux-rna.qza \
  --input-format PairedEndFastqManifestPhred33

mv paired-end-demux-rna.qza ../../

# view a summary
qiime demux summarize \
  --i-data paired-end-demux-dna.qza \
  --o-visualization paired-end-demux-dna.qzv

qiime demux summarize \
  --i-data paired-end-demux-rna.qza \
  --o-visualization paired-end-demux-rna.qzv


# Warning: the DADA2 denoising process is only applicable to a single sequencing run at a time, so we need to run this on a per sequencing run basis and then merge the results.



# DADA2 sequence quality control

# DNA
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux-dna.qza \
  --o-table table-dna.qza \
  --o-representative-sequences rep-seqs-dna.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --p-chimera-method consensus \
  --o-denoising-stats stats-dada2-dna.qza

# FeatureTable and FeatureData summarize
qiime feature-table summarize \
  --i-table table-dna.qza \
  --o-visualization table-dna.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dna.qza \
  --o-visualization rep-seqs-dna.qzv

qiime metadata tabulate \
  --m-input-file stats-dada2-dna.qza \
  --o-visualization stats-dada2-dna.qzv


# RNA dataset
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux-rna.qza \
  --o-table table-rna.qza \
  --o-representative-sequences rep-seqs-rna.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --p-chimera-method consensus \
  --o-denoising-stats stats-dada2-dna.qza

# FeatureTable and FeatureData summarize
qiime feature-table summarize \
  --i-table table-rna.qza \
  --o-visualization table-rna.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-rna.qza \
  --o-visualization rep-seqs-rna.qzv

qiime metadata tabulate \
  --m-input-file stats-dada2-rna.qza \
  --o-visualization stats-dada2-rna.qzv


# Merging denoised data, i.e. two feature table and rep_seq files
qiime feature-table merge \
  --i-tables table-dna.qza \
  --i-tables table-rna.qza \
  --o-merged-table table.qza

qiime feature-table merge-seqs \
  --i-data rep-seqs-dna.qza \
  --i-data rep-seqs-rna.qza \
  --o-merged-data rep-seqs.qza

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv


# exporting a nontaxa OTU table
qiime tools export \
  --input-path table.qza \
  --output-path exported/

# convert biom to txt
#biom convert -i nontaxonomic-otu-table/feature-table.biom -o OTU-table-nontax.tsv --to-tsv




### Taxonomic analysis ###

# assign taxonomy by silva database
qiime feature-classifier classify-sklearn \
  --i-classifier /data/p278113/QIIME2/silva-132-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# summerize taxonomy info
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Export taxonomy
qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported_taxonomy



### Remove mitochondria, cloroplasti, archaea and keep sequence assigned at phyla level (D_0_ for SILVA database) ###

# filter representative sequences

# filter feature table
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include D_0__ \
  --p-exclude archaea,eukaryota,mitochondria,chloroplast \
  --o-filtered-table table-filtered.qza

qiime feature-table filter-features \
  --i-table table-filtered.qza \
  --p-min-frequency 2 \
  --o-filtered-table table-filtered.qza

# filter representative sequences
qiime taxa filter-seqs \
  --i-sequences rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-include D_0__ \
  --p-exclude archaea,eukaryota,mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-filtered.qza

# summerize filtered FeatureTable and FeatureData
qiime feature-table summarize \
  --i-table table-filtered.qza \
  --o-visualization table-filtered.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-filtered.qza \
  --o-visualization rep-seqs-filtered.qzv

qiime tools export \
  --input-path table-filtered.qza \
  --output-path exported-table

biom convert -i exported-table/feature-table.biom -o feature-table-filtered.tsv --to-tsv


### generate a phylogenetic tree using the filtered represented sequences ###
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-filtered.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported-tree


### Alpha and beta diversity ###
# rarefy
# check --p-sampling-depth in the table-filtered.qzv
qiime feature-table rarefy \
    --i-table table-filtered.qza \
    --p-sampling-depth 15000 \
    --o-rarefied-table table-rarified.qza

# rarified FeatureTable and FeatureData summarize
qiime feature-table summarize \
  --i-table table-rarified.qza \
  --o-visualization table-rarified.qzv \
  --m-sample-metadata-file sample-metadata.tsv

# exporting a rarified feature table
qiime tools export \
  --input-path table-rarified.qza \
  --output-path exported-table/  # extracted will will do this verbosely

# convert biom to txt
biom convert -i exported-table/feature-table.biom -o feature-table-rarified-nontax.tsv --to-tsv


# diversity analysis:  # check --p-sampling-depth in the table-filtered.qzv
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-rarified.qza \
  --p-sampling-depth 15000 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results

# exporting matrix
qiime tools export \
  --input-path core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --output-path exported-matrix

mv exported-matrix/distance-matrix.tsv ./unweighted_unifrac_distance_matrix.tsv

qiime tools export \
  --input-path core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --output-path exported-matrix

mv exported-matrix/distance-matrix.tsv ./weighted_unifrac_distance_matrix.tsv

#rarefaction curve
qiime diversity alpha-rarefaction \
  --i-table table-rarified.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 15000 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv \
  --verbose

### statistics ###
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac-elevation-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/weighted-unifrac-elevation-significance.qzv \
  --p-pairwise

qiime taxa barplot \
  --i-table table-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

# Taxa summary
qiime taxa collapse \
  --i-table table-rarified.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-level-6.qza

qiime composition add-pseudocount \
  --i-table table-level-6.qza \
  --o-composition-table table-6.qza

qiime composition ancom \
  --i-table table-6.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column subject \
  --o-visualization ancom-subject.qzv
