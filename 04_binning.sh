source /apps/anaconda3/etc/profile.d/conda.sh

conda activate omics
###### Binning
#
# In this step, we try to group the assembled DNA (contigs and scaffolds) based on their genomic origin.
# If this works well, we can treat the bins as (possibly incomplete) genomes, which is pretty common nowadays.
# E.g., if you look at the NCBI genome list https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/,
# you can see a column giving the "assembly level" for every organism. If you look at some assembly
# more closely, e.g. https://www.ncbi.nlm.nih.gov/assembly/GCA_003181115.1, you will see that this genome
# assembly consists of contigs. Even though the genome representation is marked as complete, it didn't
# fulfill the RefSeq database quality creteria, because it was derived from a metagenome, so it's
# available only in the GenBank database. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965038/ gives more
# details about the criteria for a genome to make it into RefSeq, and explains some RefSeq terms.)
#
# As with the other steps, there is a plethora of tools for this task, which use several sources
# of information. Mostly, this information is:
# - Nucleotide composition of the contigs/scaffolds (e.g. frequency of k-mers, usually tetramers), e.g.
# MetaBAT, CONCOCT, MaxBin2, MetaWatt
# - Contig/scaffold coverage; this is even more powerful if multiple samples are available, where the same
# organisms occur in varying frequencies, so that every sample has its own coverage values when projected on
# one single assembly. E.g. contig_1 might have a 10x coverage in sample_1, and a 50x coverage in sample_2,
# because the corresponding organism is 5 times more abundant in the second sample. Now, if contig_2
# coverage also increased ~5x for sample_2 (and has a similar nucleotide composition), chances are high that
# those two contigs belong to the same organism; otherwise they probably belong to different organisms.
# Tools like MetaBAT or CONCOCT make use of this information, and using multiple samples (e.g. from nearby
# locations or different sampling/DNA extraction methods) can substantially improve the binning quality. 
#
# Binning based on sequence similarity to known genes usually doesn't work well, because we can encounter a lot
# of novel genes/sequences in metagenomes, and because horizontal gene transfers (HGT) can confound the results.
# However, conserved genes have low levels of HGTs and are universal, so they can help determine bin taxonomy
# and estimate genome completes; this concept is used in tools like CheckM or BUSCO (https://www.ncbi.nlm.nih.gov/pubmed/29220515).


#IN=contigs.filtered.fa
IN=megahit_out/final.contigs.fa.filtered.fasta


### Let's take a quick look at the assembly, just to get a better feeling for the data
stats.sh in=$IN


### MaxBin2
# https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html
mkdir -p maxbin_out
OUT=maxbin_out/bin
run_MaxBin.pl -contig $IN -out $OUT -reads read1.fq -reads2 read2.fq -thread 8 -markerset 40
#run_MaxBin.pl -contig $IN -out $OUT -abund abundance.txt -thread 8 -markerset 40

# Does the abundance of the bins match the 16S profile of the community?


### MetaBAT
# https://bitbucket.org/berkeleylab/metabat/wiki/Example_Large_Data
# https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices
OUT=metabat_out/bin
# We need bam files (mapping files), which can be generated according to the wiki
# or taken from a previous step, if we already did this.
#
# Align reads to assembly
#bbwrap.sh ref=$IN in=read1.fq in2=read2.fq out=aln.bam
#samtools sort aln.bam > aln.sorted.bam && rm aln.bam
#
# Use existing mapping files and convert them to bam-files
SAM=aln.sam.gz
zcat $SAM | samtools view -S -b > tmp.bam
samtools sort tmp.bam -o aln.sorted.bam && rm tmp.bam
samtools index aln.sorted.bam
mkdir -p bamdir
mv *bam *bam.bai bamdir
# (You probably should run `samtools --help` and `samtools view --help` etc. to understand
# what the commands are doing. Also, as usual after every step, you should look at the
# directory and file contents using commands like `ll -ht`, `head` and `less`, to see
# what's going on. E.g., in the course "Spezielle Werkzeuge für das QM in der Datenanalyse",
# one student complained that one of the Python scripts didn't work, while the problem was
# a broken symlink (pointing to a non-existing file). It would have been easy to spot the
# problem by looking at the directory with the symlink using `ll -h`.)
jgi_summarize_bam_contig_depths --outputDepth depth.txt --pairedContigs paired.txt --minContigLength 1000 --minContigDepth 2 bamdir/*.bam
# MetaBAT makes some adjustments to regular coverage calculation, which are described in
# the README, https://bitbucket.org/berkeleylab/metabat/src/master/README.md
metabat -i $IN -a depth.txt -o $OUT -v --saveTNF saved.tnf --saveDistance saved.dist

# Do the binning results differ between MaxBin and MetaBAT?

# There are many more binners, or a tool that can combine results from different binners:
# https://github.com/cmks/DAS_Tool, but we are not going to use it here. In
# case you ever need to do binning for your own work, you might want to try it out.
# It's important to remember that new tools/approaches do not necessarily produce better
# results than old ones, and you should probably evaluate such tools on simulated datasets
# first to see how well they work for you.


### Evaluate binning results
# Now, it's useful to assess how well the binning worked, i.e. how complete the genome
# bins are and if they contain only genomic sequences from one or from many organisms 
# (contamination). This also depends on how related the organisms in the metagenome
# are: It is probably more difficult to clearly distinguish multiple Escherichia species
# than Escherichia from an Alphaproteobacterium.
# Tools designed for this task like CheckM and BUSCO usually look at the presence and
# abundance of markergenes.

### Checkm
# https://github.com/Ecogenomics/CheckM/wiki/Workflows
conda activate checkm
IN=metabat_out/
OUT=checkm_results
checkm lineage_wf -t 20 -x fa $IN $OUT > checkm_summary.txt


### You can also get a visual overview over binning results using a tool like VizBin
# (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4305225/)


### What next
# The bins can be treated as (more or less incomplete) genomes, and courses like "Vergleichende Genomik"
# and "Genomanalyse" might come in handy now.
# You should also not forget that maybe some organisms couldn't be assembled well, and double-check
# for discrepancies between your read-based taxonomic profile and the bin taxonomy from the bins.
# E.g., if you see an organism in the read-based profiling that doesn't appear in the bins,
# maybe this organism couldn't be assembled.
# In general, the data analysis strategy depends mostly on the project goals, and you should
# try to find studies which addressed a similar type of questions and look at how they were performed.

# A nice overview of the approaches which we learned in this courses can be found in the publication
# "A review of methods and databases for metagenomic classification and assembly", which I will upload to moodle.

