source /apps/anaconda3/etc/profile.d/conda.sh

###### Taxonomic profiling
#
# The goal of this step is to get a taxonomic profile of the sequenced community:
# Who is there, and how many of them?
#
# "Overall, pairing tools with different classification strategies (k-mer, alignment, marker) can combine their respective advantages." ("Comprehensive benchmarking and ensemble approaches for metagenomic classifiers")
#
# Different approaches are based on:
# - complete genomic sequence -> e.g. KrakenUniq (unique k-mers)
# - universal tax. markers (rRNA) -> e.g. Metaxa2
# - single-copy universal tax. markers (proteins) -> e.g. mOTUs2
# - clade-specific tax. markers -> e.g. MetaPhlAn2



### Preparations
# Let's create symbolic links, so it'll be easy to repeat the analysis using other
# input files without having to change the file names throughout the script.
ln -s data/ncbi_lc_R1.fastq read1.fq
ln -s data/ncbi_lc_R2.fastq read2.fq
#ln -s other_data/another_dataset_R1.fq read1.fq
#ln -s other_data/another_dataset_R2.fq read1.fq

# For testing tools and trying out things, it's a good idea to subsample our data
# and use it instead of the original data files to prevent unnecessary waiting times.
# E.g., we can subsample to 100,000 reads (should be good enough to run quickly and
# still produce meaningful results)
conda activate omics

seqtk sample -s 100 read1.fq 100000 > sub1.fq
seqtk sample -s 100 read2.fq 100000 > sub2.fq
# If a tool works and does the right thing, we should run it on the complete dataset
# (unless this would take too long, then a subsample will have to do).


### Profiling based on the whole genome sequence

# KrakenUniq: https://github.com/fbreitwieser/krakenuniq/blob/master/MANUAL.md
#
# You'll see that if a closely related organism is in the database, this tool works pretty well.
# So, to properly test it, we should use completely novel genomes and/or remove close relatives
# from the database. We are not going to do this here, but the results won't be
# representative, because we used Refseq genomes both for simulating our metagenome and
# for building the database that KrakenUniq uses for classification.

conda activate omics  # set the right environment

# Prepare KrakenUniq database: This was already done, please don't do it again
#DBDIR=/mirror/kraken_db
#mkdir -p $DBDIR
# Download the taxonomy
#krakenuniq-download --db $DBDIR taxonomy
# Download complete bacterial and archaeal genomes from RefSeq (and mask low-complexity regions)
#krakenuniq-download --db $DBDIR --threads 10 --dust refseq/bacteria refseq/archaea
# Build database
#krakenuniq-build --db $DBDIR --kmer-len 31 --threads 10 --taxids-for-genomes --taxids-for-sequences

# Classify reads with KrakenUniq
DBDIR=/mirror/kraken_db
krakenuniq --db $DBDIR --threads 10 --report-file kraken_taxonomy_profile.txt --paired read1.fq read2.fq  > kraken_read_classification.tsv
# The problem we encountered last time when we tried to run it was related to a broken conda
# package (described in https://github.com/tseemann/nullarbor/issues/219).
# This should be fixed now, but it's good to know that (sometimes) things like that can happen.

# Visualize/evaluate results
# Pavian (https://github.com/fbreitwieser/pavian) is really nice for visualizing the results,
# you can run it in RStudio server:
# "pavian::runApp(port=5000)"
# (A new window should open, where you will need to load the obtained taxonomy profile to get
# a visualization of the taxonomic composition of the sample.)



### Profiling based on rRNA genes

# Metaxa2: https://microbiology.se/software/metaxa2/ 
conda activate metaxa

# Classify reads with Metaxa
OUT=metaxa_out
metaxa2 -1 read1.fq -2 read2.fq -g ssu --mode metagenome --plus T --cpu 8 --megablast T -o $OUT
# Btw, would blasting random sequences against the Metaxa database take the same amount of time?
# Or, in other words: why does this blast take so long ("this may take a long while")?

# Overview of the results
cut -f2 $OUT.taxonomy.txt | sort | uniq -c

# Filtering the Metaxa output for long alignments can improve the results by getting rid of unspecific 
# (short) alignments - but you should try it yourself and see how it impacts the results.
# Let's say, we want to include only reads that aligned over 80% of their length; for a ~300 bp long read, it's 240 bp
CUTOFF=240
cat $OUT.taxonomy.txt | awk -F'\t' -v cutoff="$CUTOFF" '$4 > cutoff {print $0}' | cut -f2 | sort | uniq -c
cat $OUT.taxonomy.txt | awk -F'\t' -v cutoff="$CUTOFF" '$4 > cutoff {print $0}' > $OUT.taxonomy.filtered.txt

# There are different ways to visualize, let's convert the taxonomic lineages to taxids and use Krona.
# There is a script for that, but to make sure that it works correctly, please double-check
# that the resulting taxids match to the taxonomic lineages given by Metaxa.
conda activate omics

python /apps/scripts/lineage2taxid.cpython-36.pyc --mail gruber.matthias@gmx.at --db /apps/etetoolkit/taxa.sqlite -i $OUT.taxonomy.filtered.txt -c2 > $OUT.taxonomy.filtered.taxids.txt

# Visualize with Krona
ktImportTaxonomy $OUT.taxonomy.filtered.taxids.txt -o $OUT.taxonomy.filtered.krona.html -q 1 -t 2 -s 5
# Now we can download $OUT.taxonomy.filtered.krona.html to our computer and open it in browser



### Profiling based on universal single-copy marker genes (proteins, not rRNA genes)
#
# Which advantages have universal single-copy genes compared to other genes? Mostly two:
# - More precise estimation of taxon abundance because of known copy number
# - But the really important reason is that horizontal gene transfer (HGT) is rare for single-copy genes.
# So, if we get a hit to a known gene in the database, we can be pretty sure that the
# query sequence originates from a related organism, and wasn't aquired via HGT from a different one.
# Also, chances are lower that it's an unspecific hit, because the gene in question is universal
# and thus probably present in the newly discovered organism, too.

# mOTUs2: https://motu-tool.org/, https://github.com/motu-tool/mOTUs_v2
# This tool was trained mostly with metagenomes from the human body and
# ocean water (https://www.nature.com/articles/s41467-019-08844-4#Sec7),
# so we might expect that it won't work too well e.g. for soil-related microbes.
# (In the paper introduction, they critize approaches like rRNA gene-based
# profiling and clade-specific profiling, and explain the advantages of their
# tool - please include this in your report.)
#
# This tool is very new and it would be really interesting to see how well it performs
# e.g. compared with Metaxa, so please record and discuss all your observations.

# Profile metagenome with mOTUs2
# (We should used quality-controlled reads as input for this tool -
# please set the file names accordingly)
OUT=motus_taxonomy_profile.txt
motus profile -f read1.fq -r read2.fq -t 30 > $OUT

# Get only lines where abundance is > 0
cat $OUT | awk -F'\t' 'BEGIN {OFS = FS} /^#/ {print $0} {if ($2 > 0) print $0}' > motus_taxonomy_profile_filtered.txt

# On https://motu-tool.org/, you can see a Krona visualization; can you produce a corresponding
# visualization from your mOTUs output?



### Profiling based on clade-specific marker genes
# (This is pretty much the exact opposite of universal marker genes)
# What are the conditions that must be met for this approach to work well?

# MetaPhlAn: https://bitbucket.org/nsegata/metaphlan/wiki/MetaPhlAn_Pipelines_Tutorial
conda activate metaphlan  # MetaPhlAn now has its own environment

# Profile metagenome with MetaPhlAne
OUT=metaphlan_taxonomy_profile.txt
rm metagenome.bowtie2.bz2 # delete bowtie2 to prevent failure
metaphlan2.py read1.fq,read2.fq --bowtie2out metagenome.bowtie2.bz2 --nproc 8 --input_type fastq > $OUT


# Visualize with Krona
metaphlan2krona.py -p $OUT -k metaphlan_krona.out
ktImportText -o $OUT.krona.html metaphlan_krona.out



### Evaluation of the profiling results
#
# Please compare the profiling results from different methods. Did you have organisms in your metagenome
# which could be detected using one method but not another? How precise were the taxonomic assignments
# and the quantifications (abundance) in the different cases? How does this compare to MG-RAST results?
# In your report, you can use the terms:
# - "false positive" for an organisms which was identified by the tool, but absent from the dataset
# - "false negative" for an organism which was in the dataset, but was missed by the tool
# Compare https://www.nature.com/articles/srep19233 for a publication, which is in part very similar to
# what we have been doing.



