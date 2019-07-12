source /apps/anaconda3/etc/profile.d/conda.sh

# Activate environment
conda activate omics  # note: the environment was renamed from "genomics" to "omics"

# Note: We won't be using the regular 01_experiment, 02_experiment etc. structure here,
# because the script contains only specific well-tested commands, and it's easier
# to keep them in one place. The script will save important results in separate directories.

# This script should be executed in an empty directory.




###### 1. Generate metagenomic reads

mkdir -p data && cd data

# https://github.com/HadrienG/InSilicoSeq
# https://insilicoseq.readthedocs.io/en/latest/iss/generate.html
# Simulate 2,000,000 Mio MiSeq reads from 2 bacterial and 1 archaeal genome, selected randomly from RefSeq
# (Output name "ncbi_lc" means: data from ncbi, "low complexity" metagenome)
iss generate --ncbi bacteria archaea --n_genomes_ncbi 2 1 --n_reads 2000000 --model MiSeq --output ncbi_lc --cpus 8

ll -ht  # look at results
# Ok, we got paired-end Illumina MiSeq reads, a multifasta file with the downloaded genomes
# and an abundance file. Here is some information on paired-end sequencing:
# https://emea.illumina.com/science/technology/next-generation-sequencing/paired-end-vs-single-read-sequencing.html

# Lets look at the downloaded sequences and the abundance file, which gives the relative abundances
# of the organisms in the simulated metagenome:
grep ">" ncbi_lc_ncbi_genomes.fasta
cat ncbi_lc_abundance.txt
# We see that the abundance file has a very simple text format with two columns: sequence
# identifier (=sequence identifier in the fasta file) and abundance.
# This is useful, because it allows us to easily include genomes and set abundances ourselves.


### Bonus section: Generate metagenomic reads from selected genomes

# Let's include some interesting genomes that we select ourselves.
# E.g. we can search https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/
# (sorted by release date and assembly completeness) and pick some recently sequenced
# genomes, or genomes from lineages we are particularly interested in.
#
# Here is a small selection of possibly interesting genomes (you are free to pick your own):
# (To download the data, we need the assembly accession identifiers)
# Duncaniella sp. C9: https://www.ncbi.nlm.nih.gov/genome/79406?genome_assembly_id=495902 
#   -> click on RefSeq accession -> look for the "Assembly" field: GCF_004803935.1
# Elizabethkingia sp. 2-6: https://www.ncbi.nlm.nih.gov/genome/36583?genome_assembly_id=516671 -> GCF_005234115.1
# Metallosphaera prunae: https://www.ncbi.nlm.nih.gov/genome/79563?genome_assembly_id=515905 -> GCA_005222525.1
# Haloprofundus sp. MHR1: https://www.ncbi.nlm.nih.gov/genome/72000?genome_assembly_id=515895 -> GCF_005155585.1
# E. coli O104:H4 (EHEC strain): https://www.ncbi.nlm.nih.gov/genome/?term=O104%3AH4 -> GCF_000299455.1

# First, we need to download the genomes.
# This can be done by hand: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/
# But it can also be done on the command line, e.g. using Biopython or a tool.
#
# Download several genome assemblies from RefSeq in FASTA format
ncbi-genome-download all --assembly-accessions GCF_004803935.1,GCF_005234115.1,GCF_005155585.1,GCF_000299455.1 -F fasta
# This tool needs to know if a sequence is not in the RefSeq, but in the Genbank database
ncbi-genome-download all --assembly-accessions GCA_005222525.1 -F fasta -s genbank

# Put the downloaded data into a separate folder
mkdir -p ncbi
mv refseq genbank ncbi

# Ok, lets see what we got
find ./ncbi -name '*fna.gz' | xargs zcat | grep ">"
# NZ_CP039833.1 Haloprofundus sp. MHR1 
# NC_018658.1 Escherichia coli O104:H4 str. 2011C-3493
# NZ_CP039547.1 Duncaniella sp. C9
# NZ_CP039929.1 Elizabethkingia sp. 2-6
# CP031156.1 Metallosphaera prunae strain Ron 12/II
# For our E.coli strain, we also get three plasmids, which we can
# a) disregard (not use them as input for InSilicoSeq)
# b) concatenate with the main chromosome
# c) explicitely set their abundances in the abundance file
# We will concatenate them with the main chromosome, just to demonstrate.

# Concatenate all E.coli sequences into a single sequence
# using the "union" tool (from the EMBOSS package)
zcat ncbi/refseq/bacteria/GCF_000299455.1/*fna.gz > Ecoli_tmp.fa
seqstats Ecoli_tmp.fa  # should contain 4 sequences (chromosome + 3 plasmids)
union -filter Ecoli_tmp.fa | gzip > Ecoli_combined.fna.gz
seqstats Ecoli_combined.fna.gz  # should contain 1 sequence
rm ncbi/refseq/bacteria/GCF_000299455.1/*fna.gz
mv Ecoli_combined.fna.gz ncbi/refseq/bacteria/GCF_000299455.1/

# Concatenate downloaded sequences into a multifasta file (-> input for InSilicoSeq)
find ncbi -name '*fna.gz' | xargs zcat > selected_ncbi_genomes.fasta

# Look at the resulting file
grep ">" selected_ncbi_genomes.fasta
seqstats selected_ncbi_genomes.fasta

# Write abundance file with some made-up abundances (can also be done in a text editor)
cat <<EOF > selected_ncbi_genomes_abundance.txt
NZ_CP039833.1	0.2
NC_018658.1	0.1
NZ_CP039547.1	0.25
NZ_CP039929.1	0.4
CP031156.1	0.05
EOF

# Simulate reads from our selected genomes and our abundance file
iss generate --genomes selected_ncbi_genomes.fasta --abundance_file selected_ncbi_genomes_abundance.txt --n_reads 2000000 --model MiSeq --output sel_ncbi_lc --cpus 8

# We successfully simulated a second metagenome from our selected genomes.
### End Bonus section




##### 2. Quality control (pre-processing, trimming/filtering)

# Set read names as variables (so you can quickly change them later if you want)
cd ..
IN=data/ncbi_lc_R1.fastq
IN2=data/ncbi_lc_R2.fastq
# If you named your output files differently, set the variables accordingly.

# Basic statistics: seqstats
seqstats $IN
seqstats $IN2

# More elaborate statistics: seqtk fqchk
seqtk fqchk $IN
seqtk fqchk $IN2
# This tool looks at the frequncy of all nucleotides at all positions, so if you
# had primers/adapters/contaminations, you should see irregularities.
# Note the avgQ and the errQ fields, which give the average read quality calculated
# in two different ways (explained in https://blog.liang2.tw/posts/2015/09/seqtk/).
# The "quality scores" are determined by the sequencer and try to reflect the
# probability of sequencing error: https://en.wikipedia.org/wiki/Phred_quality_score

# Graphical report: FastQC
mkdir -p 01_fastqc_results
fastqc --outdir 01_fastqc_results $IN $IN2
# The resulting html file can be opened and inspected in any browser.
#
# The "Per sequence quality scores" show how the quality scores are distributed,
# e.g. if you have groups of reads with with low quality scores that you can filter
# away by setting a corresponding treshold during quality filtering.
#
# The "Per base sequence content" pane gives similar information as seqtk fqchk,
# the frequencies of nucleotides at all read positions.
#
# The "Per sequence GC content" might show several peaks: this is what you would
# expect in a metagenome with multiple organisms with differing GC content.
# (https://www.omicsonline.org/open-access/the-gc-content-of-bacterial-genomes-2329-9002-2-e108.php?aid=26236)
# If this wasn't a metagenome, this would point to contamination.
#
# The "Sequence duplication level" should be low; you wouldn't expect many duplicated
# sequences in a shotgun metagenome, unless the coverage is very high (and the community
# is low complexity). High sequence duplication levels might point to some technical
# problems, e.g. adapter contamination.
#
# The "Overrepresented sequences" section checks for overrepresented sequences like adapters.
# 
# FastQC has very nice documentation, with helpful information and examples:
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/


# Pre-processing (trimming and filtering) can be done to:
# a) get rid of contaminations like adapters (if there are contaminations),
# b) merge paired-end sequences if the insert size is short and sequences overlap (https://emea.illumina.com/science/technology/next-generation-sequencing/paired-end-vs-single-read-sequencing.html, https://images.app.goo.gl/Hd3HRkM7Yn7FT2CR9) and
# c) trim/filter reads to remove low-quality bases (from 3' end) and short or
# low-quality sequences.
#
# Not all subsequent tools require sequence preprocessing, so it depends on what you want to do afterwards.

# Two examples describing metagenome preprocessing from the "Methods" section of papers:
# "Following sequencing, raw reads from each metagenome and metatranscriptome were filtered to remove reads containing adapters, reads containing > 10% ambiguous bases, and reads containing low quality bases (Q-score<= 5) over 50% of the total bases" (paper: "Wastewater treatment plant resistomes are shaped by bacterial composition, genetic exchange, and upregulated expression in the effluent microbiomes")
# "The SYNTH, HMP, MARINE, and SOIL datasets were pre-processed to remove adaptors and trim low quality segments of the reads. We used cutadapt software v 1.9.1 (Martin 2011), trimming bases with PHRED quality < 20 from 3’ end (parameter -q 20). Adaptor sequences were identified for each dataset individually, using FastQC v0.11.3 and manual reads inspection" (paper: "metaSPAdes: a new versatile metagenomic assembler")

# Example for removing adapters combined with quality-based trimming (prior to merging)
#cutadapt -q 20 --cores 4 -a GAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o read1_qf.fq -p read2_qf.fq read1.fq read2.fq

# Example for merging reads
#bbmerge.sh in=reads.fq in2=read1.fq out=merged.fq outu=unmerged.fq ihist=ihist.txt

# Filtering, trimming, removing contaminants: bbduk, trimmomatic, seqtk
# 
# Often, trimming/filtering is done based on 3 criteria:
# a) trim low-quality bases from 3' end (e.g. q<10)
# b) discard reads with average quality below some value (e.g. 10)
# c) discard reads which became too short after trimming (e.g. <50% of original length)
# You should consider more or less stringent trimming based on data quality,
# data size and which tools are used subsequently (e.g. high-quality sequences are more relevant
# for assembly than for alignment: a partly erroneous read might align, but not assemble.)

# For our reads, FastQC shows that the quality substantiall decreases at the 3' end.
# Do some basic quality trimming
# https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
bbduk.sh in1=$IN in2=$IN2 out1=$IN.bbduk.fq out2=$IN2.bbduk.fq qtrim=r trimq=10 maq=10 minlen=100
# As usual, we should inspect the output
seqstats $IN.bbduk.fq
seqstats $IN2.bbduk.fq
# It's also a good idea to run FastQC again on the filtered dataset to double-check. You can do it 
# yourself, if you want.

