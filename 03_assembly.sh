source /apps/anaconda3/etc/profile.d/conda.sh

###### Assembly
#
# In this step, we want to combine short NGS reads into longer DNA sequences (contigs and scaffolds),
# which makes gene calling, functional and taxonomic analysis and possibly even the reconstruction of 
# complete genomes easier.
#
# Suggested assemblers for metagenomic shotgun assembly: MetaSPAdes, Megahit, IDBA-UD
# (One difference which can be relevant to real-life application are memory requirements and processing time,
# e.g. Megahit vs. Meta-SPAdes.)
# Some basic information: https://www.ncbi.nlm.nih.gov/assembly/basics/


### Megahit assembler: https://github.com/voutcn/megahit/wiki/An-example-of-real-assembly
# For assembly, you should definitely use quality-controlled reads. Please modify the input files accordingly.
conda activate omics

# Assemble the metagenome using the Megahit assembler
OUT=megahit_out
megahit -1 read1.fq -2 read2.fq -o megahit_out

# Explanation of Megahit fasta headers (https://github.com/voutcn/megahit/issues/54):
# "flag is a tag to represent the connectivity of a contig in the assembly graph. flag=1 means the contig is standalone, flag=2 a looped path and flag=0 for other contigs
# multi is roughly the average kmer coverage. But the figure is not precise and if you want to quantify the coverage of a contig, I would suggest you align the reads back to the contigs by short read aligners.
# Basically, you could just ignore the header and put the contigs to the subsequent analysis."


### Spades assembler: https://github.com/ablab/spades
#OUT=metaspades_out
#metaspades.py -1 read1.fq -2 read2.fq -t 16 -o $OUT --phred-offset 33
# time: 1,000,000 reads on 8 cpus -> real: 13m (user: 70m)


### Evaluate assembly
#
# It's useful to know how well the assembly worked and to get some basic statistics
# on the assembled contigs. We are in luck, as there are automated tools for that.

# Quast: http://quast.sourceforge.net/quast
# http://quast.bioinf.spbau.ru/manual.html
conda activate omics
IN=megahit_out/final.contigs.fa
#IN=metaspades_out/contigs.fasta
seqstats $IN
quast -1 read1.fq -2 read2.fq $IN
#metaquast.py $IN

# Metaquast tries to determine closely related organisms and download them from NCBI.
# Not sure how reliably this works, and it can also take pretty long, so we can also use
# the regular quast instead (which is designed for single-genome assembly, but can still
# provide useful information). You can try both.
# Quast also has a lot of other useful options like gene calling and more.
#
# The results are saved as "report.pdf" and "report.html", they can be downloaded and 
# examined locally. (When a file is downloaded via the Jupyter server, it is opened in
# a new browser tab; this doesn't always work, so it can be saved using a right mouse
# click -> "Save as..." and the re-opened, this seems to work better. Anyway, it might
# be a better idea to use other file downloading methods instead, like Filezilla or scp
# on the command line (http://www.hypexr.org/linux_scp_help.php).


### Calculate contig coverage
# (and determine fraction of assembled reads; extract unassembled reads)
#
# The contig coverage helps to assess the sequencing depth (were also low-abundance organisms/
# variants sequenced?), the reliabiliy of the assembly (higher coverage usually means
# better assembly), and can help with the binning, because we may be able to differentiate
# between taxonomic groups with similar genome composition by their abundance/coverage.

# https://github.com/voutcn/megahit/wiki/An-example-of-real-assembly
IN=megahit_out/final.contigs.fa
#IN=metaspades_out/scaffolds.fasta

# Align reads back to contigs/scaffolds with the short read mapper bbmap/bbwrap.sh
bbwrap.sh ref=$IN in=read1.fq in2=read2.fq out=aln.sam.gz kfilter=22 subfilter=15 maxindel=80
# Output per contig coverage to cov.txt with bbmap/pileup.sh
pileup.sh in=aln.sam.gz out=cov.txt

# We can extract unassembled reads to examine them further (maybe there was an organism which
# couldn't be assembled?), but we're not going to do this here.
# Extract unmapped reads (single-end to unmapped.se.fq and paired-end to unmapped.pe.fq)
#samtools view -u -f4 aln.sam.gz | samtools bam2fq -s unmapped.se.f

# Get important contigs stats (length/coverage/GC)
cat cov.txt | awk -F'\t' 'BEGIN {OFS = FS} {print $1, $2, $3, $9}' | sed 's/ [^\t]*//' > assembly_stats.txt

# Shorter and low-coverage contigs are often of lower quality, so it's better to get rid of them.
# Filter by length (e.g. 500bp or 1kb) and coverage (e.g. 3x)
cat assembly_stats.txt | awk -F'\t' 'BEGIN {OFS = FS} {if ($2 > 3 && $3 > 500) print $0}' > assembly_stats.filtered.txt
cat assembly_stats.filtered.txt | cut -f1 | grep -v '^#ID' > assembly_stats.filtered.contignames
# seqtk subseq can filter a fasta file based on a list of sequence names
seqtk subseq $IN assembly_stats.filtered.contignames > $IN.filtered.fasta
# Compare the original and the filtered fasta files
seqstats $IN
seqstats $IN.filtered.fasta


