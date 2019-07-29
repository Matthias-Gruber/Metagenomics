# Metagenomik

<img src="https://github.com/fhwnmatt/Metagenomics/blob/master/figures/Abundance.png" title="Abundance">
<figcaption> Abundance </figcaption>

# Abstract

# Introduction

# Data generation

## Introduction

Im ersten Schritt werden mit dem Sequencing Simulator Tool iss (InSilicoSeq - https://github.com/HadrienG/InSilicoSeq) die benötigten paired-end Illumina MySeq Reads erzeugt. 

## Methods

Mit dem folgenden Befehl wurden 2 000 000 MiSeq reads von 2 bakteriellen und 1 Archae Genom - von RefSeq zufällig ausgewählt - simuliert.

```sh
iss generate --ncbi bacteria archaea --n_genomes_ncbi 2 1 --n_reads 2000000 --model MiSeq --output ncbi_lc --cpus 8
```

## Results and discussion

Die 2.000.000 Reads wurden in einem .fasta File und die zugehörige Abundance in einem .txt File gespeichert.
Das Abundance File enthält 2 Spalten - in der ersten den Sequence Identifier und in der zweiten die Abundance.

# Quality Control

## Introduction

Im nächsten Schritt wird die Qualität der Reads mit dem Tool FastQC (https://github.com/s-andrews/FastQC) überprüft und anschließend mit dem Tool bbduk (https://github.com/BioInfoTools/BBMap) ein Quality Trimming durchgeführt.

## Methods

Mit FastQC wurde die Qualität beurteilt.

```sh
fastqc --outdir 01_fastqc_results $IN $IN2
```

Das Quality Trimming erfolgte mit bbduk.

```sh
bbduk.sh in1=$IN in2=$IN2 out1=$IN.bbduk.fq out2=$IN2.bbduk.fq qtrim=r trimq=10 maq=10 minlen=100
```

## Results and discussion

Da die Qualität der Reads am 3' Ende abnimmt wurden sämtliche Reads mit einer Quality < 10 entfernt.
Ebenso wurden nur Reads berücksichtigt mit einer Mindestlänge von 100 Basen.

# Profiling of the community

## Introduction

In diesem Abschnitt wurden einige Tools herangezogen um ein taxonomisches Profiling des Datensatzes durchzuführen.
Folgende Tools wurden verwendet:
* Krakenuniq - https://github.com/fbreitwieser/krakenuniq
* Metaxa2    - https://microbiology.se/software/metaxa2/
* Motus      - https://motu-tool.org/
* Metaphlan  - http://huttenhower.sph.harvard.edu/metaphlan

Visualisierung der Ergebnisse:
* Krona      - https://github.com/marbl/Krona/wiki

## Methods

### Krakenuniq

Das Profiling basierend auf der  whole genome sequence erfolgte mit Krakenuniq.

```sh
krakenuniq --db $DBDIR --threads 10 --report-file kraken_taxonomy_profile.txt --paired read1.fq read2.fq  > kraken_read_classification.tsv
```

### Metaxa

Das Profiling basierend auf den rRNA Genen erfolgte mit Metaxa.

```sh
metaxa2 -1 read1.fq -2 read2.fq -g ssu --mode metagenome --plus T --cpu 8 --megablast T -o $OUT
```

###  Motus

Das Profiling basierend auf universal single-copy marker genes erfolgte mit motus.

```sh
motus profile -f read1.fq -r read2.fq -t 30 > $OUT
```

### Metaphlan

Das Profiling basierend auf clade-specific marker genes erfolgte mit Metaphlan.

```sh
metaphlan2.py read1.fq,read2.fq --bowtie2out metagenome.bowtie2.bz2 --nproc 8 --input_type fastq > $OUT
```

### Krona

Die Darstellung der Ergebnisse erfolgte mit Krona.

```sh
ktImportText -o $OUT.krona.html metaphlan_krona.out
```

## Results and discussion
Die Ergebnisse der einzelnen Tools sind nachfolgend angeführt:

<img src="https://github.com/fhwnmatt/Metagenomics/blob/master/figures/Metaxa.png" title="Metaxa">
<figcaption> Metaxa </figcaption>

-------------------------

<img src="https://github.com/fhwnmatt/Metagenomics/blob/master/figures/Metaphlan.png" title="Metaphlan">
<figcaption> Metaphlan </figcaption>

-------------------------

# Assembly

## Introduction

## Methods

Das Metagenom wurde mit dem Megahit assembler durchgeführt.

```sh
megahit -1 read1.fq -2 read2.fq -o megahit_out
```
Das Assembly wurde mit Quast evaluiert.

```sh
quast -1 read1.fq -2 read2.fq $IN
```

Die contig coverage wurde mit bbmap/bbwrap durchgeführt.

```sh
bbwrap.sh ref=$IN in=read1.fq in2=read2.fq out=aln.sam.gz kfilter=22 subfilter=15 maxindel=80
pileup.sh in=aln.sam.gz out=cov.txt
```

## Results and discussion


# Binning

## Introduction

## Methods

Das Binning erfolgte mit MaxBin2 und MetaBAT.

```sh
run_MaxBin.pl -contig $IN -out $OUT -reads read1.fq -reads2 read2.fq -thread 8 -markerset 40
metabat -i $IN -a depth.txt -o $OUT -v --saveTNF saved.tnf --saveDistance saved.dist
```

Mit Checkm wurden die Ergebnisse evaluiert.

```sh
checkm lineage_wf -t 20 -x fa $IN $OUT > checkm_summary.txt
```

## Results and discussion

# References

Breitwieser, F. P., Lu, J., & Salzberg, S. L. (2017). A review of methods and databases for metagenomic classification and assembly. Briefings in Bioinformatics, (June), 1–15. https://doi.org/10.1093/bib/bbx120

