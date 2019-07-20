# Metagenomik

<img src="https://github.com/fhwnmatt/Metagenomics/blob/master/figures/Abundance.png" title="Abundance">
<figcaption> Abundance </figcaption>

# Abstract

# Introduction

# Data generation

## Introduction

## Methods

Mit dem folgenden Befehl wurden 2 000 000 MiSeq reads von 2 bakteriellen und 1 Archae Genom - von RefSeq zufällig ausgewählt - simuliert.

```sh
iss generate --ncbi bacteria archaea --n_genomes_ncbi 2 1 --n_reads 2000000 --model MiSeq --output ncbi_lc --cpus 8
```

## Results and discussion

# Quality Control

## Introduction

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

# Profiling of the community

## Introduction

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

# Assembly

# Binning

# References

Breitwieser, F. P., Lu, J., & Salzberg, S. L. (2017). A review of methods and databases for metagenomic classification and assembly. Briefings in Bioinformatics, (June), 1–15. https://doi.org/10.1093/bib/bbx120

