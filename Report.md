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

## Results and discussion

# Assembly

# Binning

# References

