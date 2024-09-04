<p align="center">
    <a href="https://github.com/maxibor/floria-strainer/actions/workflows/ci.yaml/badge.svg"><img src="https://github.com/maxibor/floria-strainer/actions/workflows/ci.yaml/badge.svg"></a>
    <br>
   <img src="https://github.com/maxibor/floria-strainer/raw/master/assets/img/floria_strainer_logo.png" width="400">
</p>

---

**floria-strainer**: strains out genomes from [floria](https://github.com/bluenote-1577/floria).

## Introduction

Given the output of the strain haplotyping software [floria](https://github.com/bluenote-1577/floria)[^1] , **floria-strainer** computes the allele frequency at each variable position of each haploset identified by floria to cluster them into the different strains composing the mixture, using a Gaussian Mixture Model. 

## Install

```bash
pip install floria-strainer
```

## Quick start

Running floria-strainer on the provided test data

```bash
$ floria-strainer -b tests/data/test_short.bam tests/data/floria_out_dir
INFO - Writing the straining informations to floria_strained.strained.csv.
INFO - 2 strains were found: 0, 1
INFO - Writing the BAM file in tag mode to floria_strained.bam.
```

## Running floria and floria-strainer on your date

### 1. Variant calling 

The floria authors [recommend](https://phase-doc.readthedocs.io/en/latest/how-to-guides/htg3.html#short-reads) [freebayes](https://github.com/freebayes/freebayes)[^2] to call the variants from your alignment files. In this step, Freebayes is used in "naive" mode, just to count all the "alleles" present at each variable position, using the `--pooled-continous` option.

```bash
freebayes \
	-f reference.fa \
	--pooled-continuous \
	alignment.bam > freevayes_output.vcf
```

### 2. Running floria

```bash
floria \
	-b alignment.bam \
	-v freevayes_output.vcf \
	-r reference.fa \
```

The floria output will be written by default to the `floria_out_dir` directory

### 3. Running floria-strainer 


#### floria-strainer `tag` mode

```bash
floria-strainer \
	--bam alignment.bam \
	-m tag \
	floria_out_dir 
```

floria-stainer in the `tag` mode will produce a bam file (by default, `floria_strained.bam`) which contains the strain-phasing information in the `ST` tag, and the haploset information in the `HP` tag.

#### floria-strainer `split` mode

```bash
floria-strainer \
	--bam alignment.bam \
	-m split \
	floria_out_dir 
```

floria-stainer in the `split` mode will produce one bam file per detected strain (by default `floria_strained.{strain_number}.bam`). Reads not being part of clustered to a specific strain will be written to all bam files.

## Visualizing the strain clustering method


- In this example, short-reads are aligned to a reference genome (fig 1). At each variable position, we can observe a 2 different alleles, which in the case of this haploid organism, corresponds to a mixture of 2 different strains.

- floria-strainer takes the output of floria, with reads having been assigned to a haploset (fig 2) based on the MEC criteria. <details>
  <summary>Expand to see a high level overview of floria</summary>
  
  The haplotyping process floria is conceptually similar to the one of *de novo* assembly of reads into contigs, where instead of contigs, read haplosets are the results of the flow (ie. least costly path) optimization through the different reads represented as a graph. More details in the floria article [^1].
</details>
  

- Based on the average allele frequency of each haploset, reads are clustered in the different strains. In this example (fig 3), there are two strains of minor and major allele frequency, which floria-strainer clustered in strain 0 and strain 1.

Reads that weren't assigned to any haploset by floria, or whose haploset do not cluster well enough are not assigned to any strain. They are considered to be shared by the different strains present in the alignment.

<img src="https://github.com/maxibor/floria-strainer/raw/master/assets/img/igv_no_tag.png" width=70%>  

Fig 1: Reads aligned to the reference genome, visualized in IGV. The top track represents the reference genome, with variants indicated in the different colors.

<img src="https://github.com/maxibor/floria-strainer/raw/master/assets/img/igv_hp_tag.png" width=70%>  

Fig 2: Reads are grouped and colored by the `HP` **HaPloset** tag as annotated by Floria.

<img src="https://github.com/maxibor/floria-strainer/raw/master/assets/img/igv_st_tag.png" width=70%>  

Fig 3: Reads are grouped and colored by the `ST` **STrain** tag as annotated by floria-strainer.

## Working with Ancient DNA (aDNA)

In the case of very short DNA fragments typically observed in aDNA samples, the following parameter set have been found to work well. Tweak them for your own needs, and experiment as needed:

#### For freebayes variant-calling

```bash
freebayes \
	-f reference.fa \
	--min-base-quality 20 \
	--min-mapping-quality 30 \
	--min-coverage 6 \
	--limit-coverage 30 \
	--pooled-continuous \
	--use-best-n-alleles 3 \
  alignment.bam > freebayes_output.vcf
```

#### For floria

```bash
floria \
	-b alignment.bam \
  -r reference.fa \
	-v freebayes_output.vcf \
	--snp-density 0.00001 \
	--snp-count-filter 10 \
	--ploidy-sensitivity 3 \
	--threads 20
```

#### For floria-strainer (in `tag` or `split` mode)

```bash
floria-strainer \
	--bam alignment.bam \
	--hapq-cut 0 \
	floria_out_dir
```

## Help

```bash
$floria-strainer --help
                                                                                                                                                                            
 Usage: floria-strainer [OPTIONS] FLORIA_OUTDIR                                                                                                                             
                                                                                                                                                                            
 Strain the haplotypes in the floria output directory.                                                                                                                      
 Author: Maxime Borry                                                                                                                                                       
                                                                                                                                                                            
╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│    --version                      Show the version and exit.                                                                                                             │
│    --nb-strains  -n  INTEGER      Number of strains to keep. If 0, the number of strains will be determined by the mean floria average strain count with HAPQ > 15.      │
│                                   [default: 0]                                                                                                                           │
│    --hapq-cut    -h  INTEGER      Minimum HAPQ threshold [default: 15]                                                                                                   │
│    --sp-cut      -s  FLOAT        Minimum strain clustering probability threshold [default: 0.5]                                                                         │
│ *  --bam         -b  PATH         Input BAM file [required]                                                                                                              │
│    --mode        -m  [tag|split]  BAM output mode. Tag: add ST (strain) tags to the reads. Split: split the reads in different BAM files per strain. [default: tag]      │
│ *  --basename    -o  TEXT         Output file basaneme [default: floria_strained] [required]                                                                             │
│    --help                         Show this message and exit.                                                                                                            │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

## Tests

After cloning/pulling repo

```bash
$ pip install poetry pytest
$ poetry run pytest -vv
```

[^1]: [Floria: fast and accurate strain haplotyping in metagenomes](https://doi.org/10.1093/bioinformatics/btae252) 
[^2]: [Haplotype-based variant detection from short-read sequencing](https://doi.org/10.48550/arXiv.1207.3907)