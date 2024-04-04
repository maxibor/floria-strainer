<p align="center">
    <a href="https://github.com/maxibor/floria-strainer/actions/workflows/ci.yaml/badge.svg"><img src="https://github.com/maxibor/floria-strainer/actions/workflows/ci.yaml/badge.svg"></a>
   <img src="assets/img/floria_strainer_logo.png" width="400">
</p>

---

**floria-strainer**: strains out genomes from [floria](https://github.com/bluenote-1577/floria).

## Introduction

Given the output of the Strain Haplotyping software [floria](https://github.com/bluenote-1577/floria)[^1] , **floria-strainer** computes the allele frequency at each variable position of each haploset identified by floria to cluster them into the different mixtures of strains using a Gaussian Mixture Model.

[^1]: [Floria: Fast and accurate strain haplotyping in metagenomes](https://www.biorxiv.org/content/10.1101/2024.01.28.577669v1.full)  

## Install

TBD

## Help

```bash
floria-strainer --help

 Usage: floria-strainer [OPTIONS] FLORIA_OUTDIR

 Strain the haplotypes in the floria output directory.
 Author: Maxime Borry

╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│    --version                      Show the version and exit.                                                                    │
│    --nb-strains  -n  INTEGER      Number of strains to keep. If 0, the number of strains will be determined by the mean floria  │
│                                   average strain count with HAPQ > 15.                                                          │
│                                   [default: 0]                                                                                  │
│    --hapq-cut    -h  INTEGER      Minimum HAPQ threshold [default: 15]                                                          │
│    --sp-cut      -s  FLOAT        Minimum strain clustering probability threshold [default: 0.5]                                │
│ *  --bam         -b  PATH         Input BAM file [required]                                                                     │
│    --mode        -m  [tag|split]  BAM output mode. Tag: add ST (strain) tags to the reads. Split: split the reads in different  │
│                                   BAM files per strain.                                                                         │
│                                   [default: tag]                                                                                │
│ *  --basename    -o  TEXT         Output file basanem [default: floria_strained] [required]                                     │
│    --help                         Show this message and exit.                                                                   │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```