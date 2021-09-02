# Index Genomes (Nextflow)
Generate a comprehensive set of indexes for a single reference genome

![index-genomes-nf](https://github.com/FredHutch/index-genomes-nf/actions/workflows/test.yaml/badge.svg)

## What is a genome index?

For many different applications in bioinformatic analysis of genome-associated
datasets (whole-genome sequencing, transcriptional profiling, single-cell
analysis, etc.), specialized algorithms are employed to find the potential source
of individual small genome fragments from a much larger genome reference. Depending
on the approach of different algorithms, that reference genome must be pre-processed
to create an 'index' -- a compressed data structure which can be more efficiently
searched by the algorithm of interest.

## What does this workflow do?

Given a single genome reference sequence in FASTA format, this workflow will
generate the indexes for a variety of alignment algorithms listed below.
Each index will be saved to a separate sub-folder within the output directory,
named for the tool which was used to generate it. Any logs produced by the
indexing process will also be written to that folder.

## What do I need to specify?

In addition to the (1) reference genome in FASTA format and (2) the directory
to which all outputs should be written, you may optionally specify:

3. The amount of CPU and memory to allot each process
4. The compute resources to use (local computer, SLURM, cloud, etc.)

## Included tools

The tools which are currently supported by this workflow are:

- Bismark
- BLAST
- Bowtie2
- Bowtie
- BWA
- STAR
- STAR2

## How do I run the workflow?

The workflow is written using the
[Nextflow](https://www.nextflow.io/docs/latest/index.html) workflow management
system. For Fred Hutch researchers, documentation on how to get started using
Nextflow appropriately can be found on the
[Fred Hutch SciWiki](https://sciwiki.fredhutch.org/hdc/hdc_workflows/). For
other users, please consult the excellent
[Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) to
familiarize yourself with installing Nextflow.

Once you have Nextflow set up appropriately, this workflow can be run with
a command like:

```#!/bin/bash
nextflow run FredHutch/index-genomes-nf --input input.fasta.gz --output output_folder --gtf input.gtf.gz
```

For a complete set of documentation on workflow usage, please run:

```#!/bin/bash
nextflow run FredHutch/index-genomes-nf --help
```

## Specifying compute resources

To set the amount of CPUs and memory which will be used by each process,
use the flags for `--cpus` (default: 1) and `--ram_gb` (default: 4). The number
of CPUs allotted will also be passed to the algorithm being run to enable
multithreading when available.

To run this workflow on a high-performance computing system, use the appropriate
[Nextflow configuration](https://www.nextflow.io/docs/latest/config.html#config-scopes)
at runtime. Alternatively, users at Fred Hutch can follow
[institution-specific documentation](https://sciwiki.fredhutch.org/hdc/workflows/running/on_gizmo/).

## How do we know what version is being run?

One of the hardest things in bioinformatics is knowing exactly which version
of a software tool is going to be run. For some of these alignment algorithms,
the differences between software versions can be very important. To help provide
some control over these software versions, the workflow is set up to run each
tool using a publicly available software 'container' using either Docker or
Singularity. However, it is very important to note that the user is responsible
for setting up Nextflow to use the appropriate tool, either
[Docker](https://www.nextflow.io/docs/latest/docker.html) or 
[Singularity](https://www.nextflow.io/docs/latest/singularity.html). We took
this approach because users often operate under very different restrictions.
For example, users of institutional SLURM clusters may be able to use Singularity
while being disallowed Docker (because it requires root access). Users of cloud
computing resources will use Docker images natively as part of configuring those
resources. In other words, we hope that the user will not find it too difficult
to get Nextflow configured appropriately for their considerations, and then this
(and all) Nextflow workflows can run as desired.

## Using specific software versions

The recommended approach for specifying a non-default version of a software
package is to use the flags --tool_repository and --tool_tag. For a given
tool (e.g. `bowtie2`), the software image that will be used can be customized
following the pattern `"${tool_repository}:${tool_tag}"`. Most of the tools
used here have been provided by the community-driven
[BioContainers](https://biocontainers.pro/) resource. Some examples of repositories
which have a range of versions to select from can be seen here for
[bwa](https://quay.io/repository/biocontainers/bwa?tab=tags),
[bowtie2](https://quay.io/repository/biocontainers/bowtie2?tab=tags), and
[blast](https://quay.io/repository/biocontainers/blast?tab=tags).

## Gene annotations

The STAR tool requires that a set of gene annotations will be provided in
GTF format. This file can be specified for the workflow with the `--gtf`
flag.

## Getting help

This workflow has been written by the Fred Hutch Data Core. For help with
getting started or requests for additional tools to be added, please get
in touch by email at hutchdatacore (Fred Hutch).
