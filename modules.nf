
process bismark {
    container "${params.container__bismark}"
    cpus params.cpus
    memory "${params.ram_gb} GB"
    publishDir "${params.output}/bismark/", mode: "copy", overwrite: true

    input:
        file "input.fasta.gz"

    output:
        file "*"

"""#!/bin/bash

set -e

ls -lahtr

bismark_genome_preparation \
    --verbose ./ \
    2>&1 > bismark.index.log
"""
}

process blast {
    container "${params.container__blast}"
    cpus params.cpus
    memory "${params.ram_gb} GB"
    publishDir "${params.output}/blast/", mode: "copy", overwrite: true

    input:
        file "input.fasta.gz"

    output:
        file "*"

"""#!/bin/bash

set -e

# Decompress the input
gunzip -c input.fasta.gz > input.fasta

makeblastdb \
    -in input.fasta \
    -dbtype nucl \
    2>&1 > BLAST.index.log

rm input.fasta
"""
}

process bowtie2 {
    container "${params.container__bowtie2}"
    cpus params.cpus
    memory "${params.ram_gb} GB"
    publishDir "${params.output}/bowtie2/", mode: "copy", overwrite: true

    input:
        file "input.fasta.gz"

    output:
        file "*"

"""#!/bin/bash

set -e

bowtie2-build \
    input.fasta.gz \
    index \
    2>&1 > bowtie2.index.log
"""
}

process bowtie {
    container "${params.container__bowtie}"
    cpus params.cpus
    memory "${params.ram_gb} GB"
    publishDir "${params.output}/bowtie/", mode: "copy", overwrite: true

    input:
        file genome_fasta

    output:
        file "*"

"""#!/bin/bash

set -e

bowtie-build \
    "${genome_fasta}" \
    index \
    2>&1 > bowtie.index.log
"""
}

process bwa {
    container "${params.container__bwa}"
    cpus params.cpus
    memory "${params.ram_gb} GB"
    publishDir "${params.output}/bwa/", mode: "copy", overwrite: true

    input:
        file genome_fasta

    output:
        file "*"

"""#!/bin/bash

set -e

bwa index \
    "${genome_fasta}" \
    2>&1 > BWA.index.log
"""
}

process star {
    container "${params.container__star}"
    cpus params.cpus
    memory "${params.ram_gb} GB"
    publishDir "${params.output}/star/", mode: "copy", overwrite: true

    input:
        file "input.fasta.gz"
        file "input.gtf.gz"

    output:
        file "*"

"""#!/bin/bash

set -e

# Decompress the reference
gunzip -c input.fasta.gz > input.fasta
gunzip -c input.gtf.gz > input.gtf

mkdir output

STAR \
    --runMode genomeGenerate \
    --runThreadN ${params.cpus} \
    --genomeDir ./ \
    --genomeFastaFiles input.fasta \
    --sjdbGTFfile input.gtf \
    2>&1 > STAR.index.log

# Clean up the temporary decompressed files
rm input.fasta
rm input.gtf

"""
}

process star2 {
    container "${params.container__star2}"
    cpus params.cpus
    memory "${params.ram_gb} GB"
    publishDir "${params.output}/star2/", mode: "copy", overwrite: true

    input:
        file "input.fasta.gz"
        file "input.gtf.gz"

    output:
        file "*"

"""#!/bin/bash

set -e

# Decompress the reference
gunzip -c input.fasta.gz > input.fasta
gunzip -c input.gtf.gz > input.gtf

mkdir output

STAR \
    --runMode genomeGenerate \
    --runThreadN ${params.cpus} \
    --genomeDir ./ \
    --genomeFastaFiles input.fasta \
    --sjdbGTFfile input.gtf \
    2>&1 > STAR.index.log

# Clean up the temporary decompressed files
rm input.fasta
rm input.gtf

"""
}
