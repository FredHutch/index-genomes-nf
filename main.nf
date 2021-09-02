#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.input = false
params.gtf = false
params.output = false
params.cpus = 1
params.ram_gb = 4

// Containers used to run each task
params.bismark_repository = "quay.io/biocontainers/bismark"
params.bismark_tag = "0.23.1--hdfd78af_0"

params.blast_repository = "quay.io/biocontainers/blast"
params.blast_tag = "2.12.0--pl5262h3289130_0"

params.bowtie2_repository = "quay.io/biocontainers/bowtie2"
params.bowtie2_tag = "2.4.4--py37h13ad519_0"

params.bowtie_repository = "quay.io/biocontainers/bowtie"
params.bowtie_tag = "1.3.0--py38hcf49a77_2"

params.bwa_repository = "quay.io/biocontainers/bwa"
params.bwa_tag = "0.7.17--h5bf99c6_8"

params.star_repository = "quay.io/biocontainers/star"
params.star_tag = "2.7.9a--h9ee0642_0"

params.star2_repository = "quay.io/biocontainers/star"
params.star2_tag = "2.7.9a--h9ee0642_0"



// Import the processes to run in this workflow
include {
    bismark;
    blast;
    bowtie2;
    bowtie;
    bwa;
    star;
    star2;
} from './modules' params(
    cpus: params.cpus,
    ram_gb: params.ram_gb,
    output: params.output,
    container__bismark: "${params.bismark_repository}:${params.bismark_tag}",
    container__blast: "${params.blast_repository}:${params.blast_tag}",
    container__bowtie2: "${params.bowtie2_repository}:${params.bowtie2_tag}",
    container__bowtie: "${params.bowtie_repository}:${params.bowtie_tag}",
    container__bwa: "${params.bwa_repository}:${params.bwa_tag}",
    container__star: "${params.star_repository}:${params.star_tag}",
    container__star2: "${params.star2_repository}:${params.star2_tag}",
)

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/index-genomes-nf <ARGUMENTS>

    Required Arguments:
      --input               Genome sequence in FASTA format
      --output              Base directory for writing outputs

    Optional Arguments:
      --gtf                 Path to gene annotations in GTF format (required for STAR)
      --cpus                Number of CPUs to use for each process (default: 1)
      --ram_gb              Amount of RAM (GB) to use for each process (default: 4)

    """.stripIndent()
}


workflow {

    // Show help message if the user specifies the --help flag at runtime
    if (params.help){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    // The user must specify each of the required arguments
    if (!params.output || !params.input){
        log.info"""

        -----------------------
        MISSING REQUIRED INPUTS
        -----------------------

        """.stripIndent()
        helpMessage()

        // Exit out and do not run anything else
        exit 0
    }

    // Make a channel with the input genome
    input_genome = Channel.fromPath(params.input)

    // Run each tool
    bismark(input_genome)
    blast(input_genome)
    bowtie2(input_genome)
    bowtie(input_genome)
    bwa(input_genome)

    // If a GTF annotation was provided
    if ( params.gtf ) {

        // Make a channel with the input GTF annotations
        input_gtf = Channel.fromPath(params.gtf)

        star(input_genome, input_gtf)
        star2(input_genome, input_gtf)
    }

}
