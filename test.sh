#!/bin/bash

set -e

NXF_VER=21.04.1 \
nextflow \
    run \
        main.nf \
        -w work/ \
        -with-docker ubuntu:latest \
        --input test_data/GCF_000819615.1_ViralProj14015_genomic.fna.gz \
        --gtf test_data/GCF_000819615.1_ViralProj14015_genomic.gtf.gz \
        --output test_output \
        -with-report test_output/report.html \
        -resume
