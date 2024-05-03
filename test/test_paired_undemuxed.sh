#!/usr/bin/env bash

mkdir -p paired_test_undemux
cd paired_test_undemux

../../eDNAFlow.nf \
  --paired \
  --ignore-blast-env \
  --reads ../fastq/paired_undemux \
  --barcode '../pe_bc_*' \
  --blast-db ../blastdb/paired_undemuxed \
  --collapse-taxonomy
