#!/usr/bin/env bash

mkdir -p single_test_demux
cd single_test_demux

../../eDNAFlow.nf \
  --single \
  --ignore-blast-env \
  --illumina-demultiplexed \
  --reads '../fastq/single_demuxed/*.fastq' \
  --barcode ../single_demuxed_barcode.tab \
  --blast-db ../blast/single_demuxed \
  --collapse-taxonomy
