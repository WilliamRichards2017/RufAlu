#!/bin/bash

cat /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/input/input.txt | while read line
do
    f=$(echo $line | cut -d " " -f 3)
    ./../src/aluDetect /scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/"${f}".bam /scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/"${f}".bam.generator.V2.overlap.hashcount.fastq  /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/aluList/primate_non-LTR_Retrotransposon.fasta /uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa &

done
 

