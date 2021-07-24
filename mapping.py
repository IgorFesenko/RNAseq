#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

dir_name = r'~/trimmed'

names = [n for n in os.listdir(dir_name)]

corr_name = set([i.rsplit('_', maxsplit=1)[0] for i in names])

#print(*corr_name)

for i in corr_name:
    R1 = f"{i}_R1.fastq.gz"
    R2 = f"{i}_R2.fastq.gz"
    print(i)
    os.system(f"hisat2 --no-softclip --summary-file ./bam/{i}.log -x ./genome/Homo_sapiens.GRCh38.dna_rm -1 ./trimmed/{R1} -2 ./trimmed/{R2} | samtools view -Sb - > ./bam/{i}.bam")
    os.system(f"samtools sort -o ./bam/{i}.bam  ./bam/{i}.bam")
    os.system(f"samtools index ./bam/{i}.bam")
    print(f"DONE: {i}")
    print("____________")