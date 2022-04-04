#!/usr/bin/bash
# -*- coding: utf-8 -*-



for i in ~/fq_data/*_1.fastq;

do
   base=$(basename ${i} _1.fastq)
   echo $base
   
   java -jar trimmomatic-0.39.jar PE ~/fq_data/${i} ~/fq_data/${base}_2.fastq \
                ~/trimmed/${base}_1.fastq ~/trimmed/${base}_1un.trim.fastq \
                ~/trimmed/${base}_2.fastq ~/trimmed/${base}_2un.trim.fastq \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
done
