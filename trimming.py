#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

dir_name = r'/home/admin_moss/fq_data'

names = [n for n in os.listdir(dir_name)]

corr_name = set([i.rsplit('_', maxsplit=1)[0] for i in names])

#print(*corr_name)

for i in corr_name:
    R1 = f"{i}_R1.fastq.gz"
    R2 = f"{i}_R2.fastq.gz"
    os.system(f"java -jar trimmomatic-0.39.jar PE {dir_name}/{R1} {dir_name}/{R2} ~/trimmed/{R1} ~/trimmed/un_{R1} ~/trimmed/{R2} ~/trimmed/un_{R2} ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:true TRAILING:3 MINLEN:36")