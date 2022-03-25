#!/bin/bash


# prepearing

cd fq_data 
scp ./*.gz ****@******:./fq_data # копируем файлы на рабочий сервер

pip3 install multiqc # устанавливаем multiqc
mkdir bin # создаем папку для хранения нужных программ 
mkdir genome #папка для хранения генома и индексов
mkdir fq_data # для fastaqc
mkdir trimmed # для триммированных fastaqc
mkdir bam #хранение bam файлов
mkdir fastqc_res # файлы fastqc
mkdir fastqc_trimmed


cd bin # переходим в bin

# скачиваем и разархивируем Trimmomatic 
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

# скачиваем и разархивируем fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip #скачиваем fastqc
unzip fastqc_v0.11.9.zip # разархивируем его
chmod +x ~/bin/FastQC/fastqc # делаем файл fastqc исполняемым (чтобы его можно было запустить как программу)

#скачиваем hisat2
wget -O hisat2.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download # скачиваем hisat2
unzip hisat2.zip # разархивирeем его


cd ~ # выходим в домашнюю директорию



fq=~/fq_data


fastqc -o fastqc_res $fq/* # запускаем fastqc
multiqc ./fastqc_res -o ./fastqc_res # объединяем отчеты при помощи multiqc

#триммируем данные - запускаем скрипт
cd bin/Trimmomatic-0.39
nohup python3 trimming.py &

# оцениваем качество
cd ~ # возвращаемся в родительский каталог
fastqc -o fastqc_trimmed ~trimmed/* # запускаем fastqc
multiqc ./fastqc_trimmed -o ./fastqc_trimmed # объединяем отчеты при помощи multiqc

# скачиваем геном и аннотацию
cd genome
wget -c http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -c http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz # разархивируем
gunzip Homo_sapiens.GRCh38.104.gtf.gz # разархивируем

# индексируем, hisat2 установлен на сервере
hisat2-build ~/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./Homo_sapiens.GRCh38.dna &

# запускаем картирование
cd ~
nohup python3 mapping.py &

# запускаем featureCounts
featureCounts -T 20 -a ./genome/Homo_sapiens.GRCh38.104.gtf -o read_counts ./bam/*.bam

# копируем файлы на рабочий компьютер для анализа в R-studio
