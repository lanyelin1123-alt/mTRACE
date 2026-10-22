wd=/path/to/project_work_folder
db=$wd/split
index=/path/to/hisat2_index
fa=/path/to/genome_fasta

###rawdata qc & trim and cleandata qc
mkdir -p $wd/result $wd/result/rawdataqc $wd/result/cleandata $wd/result/cleandataqc

cd $db
for file in *_R1.fq.gz
do
{
fastqc -o ${wd}/result/rawdataqc -t 16 $file ${file/R1/R2}
cutadapt -a AGATCGGAAGAGC -A CTGTCTCTTATACAC  -m 15 -n 2 -j 16 -o ${wd}/result/cleandata/${file%%R1.fq.gz}R1.clean_umi.fq.gz -p ${wd}/result/cleandata/${file%%R1.fq.gz}R2.clean_umi.fq.gz $file ${file/R1/R2} 
fastp -i ${wd}/result/cleandata/${file%%R1.fq.gz}R1.clean_umi.fq.gz  -o $wd/result/cleandata/${file%%R1.fq.gz}R1.clean.fq.gz -I ${wd}/result/cleandata/${file%%R1.fq.gz}R2.clean_umi.fq.gz  -O $wd/result/cleandata/${file%%R1.fq.gz}R2.clean.fq.gz -t 6 -T 10 -w 16 -A -Q -l 15 -U --umi_loc=read1 --umi_len=10 --umi_prefix UMI 
fastqc -o ${wd}/result/cleandataqc -t 16 $wd/result/cleandata/${file%%R1.fq.gz}R1.clean.fq.gz $wd/result/cleandata/${file%%R1.fq.gz}R2.clean.fq.gz
} 
done
multiqc ${wd}/result/rawdataqc -o $wd/result/multiqc -n raw_data
multiqc ${wd}/result/cleandataqc -o $wd/result/multiqc -n clean_data
echo "Finish data clean"

###align clean reads to hg38 genome
mkdir -p $wd/result/bam $wd/result/bam/dedup $wd/result/bw
cd $wd/result/cleandata
for file in *R1.clean.fq.gz
do
{
echo ${file%%_R1.clean.fq.gz}
##alignment
hisat2 -p 16 -x $index -1 $file -2 ${file%%R1.clean.fq.gz}R2.clean.fq.gz --rna-strandness FR --trim5 2 --trim3 2 | samtools view -@ 16 -bS -f 2 -F 256 | samtools sort -@ 16 -o ${wd}/result/bam/${file%%_R1.clean.fq.gz}.bam
##PCR deduplication
gencore -i ${wd}/result/bam/${file%%_R1.clean.fq.gz}.bam -r $fa -u UMI --moderate_qual 20 -h $wd/result/bam/dedup/${file%%_R1.clean.fq.gz}_dedup.html | samtools sort -@ 16 -o $wd/result/bam/dedup/${file%%_R1.clean.fq.gz}.bam
samtools index -@ 16 $wd/result/bam/dedup/${file%%_R1.clean.fq.gz}.bam
##bam2bw
bamCoverage -b $wd/result/bam/dedup/${file%%_R1.clean.fq.gz}.bam -p 16 --binSize 15 -o $wd/result/bw/${file%%_R1.clean.fq.gz}.bw
} 
done
echo "finish"

