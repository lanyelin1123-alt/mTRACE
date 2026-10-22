#!/bin/bash
wd=/path/to/project_work_folder
mkdir -p $wd/result/bam/stranded_excldWT
cd $wd/result/bam/exclude_wt
stranded_dir=$wd/result/bam/stranded_excldWT

##split stranded bam files
for file in *.bam
do
{
bn=${file%%.bam*}
echo processing $bn
samtools view -@ 16 -b -f 99  $file | samtools sort -o $stranded_dir/${bn}_99.bam -@ 16 
samtools view -@ 16 -b -f 147 $file | samtools sort -o $stranded_dir/${bn}_147.bam -@ 16
samtools view -@ 16 -b -f 163 $file | samtools sort -o $stranded_dir/${bn}_163.bam -@ 16
samtools view -@ 16 -b -f 83  $file | samtools sort -o $stranded_dir/${bn}_83.bam -@ 16

echo merging tmp BAM files
samtools merge -@ 16 $stranded_dir/${bn}_pos_sort.bam $stranded_dir/${bn}_99.bam $stranded_dir/${bn}_147.bam 
samtools merge -@ 16 $stranded_dir/${bn}_neg_sort.bam $stranded_dir/${bn}_83.bam $stranded_dir/${bn}_163.bam 
samtools index -@ 16 $stranded_dir/${bn}_pos_sort.bam
samtools index -@ 16 $stranded_dir/${bn}_neg_sort.bam

echo cleaning up
rm $stranded_dir/${bn}_99.bam
rm $stranded_dir/${bn}_147.bam
rm $stranded_dir/${bn}_163.bam
rm $stranded_dir/${bn}_83.bam
}
done

##switch R2 flag to R1 flag
cd $stranded_dir
mkdir -p switch
iswitch_dir=switch

for neg_file in *neg_sort.bam
do
{
###neg_strand:R2toR1(163to99)
samtools view -h ${neg_file} -@ 2| awk 'BEGIN{FS=OFS="\t"}{if ($2==163) $2=99}{if ($2==147) $2=83}1' | samtools view -bS -@ 2 > $switch_dir/${neg_file%%.bam*}_switch.bam 
}&
done
wait

for pos_file in *pos_sort.bam
do
{
###pos_strand:R2toR1(147to83)
samtools view -h ${pos_file} -@ 2| awk 'BEGIN{FS=OFS="\t"}{if ($2==163) $2=99}{if ($2==147) $2=83}1' | samtools view -bS -@ 2 > $switch_dir/${pos_file%%.bam*}_switch.bam 
} & 
done
wait
echo "finish"
