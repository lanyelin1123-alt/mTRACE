wd=/path/to/project_work_folder
highcon_peak_dir=$wd/result/peak/macs2/highcon
bam_dir=$wd/result/bam/dedup
exclude_bam_dir=$wd/result/bam/exclude_wt

mkdir -p $exclude_bam_dir  $exclude_bam_dir/removed

wt_region_bed=$highcon_peak_dir/wt_replicate_merge.bed

cd $bam_dir
for file in *.bam
do
{
bn=${file%%.bam}
samtools view -@ 24 -b -L $wt_region_bed -U $exclude_bam_dir/${bn}.excldWT.bam $bam_dir/${file} -o $exclude_bam_dir/removed/${bn}.removed.bam
samtools index -@ 24 $exclude_bam_dir/${bn}.excldWT.bam
}
done
