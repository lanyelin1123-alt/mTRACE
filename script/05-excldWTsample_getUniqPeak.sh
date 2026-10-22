wd=/path/to/project_work_folder
uniq_peak_folder=${wd}/result/peak/macs2_excldWT/sample_uniqPeak_f50
highcon_peak_folder=${wd}/result/peak/macs2_excldWT/highcon

idr_filter_folder=$highcon_peak_folder/idr_filter
idr_folder=$highcon_peak_folder/idr
m6Acov_folder=$wd/result/m6Acov
mkdir -p $m6Acov_folder $uniq_peak_folder

delete_peak=(/path/to/control_sample.idr.bed files)
###e.g. delete_peak=($idr_folder/ecog2_dox.idr.bed $idr_folder/ecog2_con.idr.bed $idr_folder/rbfox2_ecog2_con.idr.bed $idr_folder/wt.idr.bed)
delete_idrPeak=(path/to/control_sample.idrPeak.sort.bed files)
###e.g. delete_idrPeak=($idr_folder/ecog2_dox.idrPeak.sort.bed $idr_folder/ecog2_con.idrPeak.sort.bed $idr_folder/rbfox2_ecog2_con.idrPeak.sort.bed $idr_folder/wt.idrPeak.sort.bed)
sample_group=(prefix of experimental sample)
###e.g. sample_group=(rbfox2_ecog2_dox)
for name in ${sample_group[@]}
do
{
bedtools intersect -a $idr_filter_folder/${name}_filter6.idr.bed -b ${delete_peak[@]} -f 0.5 -v -s -wa  | sort -k1,1 -k2,2n | uniq  > $uniq_peak_folder/${name}.idr.uniq.bed 
bedtools intersect -a $idr_filter_folder/${name}_filter6.idrPeak.sort.bed -b ${delete_idrPeak[@]} -f 0.5 -v -s -wa | sort -k1,1 -k2,2n | uniq  > $uniq_peak_folder/${name}.idrPeakSort.uniq.bed
awk 'BEGIN{OFS="\t"}{print $1,int(($13+$16+$17+$20)/2),int(($13+$16+$17+$20)/2)+1,"'${name}'",".",$6}' $uniq_peak_folder/${name}.idrPeakSort.uniq.bed > $uniq_peak_folder/${name}.idrPeakSort.uniq.summit.bed
}
done

###bedtools multicov get peak count
bam_dir=$wd/result/bam/exclude_wt
bam_list=(/path/to/control_and_experimental_ip_group_bam files)
###e.g. bam_list=($bam_dir/ecog2_dox_rep1_ip.excldWT.bam $bam_dir/ecog2_dox_rep2_ip.excldWT.bam $bam_dir/rbfox2_ecog2_dox_rep1_ip.excldWT.bam $bam_dir/rbfox2_ecog2_dox_rep2_ip.excldWT.bam)

bedtools multicov -bams ${bam_list[@]} -bed $uniq_peak_folder/experimental_sample.idr.uniq.bed -s > $m6Acov_folder/experimental_sample_uniqIdr_ip_multicov.txt 
awk 'BEGIN{OFS="\t"}{print $1"|"$2"|"$3"|"$4"|"$5"|"$6,$7,$8,$9,$10,$11,$12}' $m6Acov_folder/experimental_sample_uniqIdr_ip_multicov.txt > $m6Acov_folder/experimental_sample_uniqIdr_ip_multicov.std.txt


