wd=wd=/path/to/project_work_folder
highcon_dir=${wd}/result/peak/macs2_excldWT/highcon
idr_dir=${highcon_dir}/idr
mkdir -p $idr_dir

sample_group=(prefix of sample replicate1 peak file)
##e.g. sample_group=(ecog2_dox_rep1 rbfox2_ecog2_dox_rep1 wt_rep1)

for file in ${sample_group[@]}
do
{
prefix=${file}
###Sort peak by -log10(p-value)
sort -k8,8nr $highcon_dir/${prefix}_ip.excldWT_highcon_ss.narrowPeak > $idr_dir/${prefix}.sortp.narrowPeak
sort -k8,8nr $highcon_dir/${prefix/rep1/rep2}_ip.excldWT_highcon_ss.narrowPeak > $idr_dir/${prefix/rep1/rep2}.sortp.narrowPeak
###idr reproducible peak
/public1/work/yelin/software/idr-2.0.2/bin/idr --samples $idr_dir/${prefix}.sortp.narrowPeak $idr_dir/${prefix/rep1/rep2}.sortp.narrowPeak --input-file-type narrowPeak --rank p.value --output-file $idr_dir/${prefix/_rep1/}.idrPeak.bed --log-output-file $idr_dir/${prefix/_rep1/}.idr.log
sort -k1,1 -k2,2n $idr_dir/${prefix/_rep1/}.idrPeak.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' > $idr_dir/${prefix/_rep1/}.idr.bed
sort -k1,1 -k2,2n $idr_dir/${prefix/_rep1/}.idrPeak.bed > $idr_dir/${prefix/_rep1/}.idrPeak.sort.bed
}
done

idr_filter_dir=${highcon_dir}/idr_filter
mkdir -p $idr_filter_dir

sample_group=(prefix of sample names)
##e.g. sample_group=(ecog2_dox rbfox2_ecog2_dox wt)
for file in ${sample_group[@]}
do
{
###filter by peak signal-value (fold-enrichment)
awk '$7>=6' $highcon_dir/${file}_rep1_ip.excldWT_highcon_ss.narrowPeak > $highcon_dir/${file}_rep1_filter6.narrowPeak
awk '$7>=6' $highcon_dir/${file}_rep2_ip.excldWT_highcon_ss.narrowPeak > $highcon_dir/${file}_rep2_filter6.narrowPeak
###sort peak by -log10(p-value)
sort -k8,8nr $highcon_dir/${file}_rep1_filter6.narrowPeak > $idr_filter_dir/${file}_rep1_filter6.sortp.narrowPeak
sort -k8,8nr $highcon_dir/${file}_rep2_filter6.narrowPeak > $idr_filter_dir/${file}_rep2_filter6.sortp.narrowPeak
###idr reproducible filtered-peak
/public1/work/yelin/software/idr-2.0.2/bin/idr --samples $idr_filter_dir/${file}_rep1_filter6.sortp.narrowPeak $idr_filter_dir/${file}_rep2_filter6.sortp.narrowPeak --input-file-type narrowPeak --rank p.value --output-file $idr_filter_dir/${file}_filter6.idrPeak.bed --log-output-file $idr_filter_dir/${file}_filter6.idr.log
###sort peak by genome position
sort -k1,1 -k2,2n $idr_filter_dir/${file}_filter6.idrPeak.bed > $idr_filter_dir/${file}_filter6.idrPeak.sort.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' $idr_filter_dir/${file}_filter6.idrPeak.sort.bed > $idr_filter_dir/${file}_filter6.idr.bed
}
done

