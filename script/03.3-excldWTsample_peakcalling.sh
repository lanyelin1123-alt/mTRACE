#peak calling
#set -e 
wd=/path/to/project_work_folder
excldWT_bam_dir=$wd/result/bam/exclude_wt
excldWT_peak=$wd/result/peak/macs2_excldWT
mkdir -p $excldWT_peak $excldWT_peak/highcon
peak_dir=$wd/result/peak

cd $wd/result/bam/stranded_excldWT/switch
for file in *ip*bam
do
{
prefix=${file%%_sort_switch.bam}
###extsize 150; qvalue 0.05; fold_change 1;keepdup all
macs2 callpeak -t $file -c ${file/_ip/_input} --keep-dup all -n ${prefix} --outdir $peak_dir/macs2_excldWT -f BAM -g hs --nomodel --extsize 150 --fe-cutoff 1 
###filter high-confidence peaks (signalValue & -log10(qvalue))
awk '$7 >= 2 && $9 > 4' $peak_dir/macs2_excldWT/${prefix}_peaks.narrowPeak > $peak_dir/macs2_excldWT/highcon/${prefix}_highcon_peaks.narrowPeak
} 
done
echo "macs2 callpeak finish"

###add strand information to the callpeak results
cd $peak_dir/macs2_excldWT/highcon
for file in *neg_highcon_peaks.narrowPeak
do
{
awk 'BEGIN{OFS=FS="\t"}{if ($6==".") $6="-"}1' ${file} > ${file%%_neg*}_tmp_stranded.highcon.narrowpeak
awk 'BEGIN{OFS=FS="\t"}{if ($6==".") $6="+"}1' ${file/neg/pos} >> ${file%%_neg*}_tmp_stranded.highcon.narrowpeak
}
done

for i in *_tmp_stranded.highcon.narrowpeak
do
{
sort -k1,1 -k2,2n ${i} > ${i%%_tmp*}_highcon_ss.narrowPeak
awk 'BEGIN{OFS=FS="\t"}{print $1,$2+$10,$2+$10+1,$4,$5,$6}' ${i%%_tmp*}_highcon_ss.narrowPeak |sort -k1,1 -k2,2n > ${i%%_tmp*}_highcon_ss.summit.bed
}
done
echo "summit merge"

rm *_tmp_stranded.highcon.narrowpeak
