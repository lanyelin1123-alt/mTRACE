###annotate m6A sites to region
wd=/path/to/project_work_folder
hg38_features=/path/to/reference/hg38_features_IntergenicExcldTSS500bp.bed
hg38_gene=/path/to/reference/hg38_gene.bed

annot2feature_dir=$wd/result/annot2feature
annot2gene_dir=$wd/result/annot2gene
filter_dir=$wd/result/excldPeak_filtered_by_limma
uniq_peak_dir=$wd/result/peak/macs2_excldWT/sample_uniqPeak_f50

mkdir -p $annot2feature_dir $annot2gene_dir $filter_dir

sample_group=(prefix of experimental sample)
###e.g. sample_group=(rbfox2_ecog2_dox)
for name in ${sample_group[@]} 
do
{
###extract peak summits (filtered by limma)
i=log2FCquant75
sort -k1,1 -k2,2n $filter_dir/${name}_${i}.bed > $filter_dir/${name}_${i}.sort.bed
bedtools intersect -a $uniq_peak_dir/${name}.idrPeakSort.uniq.summit.bed -b $filter_dir/${name}_${i}.sort.bed -s -wa | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"summit_"NR,$5,$6}' > $filter_dir/${name}_${i}.summit.bed

###anno2geneFeature
bedtools intersect -a $filter_dir/${name}_${i}.summit.bed  -b $hg38_features -s -wa -wb > $annot2feature_dir/${name}_${i}.genic.anno.tmp.bed
awk 'BEGIN{OFS="\t"}ARGIND==1{a[$4]=$0}ARGIND==2{if(!($4 in a)) {print $0,"Intergenic","Intergenic","Intergenic","Intergenic","Intergenic","Intergenic","Intergenic"}}' $annot2feature_dir/${name}_${i}.genic.anno.tmp.bed $filter_dir/${name}_${i}.summit.bed > $annot2feature_dir/${name}_${i}.intergenic.anno.tmp.bed
cat $annot2feature_dir/${name}_${i}.genic.anno.tmp.bed $annot2feature_dir/${name}_${i}.intergenic.anno.tmp.bed | sort -k1,1 -k2,2n > $annot2feature_dir/${name}_${i}.anno.bed
bedtools merge -i $annot2feature_dir/${name}_${i}.anno.bed -s -c 4,10,13 -o first,distinct,distinct | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5,$6,$4,$7}' > $annot2feature_dir/${name}_${i}.annoMerge.bed
python parse_m6Aannotation.step06.py -i $annot2feature_dir/${name}_${i}.annoMerge.bed -o $annot2feature_dir/${name}_${i}.annoParse.bed

####anno2gene
bedtools intersect -a $filter_dir/${name}_${i}.summit.bed -b $hg38_gene -s -wa -wb > $annot2gene_dir/${name}_${i}.geneAnno.bed
awk 'BEGIN{OFS="\t"}{print $10}' $annot2gene_dir/${name}_${i}.geneAnno.bed | sort | uniq > $annot2gene_dir/${name}_${i}.geneList
}
done

rm $annot2feature_dir/*.tmp.bed
