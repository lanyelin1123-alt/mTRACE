wd=/path/to/project_work_folder
bam_dir=${wd}/result/bam/dedup
excldWT_bam_dir=${wd}/result/bam/exclude_wt
group=Result_name
ref=/path/to/genome.gtf

featurecount_folder=$wd/result/featurecount
mkdir -p $featurecount_folder

cd $bam_dir
featureCounts -s 1 -p -T 24 -t exon -g gene_id -a $ref -o $featurecount_folder/${group}_incldWT_exon.count *input.bam *ip.bam 2> $featurecount_folder/${group}_incldWT_exon.log 

cd $excldWT_bam_dir
featureCounts -s 1 -p -T 24 -t exon -g gene_id -a $ref -o $featurecount_folder/${group}_excldWT_exon.count *input.excldWT.bam *ip.excldWT.bam 2> $featurecount_folder/${group}_excldWT_exon.log
