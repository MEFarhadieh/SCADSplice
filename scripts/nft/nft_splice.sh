### Following codes were used FASTQ files to identify ASE & DTU & APA in nft
# barcodes filtering

ls ~/SCADSplice/scrnaseq/nft/trimed_data/*_2.fastq.gz | awk -F'/' '{print$6}' | awk -F'_' '{print$1}' >um.ID.txt
cat ~/um.ID.txt | while read id
do
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
 --stdin ~/SCADSplice/scrnaseq/nft/trimed_data/${id}_2.fastq.gz \
 --stdout ~/SCADSplice/scrnaseq/nft/trimed_data/${id}_extracted_2.fastq.gz \
 --read2-in ~/SCADSplice/scrnaseq/nft/trimed_data/${id}_3.fastq.gz \
 --read2-out=~/SCADSplice/scrnaseq/nft/trimed_data/${id}_extracted_3.fastq.gz \
 --whitelist=~/SCADSplice/scrnaseq/nft/results/metadata/nft_whitelist/${id}.txt
done

# alignment

ls ~/SCADSplice/scrnaseq/nft/trimed_data/*_extracted_3.fastq.gz | awk -F'/' '{print$6}' | awk -F'_' '{print$1}' >star.ID.txt
cat ~/star.ID.txt | while read id
do
STAR --runThreadN 12 \
--genomeDir ~/SCADSplice/references/SICILIAN_human_hg38_Refs/star_ref_file/hg38_ERCC_STAR_2.7.5.a \
--readFilesIn ~/SCADSplice/scrnaseq/nft/trimed_data/${id}_extracted_3.fastq.gz \
--readFilesCommand zcat \
--twopassMode Basic \
--alignIntronMax 1000000 \
--outFileNamePrefix ~/SCADSplice/scrnaseq/nft/results/${id}/2 \
--outSAMtype BAM Unsorted \
--outSAMattributes All \
--chimOutType WithinBAM SoftClip Junctions \
--chimJunctionOverhangMin 10 \
--chimSegmentReadGapMax 0 \
--chimOutJunctionFormat 1 \
--chimSegmentMin 12 \
--chimScoreJunctionNonGTAG -4 \
--chimNonchimScoreDropMin 10 \
--quantMode GeneCounts \
--sjdbGTFfile ~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf \
--outReadsUnmapped Fastx
done

# splice junction detection

git clone https://github.com/salzmanlab/SICILIAN.git
cd ~/SICILIAN

cat ~/star.ID.txt | while read id
do
python3 scripts/light_class_input.py --outpath ~/SCADSplice/scrnaseq/nft/results/${id}/ \
--gtf ~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf \
--annotator ~/SCADSplice/references/SICILIAN_human_hg38_Refs/annotator_file/hg38_refseq.pkl \
--bams ~/SCADSplice/scrnaseq/nft/results/${id}/2Aligned.out.bam \
--UMI_bar \
--stranded_library
done

cat ~/star.ID.txt | while read id
do
Rscript scripts/GLM_script_light.R ~/SCADSplice/scrnaseq/nft/results/${id}/ \
~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf  1  1  1 \
~/SCADSplice/references/SICILIAN_human_hg38_Refs/domain_file/ucscGenePfam.txt \
~/SCADSplice/references/SICILIAN_human_hg38_Refs/exon_pickle_file/hg38_refseq_exon_bounds.pkl \
~/SCADSplice/references/SICILIAN_human_hg38_Refs/splice_pickle_file/hg38_refseq_splices.pkl
done

cat ~/star.ID.txt | while read id
do
Rscript scripts/consolidate_GLM_output_files.R ~/SCADSplice/scrnaseq/nft/results/${id}/ \
${id} \
~/SCADSplice/references/SICILIAN_human_hg38_Refs/exon_pickle_file/hg38_refseq_exon_bounds.pkl \
~/SCADSplice/references/SICILIAN_human_hg38_Refs/splice_pickle_file/hg38_refseq_splices.pkl  1
done

cat ~/star.ID.txt | while read id
do
python3 scripts/Process_CI_10x.py -d ~/SCADSplice/scrnaseq/nft/results/${id}/ \
-o ${id} \
-g ~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf \
-e ~/SCADSplice/references/SICILIAN_human_hg38_Refs/exon_pickle_file/hg38_refseq_exon_bounds.pkl \
-s ~/SCADSplice/references/SICILIAN_human_hg38_Refs/splice_pickle_file/hg38_refseq_splices.pkl
done

cat ~/star.ID.txt | while read id
do
Rscript scripts/post_processing.R ~/SCADSplice/scrnaseq/nft/results/${id}/ ${id}  1
done
