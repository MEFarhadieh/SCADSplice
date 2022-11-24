### Following codes were used after triming FASTQ files to identify ASE in Astrocytes
# alignment

ls ~/SCADSplice/bulk/trimed_data/astro/*.gz | awk -F'/' '{print$6}' | awk -F'.' '{print$1}' >star.ID.txt

cat ~/star.ID.txt | while read id
do
STAR --runThreadN 4 --genomeDir ~/SCADSplice/references/SICILIAN_human_hg38_Refs/star_ref_file/hg38_ERCC_STAR_2.7.5.a \
 --readFilesIn ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix ~/SCADSplice/bulk/results/astro/${id} \
 --outBAMsortingThreadN 4
done

# sudo_alignment

cat ~/star.ID.txt |while read id
do
salmon quant -l A -i ~/SCADSplice/references/index/ -r ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz --seqBias --gcBias --validateMappings --rangeFactorizationBins 4 --threads 8 -o ~/SCADSplice/bulk/results/astro/salmon/${id}.quant
done

#rmats <conda install -y rmats rmats2sashimiplot>

cd ~/SCADSplice/bulk/results/astro/splicing/
ln -s ~/SCADSplice/bulk/results/astro/*.out.bam ./
ln -s ~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf ./
nano b1.txt
###{
./AAD1Aligned.sortedByCoord.out.bam,./AAD2Aligned.sortedByCoord.out.bam,./AAD3Aligned.sortedByCoord.out.bam,./AAD4Aligned.sortedByCoord.out.bam,./AAD5Aligned.sortedByCoord.out.bam,./AAD6Aligned.sortedByCoord.out.bam
###}
nano b2.txt
###{
./ACt1Aligned.sortedByCoord.out.bam,./ACt2Aligned.sortedByCoord.out.bam,./ACt3Aligned.sortedByCoord.out.bam,./ACt4Aligned.sortedByCoord.out.bam,./ACt5Aligned.sortedByCoord.out.bam,./ACt6Aligned.sortedByCoord.out.bam
###}

rmats.py --b1 ./b1.txt --b2 ./b2.txt --gtf ./grch38_known_genes.gtf -t single --readLength 50 --nthread 8 --od ./output --tmp ./tmp_output

#rmats2sashimiplot
#filter MXE events with FDR<0.05 and deltaPSI>0.1 to sig_MXE.txt

cat MXE.MATS.JC.txt | awk 'NR==1'>sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25>0.1)print$0}'>> sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25<(-0.1))print$0}'>> sig.MXE.txt

rmats2sashimiplot \
--b1 ../AAD1Aligned.sortedByCoord.out.bam,../AAD2Aligned.sortedByCoord.out.bam,../AAD3Aligned.sortedByCoord.out.bam,../AAD4Aligned.sortedByCoord.out.bam,../AAD5Aligned.sortedByCoord.out.bam,../AAD6Aligned.sortedByCoord.out.bam \
--b2 ../ACt1Aligned.sortedByCoord.out.bam,../ACt2Aligned.sortedByCoord.out.bam,../ACt3Aligned.sortedByCoord.out.bam,../ACt4Aligned.sortedByCoord.out.bam,../ACt5Aligned.sortedByCoord.out.bam,../ACt6Aligned.sortedByCoord.out.bam \
-t MXE -e sig.MXE.txt \
--l1 AAD --l2 ACt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/astro/Astro_MXE_plots/

#filter A3SS events with FDR<0.05 and deltaPSI>0.1 to sig_A3SS.txt

cat A3SS.MATS.JC.txt | awk 'NR==1'>sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A3SS.txt

rmats2sashimiplot \
--b1 ../AAD1Aligned.sortedByCoord.out.bam,../AAD2Aligned.sortedByCoord.out.bam,../AAD3Aligned.sortedByCoord.out.bam,../AAD4Aligned.sortedByCoord.out.bam,../AAD5Aligned.sortedByCoord.out.bam,../AAD6Aligned.sortedByCoord.out.bam \
--b2 ../ACt1Aligned.sortedByCoord.out.bam,../ACt2Aligned.sortedByCoord.out.bam,../ACt3Aligned.sortedByCoord.out.bam,../ACt4Aligned.sortedByCoord.out.bam,../ACt5Aligned.sortedByCoord.out.bam,../ACt6Aligned.sortedByCoord.out.bam \
-t A3SS -e sig.A3SS.txt \
--l1 AAD --l2 ACt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/astro/Astro_A3SS_plots/

#filter A5SS events with FDR<0.05 and deltaPSI>0.1 to sig_A5SS.txt

cat A5SS.MATS.JC.txt | awk 'NR==1'>sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A5SS.txt

rmats2sashimiplot \
--b1 ../AAD1Aligned.sortedByCoord.out.bam,../AAD2Aligned.sortedByCoord.out.bam,../AAD3Aligned.sortedByCoord.out.bam,../AAD4Aligned.sortedByCoord.out.bam,../AAD5Aligned.sortedByCoord.out.bam,../AAD6Aligned.sortedByCoord.out.bam \
--b2 ../ACt1Aligned.sortedByCoord.out.bam,../ACt2Aligned.sortedByCoord.out.bam,../ACt3Aligned.sortedByCoord.out.bam,../ACt4Aligned.sortedByCoord.out.bam,../ACt5Aligned.sortedByCoord.out.bam,../ACt6Aligned.sortedByCoord.out.bam \
-t A5SS -e sig.A5SS.txt \
--l1 AAD --l2 ACt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/astro/Astro_A5SS_plots/

#filter SE events with FDR<0.05 and deltaPSI>0.1 to sig_SE.txt

cat SE.MATS.JC.txt | awk 'NR==1'>sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.SE.txt

rmats2sashimiplot \
--b1 ../AAD1Aligned.sortedByCoord.out.bam,../AAD2Aligned.sortedByCoord.out.bam,../AAD3Aligned.sortedByCoord.out.bam,../AAD4Aligned.sortedByCoord.out.bam,../AAD5Aligned.sortedByCoord.out.bam,../AAD6Aligned.sortedByCoord.out.bam \
--b2 ../ACt1Aligned.sortedByCoord.out.bam,../ACt2Aligned.sortedByCoord.out.bam,../ACt3Aligned.sortedByCoord.out.bam,../ACt4Aligned.sortedByCoord.out.bam,../ACt5Aligned.sortedByCoord.out.bam,../ACt6Aligned.sortedByCoord.out.bam \
-t SE -e sig.SE.txt \
--l1 AAD --l2 ACt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/astro/Astro_SE_plots/

#filter RI events with FDR<0.05 and deltaPSI>0.1 to sig_RI.txt

cat RI.MATS.JC.txt | awk 'NR==1'>sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.RI.txt

rmats2sashimiplot \
--b1 ../AAD1Aligned.sortedByCoord.out.bam,../AAD2Aligned.sortedByCoord.out.bam,../AAD3Aligned.sortedByCoord.out.bam,../AAD4Aligned.sortedByCoord.out.bam,../AAD5Aligned.sortedByCoord.out.bam,../AAD6Aligned.sortedByCoord.out.bam \
--b2 ../ACt1Aligned.sortedByCoord.out.bam,../ACt2Aligned.sortedByCoord.out.bam,../ACt3Aligned.sortedByCoord.out.bam,../ACt4Aligned.sortedByCoord.out.bam,../ACt5Aligned.sortedByCoord.out.bam,../ACt6Aligned.sortedByCoord.out.bam \
-t RI -e sig.RI.txt \
--l1 AAD --l2 ACt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/astro/Astro_RI_plots/

### Following codes were used after triming FASTQ files to identify ASE in Endothelial
# alignment

ls ~/SCADSplice/bulk/trimed_data/endo/*.gz | awk -F'/' '{print$6}' | awk -F'.' '{print$1}' >star.ID.txt

cat ~/star.ID.txt | while read id
do
STAR --runThreadN 4 --genomeDir ~/SCADSplice/references/SICILIAN_human_hg38_Refs/star_ref_file/hg38_ERCC_STAR_2.7.5.a \
 --readFilesIn ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix ~/SCADSplice/bulk/results/endo/${id} \
 --outBAMsortingThreadN 4
done

# sudo_alignment

cat ~/star.ID.txt |while read id
do
salmon quant -l A -i ~/SCADSplice/references/index/ -r ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz --seqBias --gcBias --validateMappings --rangeFactorizationBins 4 --threads 8 -o ~/SCADSplice/bulk/results/endo/salmon/${id}.quant
done

#rmats <conda install -y rmats rmats2sashimiplot>

cd ~/SCADSplice/bulk/results/endo/splicing/
ln -s ~/SCADSplice/bulk/results/endo/*.out.bam ./
ln -s ~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf ./
nano b1.txt
###{
./EAD1Aligned.sortedByCoord.out.bam,./EAD2Aligned.sortedByCoord.out.bam,./EAD3Aligned.sortedByCoord.out.bam,./EAD4Aligned.sortedByCoord.out.bam,./EAD5Aligned.sortedByCoord.out.bam,./EAD6Aligned.sortedByCoord.out.bam
###}
nano b2.txt
###{
./ECt1Aligned.sortedByCoord.out.bam,./ECt2Aligned.sortedByCoord.out.bam,./ECt3Aligned.sortedByCoord.out.bam,./ECt4Aligned.sortedByCoord.out.bam,./ECt5Aligned.sortedByCoord.out.bam,./ECt6Aligned.sortedByCoord.out.bam
###}

rmats.py --b1 ./b1.txt --b2 ./b2.txt --gtf ./grch38_known_genes.gtf -t single --readLength 50 --nthread 8 --od ./output --tmp ./tmp_output

#rmats2sashimiplot
#filter MXE events with FDR<0.05 and deltaPSI>0.1 to sig_MXE.txt

cat MXE.MATS.JC.txt | awk 'NR==1'>sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25>0.1)print$0}'>> sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25<(-0.1))print$0}'>> sig.MXE.txt

rmats2sashimiplot \
--b1 ../EAD1Aligned.sortedByCoord.out.bam,../EAD2Aligned.sortedByCoord.out.bam,../EAD3Aligned.sortedByCoord.out.bam,../EAD4Aligned.sortedByCoord.out.bam,../EAD5Aligned.sortedByCoord.out.bam,../EAD6Aligned.sortedByCoord.out.bam \
--b2 ../ECt1Aligned.sortedByCoord.out.bam,../ECt2Aligned.sortedByCoord.out.bam,../ECt3Aligned.sortedByCoord.out.bam,../ECt4Aligned.sortedByCoord.out.bam,../ECt5Aligned.sortedByCoord.out.bam,../ECt6Aligned.sortedByCoord.out.bam \
-t MXE -e sig.MXE.txt \
--l1 EAD --l2 ECt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/endo/endo_MXE_plots/

#filter A3SS events with FDR<0.05 and deltaPSI>0.1 to sig_A3SS.txt

cat A3SS.MATS.JC.txt | awk 'NR==1'>sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A3SS.txt

rmats2sashimiplot \
--b1 ../EAD1Aligned.sortedByCoord.out.bam,../EAD2Aligned.sortedByCoord.out.bam,../EAD3Aligned.sortedByCoord.out.bam,../EAD4Aligned.sortedByCoord.out.bam,../EAD5Aligned.sortedByCoord.out.bam,../EAD6Aligned.sortedByCoord.out.bam \
--b2 ../ECt1Aligned.sortedByCoord.out.bam,../ECt2Aligned.sortedByCoord.out.bam,../ECt3Aligned.sortedByCoord.out.bam,../ECt4Aligned.sortedByCoord.out.bam,../ECt5Aligned.sortedByCoord.out.bam,../ECt6Aligned.sortedByCoord.out.bam \
-t A3SS -e sig.A3SS.txt \
--l1 EAD --l2 ECt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/endo/endo_A3SS_plots/

#filter A5SS events with FDR<0.05 and deltaPSI>0.1 to sig_A5SS.txt

cat A5SS.MATS.JC.txt | awk 'NR==1'>sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A5SS.txt

rmats2sashimiplot \
--b1 ../EAD1Aligned.sortedByCoord.out.bam,../EAD2Aligned.sortedByCoord.out.bam,../EAD3Aligned.sortedByCoord.out.bam,../EAD4Aligned.sortedByCoord.out.bam,../EAD5Aligned.sortedByCoord.out.bam,../EAD6Aligned.sortedByCoord.out.bam \
--b2 ../ECt1Aligned.sortedByCoord.out.bam,../ECt2Aligned.sortedByCoord.out.bam,../ECt3Aligned.sortedByCoord.out.bam,../ECt4Aligned.sortedByCoord.out.bam,../ECt5Aligned.sortedByCoord.out.bam,../ECt6Aligned.sortedByCoord.out.bam \
-t A5SS -e sig.A5SS.txt \
--l1 EAD --l2 ECt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/endo/endo_A5SS_plots/

#filter SE events with FDR<0.05 and deltaPSI>0.1 to sig_SE.txt

cat SE.MATS.JC.txt | awk 'NR==1'>sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.SE.txt

rmats2sashimiplot \
--b1 ../EAD1Aligned.sortedByCoord.out.bam,../EAD2Aligned.sortedByCoord.out.bam,../EAD3Aligned.sortedByCoord.out.bam,../EAD4Aligned.sortedByCoord.out.bam,../EAD5Aligned.sortedByCoord.out.bam,../EAD6Aligned.sortedByCoord.out.bam \
--b2 ../ECt1Aligned.sortedByCoord.out.bam,../ECt2Aligned.sortedByCoord.out.bam,../ECt3Aligned.sortedByCoord.out.bam,../ECt4Aligned.sortedByCoord.out.bam,../ECt5Aligned.sortedByCoord.out.bam,../ECt6Aligned.sortedByCoord.out.bam \
-t SE -e sig.SE.txt \
--l1 EAD --l2 ECt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/endo/endo_SE_plots/

#filter RI events with FDR<0.05 and deltaPSI>0.1 to sig_RI.txt

cat RI.MATS.JC.txt | awk 'NR==1'>sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.RI.txt

rmats2sashimiplot \
--b1 ../EAD1Aligned.sortedByCoord.out.bam,../EAD2Aligned.sortedByCoord.out.bam,../EAD3Aligned.sortedByCoord.out.bam,../EAD4Aligned.sortedByCoord.out.bam,../EAD5Aligned.sortedByCoord.out.bam,../EAD6Aligned.sortedByCoord.out.bam \
--b2 ../ECt1Aligned.sortedByCoord.out.bam,../ECt2Aligned.sortedByCoord.out.bam,../ECt3Aligned.sortedByCoord.out.bam,../ECt4Aligned.sortedByCoord.out.bam,../ECt5Aligned.sortedByCoord.out.bam,../ECt6Aligned.sortedByCoord.out.bam \
-t RI -e sig.RI.txt \
--l1 EAD --l2 ECt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/endo/endo_RI_plots/

### Following codes were used after triming FASTQ files to identify ASE in Microglia
# alignment

ls ~/SCADSplice/bulk/trimed_data/micro/*.gz | awk -F'/' '{print$6}' | awk -F'.' '{print$1}' >star.ID.txt

cat ~/star.ID.txt | while read id
do
STAR --runThreadN 4 --genomeDir ~/SCADSplice/references/SICILIAN_human_hg38_Refs/star_ref_file/hg38_ERCC_STAR_2.7.5.a \
 --readFilesIn ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix ~/SCADSplice/bulk/results/micro/${id} \
 --outBAMsortingThreadN 4
done

# sudo_alignment

cat ~/star.ID.txt |while read id
do
salmon quant -l A -i ~/SCADSplice/references/index/ -r ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz --seqBias --gcBias --validateMappings --rangeFactorizationBins 4 --threads 8 -o ~/SCADSplice/bulk/results/micro/salmon/${id}.quant
done

#rmats <conda install -y rmats rmats2sashimiplot>

cd ~/SCADSplice/bulk/results/micro/splicing/
ln -s ~/SCADSplice/bulk/results/micro/*.out.bam ./
ln -s ~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf ./
nano b1.txt
###{
./MAD1Aligned.sortedByCoord.out.bam,./MAD2Aligned.sortedByCoord.out.bam,./MAD3Aligned.sortedByCoord.out.bam,./MAD4Aligned.sortedByCoord.out.bam,./MAD5Aligned.sortedByCoord.out.bam,./MAD6Aligned.sortedByCoord.out.bam
###}
nano b2.txt
###{
./MCt1Aligned.sortedByCoord.out.bam,./MCt2Aligned.sortedByCoord.out.bam,./MCt3Aligned.sortedByCoord.out.bam,./MCt4Aligned.sortedByCoord.out.bam,./MCt5Aligned.sortedByCoord.out.bam,./MCt6Aligned.sortedByCoord.out.bam
###}

rmats.py --b1 ./b1.txt --b2 ./b2.txt --gtf ./grch38_known_genes.gtf -t single --readLength 50 --nthread 8 --od ./output --tmp ./tmp_output

#rmats2sashimiplot
#filter MXE events with FDR<0.05 and deltaPSI>0.1 to sig_MXE.txt

cat MXE.MATS.JC.txt | awk 'NR==1'>sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25>0.1)print$0}'>> sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25<(-0.1))print$0}'>> sig.MXE.txt

rmats2sashimiplot \
--b1 ../MAD1Aligned.sortedByCoord.out.bam,../MAD2Aligned.sortedByCoord.out.bam,../MAD3Aligned.sortedByCoord.out.bam,../MAD4Aligned.sortedByCoord.out.bam,../MAD5Aligned.sortedByCoord.out.bam,../MAD6Aligned.sortedByCoord.out.bam \
--b2 ../MCt1Aligned.sortedByCoord.out.bam,../MCt2Aligned.sortedByCoord.out.bam,../MCt3Aligned.sortedByCoord.out.bam,../MCt4Aligned.sortedByCoord.out.bam,../MCt5Aligned.sortedByCoord.out.bam,../MCt6Aligned.sortedByCoord.out.bam \
-t MXE -e sig.MXE.txt \
--l1 MAD --l2 MCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/micro/micro_MXE_plots/

#filter A3SS events with FDR<0.05 and deltaPSI>0.1 to sig_A3SS.txt

cat A3SS.MATS.JC.txt | awk 'NR==1'>sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A3SS.txt

rmats2sashimiplot \
--b1 ../MAD1Aligned.sortedByCoord.out.bam,../MAD2Aligned.sortedByCoord.out.bam,../MAD3Aligned.sortedByCoord.out.bam,../MAD4Aligned.sortedByCoord.out.bam,../MAD5Aligned.sortedByCoord.out.bam,../MAD6Aligned.sortedByCoord.out.bam \
--b2 ../MCt1Aligned.sortedByCoord.out.bam,../MCt2Aligned.sortedByCoord.out.bam,../MCt3Aligned.sortedByCoord.out.bam,../MCt4Aligned.sortedByCoord.out.bam,../MCt5Aligned.sortedByCoord.out.bam,../MCt6Aligned.sortedByCoord.out.bam \
-t A3SS -e sig.A3SS.txt \
--l1 MAD --l2 MCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/micro/micro_A3SS_plots/

#filter A5SS events with FDR<0.05 and deltaPSI>0.1 to sig_A5SS.txt

cat A5SS.MATS.JC.txt | awk 'NR==1'>sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A5SS.txt

rmats2sashimiplot \
--b1 ../MAD1Aligned.sortedByCoord.out.bam,../MAD2Aligned.sortedByCoord.out.bam,../MAD3Aligned.sortedByCoord.out.bam,../MAD4Aligned.sortedByCoord.out.bam,../MAD5Aligned.sortedByCoord.out.bam,../MAD6Aligned.sortedByCoord.out.bam \
--b2 ../MCt1Aligned.sortedByCoord.out.bam,../MCt2Aligned.sortedByCoord.out.bam,../MCt3Aligned.sortedByCoord.out.bam,../MCt4Aligned.sortedByCoord.out.bam,../MCt5Aligned.sortedByCoord.out.bam,../MCt6Aligned.sortedByCoord.out.bam \
-t A5SS -e sig.A5SS.txt \
--l1 MAD --l2 MCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/micro/micro_A5SS_plots/

#filter SE events with FDR<0.05 and deltaPSI>0.1 to sig_SE.txt

cat SE.MATS.JC.txt | awk 'NR==1'>sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.SE.txt

rmats2sashimiplot \
--b1 ../MAD1Aligned.sortedByCoord.out.bam,../MAD2Aligned.sortedByCoord.out.bam,../MAD3Aligned.sortedByCoord.out.bam,../MAD4Aligned.sortedByCoord.out.bam,../MAD5Aligned.sortedByCoord.out.bam,../MAD6Aligned.sortedByCoord.out.bam \
--b2 ../MCt1Aligned.sortedByCoord.out.bam,../MCt2Aligned.sortedByCoord.out.bam,../MCt3Aligned.sortedByCoord.out.bam,../MCt4Aligned.sortedByCoord.out.bam,../MCt5Aligned.sortedByCoord.out.bam,../MCt6Aligned.sortedByCoord.out.bam \
-t SE -e sig.SE.txt \
--l1 MAD --l2 MCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/micro/micro_SE_plots/

#filter RI events with FDR<0.05 and deltaPSI>0.1 to sig_RI.txt

cat RI.MATS.JC.txt | awk 'NR==1'>sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.RI.txt

rmats2sashimiplot \
--b1 ../MAD1Aligned.sortedByCoord.out.bam,../MAD2Aligned.sortedByCoord.out.bam,../MAD3Aligned.sortedByCoord.out.bam,../MAD4Aligned.sortedByCoord.out.bam,../MAD5Aligned.sortedByCoord.out.bam,../MAD6Aligned.sortedByCoord.out.bam \
--b2 ../MCt1Aligned.sortedByCoord.out.bam,../MCt2Aligned.sortedByCoord.out.bam,../MCt3Aligned.sortedByCoord.out.bam,../MCt4Aligned.sortedByCoord.out.bam,../MCt5Aligned.sortedByCoord.out.bam,../MCt6Aligned.sortedByCoord.out.bam \
-t RI -e sig.RI.txt \
--l1 MAD --l2 MCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/micro/micro_RI_plots/

### Following codes were used after triming FASTQ files to identify ASE in Neuron
# alignment

ls ~/SCADSplice/bulk/trimed_data/neuro/*.gz | awk -F'/' '{print$6}' | awk -F'.' '{print$1}' >star.ID.txt

cat ~/star.ID.txt | while read id
do
STAR --runThreadN 4 --genomeDir ~/SCADSplice/references/SICILIAN_human_hg38_Refs/star_ref_file/hg38_ERCC_STAR_2.7.5.a \
 --readFilesIn ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate \
 --outFileNamePrefix ~/SCADSplice/bulk/results/neuro/${id} \
 --outBAMsortingThreadN 4
done

# sudo_alignment

cat ~/star.ID.txt |while read id
do
salmon quant -l A -i ~/SCADSplice/references/index/ -r ~/SCADSplice/bulk/data/${id}_trimmed.fq.gz --seqBias --gcBias --validateMappings --rangeFactorizationBins 4 --threads 8 -o ~/SCADSplice/bulk/results/neuro/salmon/${id}.quant
done

#rmats <conda install -y rmats rmats2sashimiplot>

cd ~/SCADSplice/bulk/results/neuro/splicing/
ln -s ~/SCADSplice/bulk/results/neuro/*.out.bam ./
ln -s ~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf ./
nano b1.txt
###{
./NAD1Aligned.sortedByCoord.out.bam,./NAD2Aligned.sortedByCoord.out.bam,./NAD3Aligned.sortedByCoord.out.bam,./NAD4Aligned.sortedByCoord.out.bam,./NAD5Aligned.sortedByCoord.out.bam,./NAD6Aligned.sortedByCoord.out.bam
###}
nano b2.txt
###{
./NCt1Aligned.sortedByCoord.out.bam,./NCt2Aligned.sortedByCoord.out.bam,./NCt3Aligned.sortedByCoord.out.bam,./NCt4Aligned.sortedByCoord.out.bam,./NCt5Aligned.sortedByCoord.out.bam,./NCt6Aligned.sortedByCoord.out.bam
###}

rmats.py --b1 ./b1.txt --b2 ./b2.txt --gtf ./grch38_known_genes.gtf -t single --readLength 50 --nthread 8 --od ./output --tmp ./tmp_output

#rmats2sashimiplot
#filter MXE events with FDR<0.05 and deltaPSI>0.1 to sig_MXE.txt

cat MXE.MATS.JC.txt | awk 'NR==1'>sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25>0.1)print$0}'>> sig.MXE.txt
cat MXE.MATS.JC.txt | awk -F '\t' '{if($22<0.05 && $25<(-0.1))print$0}'>> sig.MXE.txt

rmats2sashimiplot \
--b1 ../NAD1Aligned.sortedByCoord.out.bam,../NAD2Aligned.sortedByCoord.out.bam,../NAD3Aligned.sortedByCoord.out.bam,../NAD4Aligned.sortedByCoord.out.bam,../NAD5Aligned.sortedByCoord.out.bam,../NAD6Aligned.sortedByCoord.out.bam \
--b2 ../NCt1Aligned.sortedByCoord.out.bam,../NCt2Aligned.sortedByCoord.out.bam,../NCt3Aligned.sortedByCoord.out.bam,../NCt4Aligned.sortedByCoord.out.bam,../NCt5Aligned.sortedByCoord.out.bam,../NCt6Aligned.sortedByCoord.out.bam \
-t MXE -e sig.MXE.txt \
--l1 NAD --l2 NCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/neuro/neuro_MXE_plots/

#filter A3SS events with FDR<0.05 and deltaPSI>0.1 to sig_A3SS.txt

cat A3SS.MATS.JC.txt | awk 'NR==1'>sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A3SS.txt
cat A3SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A3SS.txt

rmats2sashimiplot \
--b1 ../NAD1Aligned.sortedByCoord.out.bam,../NAD2Aligned.sortedByCoord.out.bam,../NAD3Aligned.sortedByCoord.out.bam,../NAD4Aligned.sortedByCoord.out.bam,../NAD5Aligned.sortedByCoord.out.bam,../NAD6Aligned.sortedByCoord.out.bam \
--b2 ../NCt1Aligned.sortedByCoord.out.bam,../NCt2Aligned.sortedByCoord.out.bam,../NCt3Aligned.sortedByCoord.out.bam,../NCt4Aligned.sortedByCoord.out.bam,../NCt5Aligned.sortedByCoord.out.bam,../NCt6Aligned.sortedByCoord.out.bam \
-t A3SS -e sig.A3SS.txt \
--l1 NAD --l2 NCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/neuro/neuro_A3SS_plots/

#filter A5SS events with FDR<0.05 and deltaPSI>0.1 to sig_A5SS.txt

cat A5SS.MATS.JC.txt | awk 'NR==1'>sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.A5SS.txt
cat A5SS.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.A5SS.txt

rmats2sashimiplot \
--b1 ../NAD1Aligned.sortedByCoord.out.bam,../NAD2Aligned.sortedByCoord.out.bam,../NAD3Aligned.sortedByCoord.out.bam,../NAD4Aligned.sortedByCoord.out.bam,../NAD5Aligned.sortedByCoord.out.bam,../NAD6Aligned.sortedByCoord.out.bam \
--b2 ../NCt1Aligned.sortedByCoord.out.bam,../NCt2Aligned.sortedByCoord.out.bam,../NCt3Aligned.sortedByCoord.out.bam,../NCt4Aligned.sortedByCoord.out.bam,../NCt5Aligned.sortedByCoord.out.bam,../NCt6Aligned.sortedByCoord.out.bam \
-t A5SS -e sig.A5SS.txt \
--l1 NAD --l2 NCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/neuro/neuro_A5SS_plots/

#filter SE events with FDR<0.05 and deltaPSI>0.1 to sig_SE.txt

cat SE.MATS.JC.txt | awk 'NR==1'>sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.SE.txt
cat SE.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.SE.txt

rmats2sashimiplot \
--b1 ../NAD1Aligned.sortedByCoord.out.bam,../NAD2Aligned.sortedByCoord.out.bam,../NAD3Aligned.sortedByCoord.out.bam,../NAD4Aligned.sortedByCoord.out.bam,../NAD5Aligned.sortedByCoord.out.bam,../NAD6Aligned.sortedByCoord.out.bam \
--b2 ../NCt1Aligned.sortedByCoord.out.bam,../NCt2Aligned.sortedByCoord.out.bam,../NCt3Aligned.sortedByCoord.out.bam,../NCt4Aligned.sortedByCoord.out.bam,../NCt5Aligned.sortedByCoord.out.bam,../NCt6Aligned.sortedByCoord.out.bam \
-t SE -e sig.SE.txt \
--l1 NAD --l2 NCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/neuro/neuro_SE_plots/

#filter RI events with FDR<0.05 and deltaPSI>0.1 to sig_RI.txt

cat RI.MATS.JC.txt | awk 'NR==1'>sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23>0.1)print$0}'>> sig.RI.txt
cat RI.MATS.JC.txt | awk -F '\t' '{if($20<0.05 && $23<(-0.1))print$0}'>> sig.RI.txt

rmats2sashimiplot \
--b1 ../NAD1Aligned.sortedByCoord.out.bam,../NAD2Aligned.sortedByCoord.out.bam,../NAD3Aligned.sortedByCoord.out.bam,../NAD4Aligned.sortedByCoord.out.bam,../NAD5Aligned.sortedByCoord.out.bam,../NAD6Aligned.sortedByCoord.out.bam \
--b2 ../NCt1Aligned.sortedByCoord.out.bam,../NCt2Aligned.sortedByCoord.out.bam,../NCt3Aligned.sortedByCoord.out.bam,../NCt4Aligned.sortedByCoord.out.bam,../NCt5Aligned.sortedByCoord.out.bam,../NCt6Aligned.sortedByCoord.out.bam \
-t RI -e sig.RI.txt \
--l1 NAD --l2 NCt \
--group-info group.gf \
-o ~/SCADSplice/bulk/results/neuro/neuro_RI_plots/
