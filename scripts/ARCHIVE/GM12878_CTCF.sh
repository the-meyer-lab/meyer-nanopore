### Purpose: Replicating CTCF targeted methylation with DiMeLo-seq
#__author__ = "Yuri Malina"
#__contact__ = "ymalina@berkeley.edu"
#__date__ = "4/10/2023"
#__status__ = "In development"
#__version__ = "0.0.1"

### Download appropriate reference file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Create minimap index
/Data1/software/minimap2-2.20_x64-linux/minimap2 -d /Data1/reference/hg38.mmi /Data1/reference/hg38.fa


### Run for basecalling GM12878 with hg38 reference 
megalodon /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/fast5 \
--output-directory /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/basecalls/m6A/ \
--overwrite \
--write-mods-text \
--mod-motif Y A 0 \
--outputs basecalls mod_basecalls mappings mod_mappings per_read_mods mappings mods \
--mod-output-formats wiggle \
--sort-mappings \
--guppy-params "-d /Data1/software/rerio/basecall_models --chunks_per_runner 250 --gpu_runners_per_device 24 --num_callers 64 --chunk_size 500" \
--guppy-server-path /Data1/software/ont-guppy_5.0.16/bin/guppy_basecall_server \
--guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
--reference /Data1/reference/hg38.mmi \
--device cuda:all \
--processes 92 \
--mod-min-prob 0 \
--guppy-timeout 2000

# DiMeLo QC Report (Remember to conda activate dimelo!)
/Data1/software/dimelo/bin/dimelo-qc-report \
-f /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/basecalls/m6A/mappings.sorted.bam \
-s "TubeAA1_GM12878_antiCTCF" \
-o /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/analysis/ -p 46

# DiMeLo CTCF
/Data1/software/dimelo/bin/dimelo-plot-enrichment-profile \
-f /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/basecalls/m6A/mod_mappings.sorted.bam \
-s "TubeAA1_GM12878_antiCTCF_m6Athresh232" \
-b /Data1/reference/ENCFF797SDL_100.bed \
-m A \
-o /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/analysis/ \
-A 232 \
-d 0.05