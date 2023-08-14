#!/bin/bash
# Requires conda activat megalodon

declare -a fast5_pathnames=("/home/ubuntu/Data1/seq_data/TubeY9_N2_fiberseq_timec_mpx3_17_23/") # "/my/second/path" )
basecall_5mC=0

for fast5_path in "${fast5_pathnames[@]}"
    do
        ### Basecall m6A
        megalodon \
        "${fast5_path}fast5/" \
        --output-directory "${fast5_path}basecalls/m6A_nompx/" \
        --overwrite \
        --write-mods-text \
        --mod-motif Y A 0 \
        --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods mappings mods \
        --mod-output-formats wiggle \
        --guppy-params "-d /home/ubuntu/Data1/software/rerio/basecall_models/" \
        --guppy-server-path "/home/ubuntu/Data1/software/ont-guppy_5.0.16/bin/guppy_basecall_server" \
        --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
        --reference /home/ubuntu/Data1/reference/ws235.mmi \
        --device cuda:all \
        --processes 50 \
        --mod-min-prob 0 \
        --guppy-timeout 2000
        #--guppy-concurrent-reads 4 \
        
        #3) sort megalodon mod_mappings.bam file
        samtools sort -m 3800M -@ 8 \
        -o "${fast5_path}basecalls/m6A/mod_mappings.sorted.bam" \
        "${fast5_path}basecalls/m6A/mod_mappings.bam"
        
        # For viewing in samtools, all context models should have their modification tags converted using ChEBI standards
        # See section 1.7 https://samtools.github.io/hts-specs/SAMtags.pdf
        samtools view -h -@ 8 "${fast5_path}basecalls/m6A/mod_mappings.sorted.bam" | \
        sed 's/A+Y/A+a/g' | \
        sed 's/C+Z/C+m/g' | \
        samtools view -@ 8 -bh -o "${fast5_path}basecalls/m6A/mod_mappings.sorted.m6Aonly.bam"
        
        # Index BAM file
        samtools index "${fast5_path}basecalls/m6A/mod_mappings.sorted.m6Aonly.bam"
        
        if [[ "$basecall_5mC" -eq 1 ]]; then
            ### Basecall 5mC
            megalodon \
            "${fast5_path}fast5/" \
            --output-directory "${fast5_path}basecalls/5mC/" \
            --overwrite \
            --write-mods-text \
            --mod-motif Y A 0 \
            --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods mappings mods \
            --mod-output-formats wiggle \
            --guppy-params "-d /Data1/software/rerio/basecall_models/" \
            --guppy-server-path /Data1/software/ont-guppy/bin/guppy_basecall_server \
            --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg \
            --reference /Data1/reference/ws235.mmi \
            --device cuda:all \
            --processes 48 \
            --mod-min-prob 0 \
            --guppy-concurrent-reads 4 \
            --guppy-timeout 2000
            
            #3) sort megalodon mod_mappings.bam file
            samtools sort -m 3800M -@ 8 \
            -o "${fast5_path}basecalls/5mC/mod_mappings.sorted.bam" \
            "${fast5_path}basecalls/5mC/mod_mappings.bam"
            
            # For viewing in samtools, all context models should have their modification tags converted using ChEBI standards
            # See section 1.7 https://samtools.github.io/hts-specs/SAMtags.pdf
            samtools view -h -@ 8 "${fast5_path}basecalls/5mC/mod_mappings.sorted.bam" | \
            sed 's/A+Y/A+a/g' | \
            sed 's/C+Z/C+m/g' | \
            samtools view -@ 8 -bh -o "${fast5_path}basecalls/5mC/mod_mappings.sorted.5mConly.bam"
            
            # Index BAM file
            samtools index "${fast5_path}basecalls/5mC/mod_mappings.sorted.5mConly.bam"
        fi
    done
    
'''
### Code run for basecalling GM12878 with hg38 reference 
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

# DiMeLo QC Report
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
-o /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/analysis_chm13/ \
-A 232 \
-d 0.05

# DiMeLo CTCF
/Data1/software/dimelo/bin/dimelo-plot-enrichment-profile \
-f /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/basecalls/m6A_chm13/mod_mappings.sorted.bam \
-s "TubeAA1_GM12878_antiCTCF_m6Athresh190" \
-b /Data1/reference/intersection.motifs.chip.formatted.chm13v1.1.sorted.top3k.bed \
-m A \
-o /Data1/seq_data/TubeAA1_GM12878_DiMeLo_aCTCF_4_9_2023/analysis_chm13/ \
-A 190 \
-d 0.05

'''