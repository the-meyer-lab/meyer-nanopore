### Originally Authored by Nicolas Altemose 2022, modified by Yuri Malina

### m6A (note CpG calls in all context model are inaccurate!)
#1) run megalodon on gpu instance & output basecalls + mod-mappings (& mappings just in case)

# Note per-reads-mod required output otherwise weird error is received, need to figure out why: "Alphabet (ACGT) and model number of modified bases (-41) do not agree."
megalodon /Data2/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5/ --output-directory /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A --overwrite --guppy-params "-d /Data2/software/rerio/basecall_models/" --guppy-server-path /Data2/software/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --outputs basecalls mappings mod_mappings per_read_mods --reference /Data2/reference/ws235.mmi --device cuda:all --processes 92 --mod-min-prob 0

#2) run guppy_barcoder on the fastq outputs (use of GPU is optional, specified by --device)
mkdir /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/fastq
mv /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/basecalls.fastq /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/fastq
mkdir /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/demuxn
guppy_barcoder -i /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/fastq -s /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/demux --barcode_kits EXP-NBD104 --device cuda:all

#3) sort megalodon mod_mappings.bam file
cd /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/
samtools sort -m 3800M -@ 8 -o mod_mappings.sorted.bam mod_mappings.bam

# For viewing in samtools, all context models should have their modification tags converted using ChEBI standards
# See section 1.7 https://samtools.github.io/hts-specs/SAMtags.pdf
samtools view -h -@ 8 mod_mappings.01.sorted.bam | sed 's/A+Y/A+a/g' | sed 's/C+Z/C+m/g' | samtools view -@ 8 -bh -o mod_mappings.01.sorted.m6Aonly.bam

#4) pipe sorted mod_mappings and barcoding_summary info into perl script (e.g. "perl AddBarcodeToBam.pl /path/to/barcoding_summary.txt 04 05 06" would only bother to label barcodes 4,5,6 in the bam file, while putting anything else into a "0" bin)
samtools view -h mod_mappings.sorted.bam | perl AddBarcodeToBam.pl /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/demux/barcoding_summary.txt 01 02 | samtools view -b >mod_mappings_barcode.bam

#this can be trivially parellized by chromosome, for example:
# for i in chr1 chr2 chr3;do
#    samtools view -h mod_mappings.sorted.bam $i | perl AddBarcodeToBam.pl /output/path/demux/barcoding_summary.txt 04 05 06 | samtools view -b >mod_mappings_barcode.$i.bam &
# done

#5) use samtools to split bam by barcode
samtools split -@ 10 -f '%*_%!.%.' mod_mappings_barcode.bam

#6) if using IGV for visualization, you need to adjust the mod tags a little bit

#this only fixes the mA tag so can be used to effectively eliminate mC from visualization
# samtools view -h -@ 8 mod_mappings.sorted.bam | sed 's/A+Y/A+a/g' | samtools view -@ 8 -bh -o mod_mappings.sorted.m6Aonly.bam

# Index .bam file
samtools index mod_mappings_barcode.bam

#this converts both the mA and mC tags
# samtools view -h -@ 8 mod_mappings.sorted.bam | sed 's/A+Y/A+a/g' | sed 's/C+Z/C+m/g' | samtools view -@ 8 -bh -o mod_mappings.sorted.m6AplusmC.bam


### 5mC in CpG context
# Note per-reads-mod required output otherwise weird error is received, need to figure out why: "Alphabet (ACGT) and model number of modified bases (-41) do not agree."
megalodon /Data2/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5 --output-directory /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG --overwrite --guppy-params "-d /Data2/software/rerio/basecall_models/" --guppy-server-path /Data2/software/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg --outputs basecalls mappings mod_mappings --reference /Data2/reference/ws235.mmi --device cuda:all --processes 92 --mod-min-prob 0

#2) run guppy_barcoder on the fastq outputs (use of GPU is optional, specified by --device)
mkdir /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG/fastq
mv /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG/basecalls.fastq /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG/fastq
mkdir /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG/demux
guppy_barcoder -i /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG/fastq -s /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG/demux --barcode_kits EXP-NBD104 --device cuda:all

#3) sort megalodon mod_mappings.bam file
cd /Data2/seq_data/210614_Raja/megalodon_dmplx/CpG/
samtools sort -m 3800M -@ 8 -o mod_mappings.sorted.bam mod_mappings.bam

#4) pipe sorted mod_mappings and barcoding_summary info into perl script (e.g. "perl AddBarcodeToBam.pl /path/to/barcoding_summary.txt 04 05 06" would only bother to label barcodes 4,5,6 in the bam file, while putting anything else into a "0" bin)
samtools view -h mod_mappings.sorted.bam | perl AddBarcodeToBam.pl /Data2/seq_data/210614_Raja/megalodon_dmplx/m6A/demux/barcoding_summary.txt 01 02 | samtools view -b >mod_mappings_barcode.bam

#this can be trivially parellized by chromosome, for example:
# for i in chr1 chr2 chr3;do
#    samtools view -h mod_mappings.sorted.bam $i | perl AddBarcodeToBam.pl /output/path/demux/barcoding_summary.txt 04 05 06 | samtools view -b >mod_mappings_barcode.$i.bam &
# done

#5) use samtools to split bam by barcode
samtools split -@ 10 -f '%*_%!.%.' mod_mappings_barcode.bam

#m6A
megalodon /Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/fast5/ --output-directory /Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/basecalls/m6A/ --overwrite --guppy-params "-d /Data1/software/rerio/basecall_models/" --guppy-server-path /Data1/software/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods --reference /Data1/reference/ws235.mmi --device cuda:all --processes 92 --mod-min-prob 0

#5mC
megalodon /Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/fast5/ --output-directory /Data1/seq_data/TubeH1_021_SDC2-AIDpAux_Hia5_MSssI_12_19/basecalls/5mC/ --overwrite --guppy-params "-d /Data1/software/rerio/basecall_models/" --guppy-server-path /Data1/software/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg --outputs basecalls mod_basecalls mappings mod_mappings per_read_mods --reference /Data1/reference/ws235.mmi --device cuda:all --processes 80 --mod-min-prob 0
