### run guppy with demultiplexing to get barcode assignments for all reads (here it also does reference-free modification calling) ###

#1) initiate server in one screen
screen -S guppy_server
# Nic Altermose's original command:
# guppy_basecall_server -d ~/Data1/software/rerio/basecall_models/ -c ~/Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg --log_path ~/Data1/log --port 5555 --ipc_threads 16 --device cuda:all --num_callers 16 --gpu_runners_per_device 24 --chunks_per_runner 1024 --chunk_size 2000 --barcode_kits "SQK-LSK109" --trim_barcodes
# Note: the gpu configs spec'd were found to slow basecalling down. New command eleminates these, and runs with defaults.
# Note: since we are just using guppy for demultiplexing, chose a non-modified base model.
# change --barcode_kits parameter to match kit used. This can be found in your Oxford Nanopore account order invoice, or on packaging.
guppy_basecall_server -d /Data1/software/ont-guppy/data/ -c /Data1/software/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg --log_path /Data1/log --port 5555 --device cuda:all  --barcode_kit EXP-NBD104 --trim_barcodes
# ctrl + A; ctrl + D to detach from screen
# screen -ls lists all available screens
# screen -r [screen name] attaches to screen

#2) call supervisor in another screen once server is running
screen -S guppy_supervisor
# Nic Altemose's original command:
# guppy_basecaller_supervisor --num_clients 16 --input_path /Data2/fast5 --save_path /Data1 --fast5_out --bam_out --bam_methylation_threshold 0 --config res_dna_r941_min_modbases-all-context_v001.cfg --port 5555 --barcode_kits "SQK-NBD110-24" --trim_barcodes
# Note: since we are just using guppy for demultiplexing, chose a non-modified base model.
#Duration: this command takes ~3h on 8 GPUs for c.elegans
guppy_basecaller_supervisor --num_clients 15 --input_path /Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5 --save_path /Data1/seq_data/210614_Raja/fastq_3_1 --fast5_out --bam_out --port 5555 --trim_barcodes --config ./Data1/software/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg --progress_stats_frequency 30 --barcode_kits EXP-NBD104

# Our kits: Ligation Sequencing Kit SKU: SQK-LSK109 | SPOT ON FLOW CELL MK 1 R9 VERSION SKU: FLO-MIN106D | Native Barcoding Expansion 1-12 SKU: EXP-NBD104
# Use this command to determine correct config file: guppy_basecaller --print_workflows
#FLO-MIN106     SQK-LSK109                  dna_r9.4.1_450bps_hac          2021-05-17_dna_r9.4.1_minion_384_d37a2ab9

### MERGE SEQ SUMMARY FILTES ###
awk 'FNR==1 && NR!=1{next;}{print}' /Data1/seq_dat/210614_Raja/sequencing_summary_*.txt > /Data1/seq_dat/210614_Raja/sequencing_summary_ALL.txt


### SPLIT FAST5 FILES BY BARCODE ###
#you could choose to split the raw input fast5s or the fast5s output by guppy (which contain base/mod call info)
#this submits the commands to run in the background to parallelize by barcode (by using & at the end of the line)
# this command takes <3h for C Elegans
# Note: replace indicese with your barcode numbers of interest
for i in 01 02; do
	grep "barcode$i" sequencing_summary_ALL.txt | cut -f 2 > barcode$i.reads.txt
done

mkdir /Data1/fast5/split
for i in 01 02; do
  mkdir 20210614_2006_MN36964_FAQ21673_44f55e74/fast5/split/barcode$i
  fast5_subset -i 20210614_2006_MN36964_FAQ21673_44f55e74/fast5/ -s 20210614_2006_MN36964_FAQ21673_44f55e74/fast5/split/barcode$i/ -l fastq_3_1/barcode$i.reads.txt -n 4000 -t 7 &
done
#pkill -f fast5 to kill process

### RUN MEGALODON FOR EACH BARCODE FOR MOD BASE CALLS###
#Note update target input (ref) and output files!
#This command takes ~3h on 8GPUs for C Elegans.
mkdir /Data1/megalodon

# for c elegans, m6A:
for i in 01; do
	mkdir /Data2/megalodon_yeast/barcode$i
	megalodon /Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5/split/barcode$i --output-directory /Data1/seq_data/210614_Raja/megalodon/barcode$i --overwrite --guppy-params "-d /opt/ont/guppy/data" --guppy-server-path /usr/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg  --outputs mappings mod_mappings per_read_mods --write-mods-text --sort-mappings --reference /Data1/reference/ws235.mmi --device cuda:all --processes 92 --mod-min-prob 0
done
# Note per_read_mods outputs a db file with per-read per-position modified base scores.

# for C Elegans CpG 5mC
megalodon /Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5/split/barcode01 --output-directory /Data1/seq_data/210614_Raja/megalodon/barcode01_CpG --overwrite --guppy-params "-d /Data1/software/rerio/basecall_models/" --guppy-server-path /usr/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg  --outputs mappings mod_mappings per_read_mods --write-mods-text --sort-mappings --reference /Data1/reference/ws235.mmi --device cuda:all --processes 92 --mod-min-prob 0

# for yeast 6mA
megalodon /Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5/split/barcode02 --output-directory /Data1/seq_data/210614_Raja/megalodon/barcode02_m6A --overwrite --guppy-params "-d /Data1/software/rerio/basecall_models/" --guppy-server-path /usr/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --outputs mappings mod_mappings per_read_mods --write-mods-text --sort-mappings --reference /Data1/reference/sacCer3.fa --device cuda:all --processes 92 --mod-min-prob 0

# for yeast CpG 5mC
megalodon /Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5/split/barcode02 --output-directory /Data1/seq_data/210614_Raja/megalodon/barcode02_CpG --overwrite --guppy-params "-d /Data1/software/rerio/basecall_models/" --guppy-server-path /usr/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg --outputs mappings mod_mappings per_read_mods --write-mods-text --sort-mappings --reference /Data1/reference/sacCer3.fa --device cuda:all --processes 92 --mod-min-prob 0

# Create bed files for methylartist viz

# Create window labels for violin plot
bedtools makewindows -g ce11_size.txt -w 10000 > ce11_10kb_windows.bed
bedtools makewindows -g saccer3_size.txt -w 10000 > saccer3_10kb_windows.bed

# genome.txt file with the following format
# chrom	size
# chrV	20924180
# These can be found in per_read_mod_base_calls.db.chrm table
>sqlite3 per_read_modified_base_calls.db
sqlite> .headers on
sqlite> .mode tabs
sqlite> .output size.txt
sqlite> select chrm as chrom, chrm_len as size from chrm;
sqlite> .quit



# Add name column with name of chromsome
awk -v f2=ce11_10k_windows.txt ' { c = $1; getline < f2; print $0, c; } ' ce11_10kb_windows.txt  >ce11-10kb_windows_labels.bed
awk -v f2=sacCer3_10k_windows.txt ' { c = $1; getline < f2; print $0, c; } ' sacCer3_10kb_windows.txt  > sacCer3-10kb_windows_labels.bed

#Generate methylartist .megalodon database for use in segmeth and region plotting
methylartist db-megalodon -m /Data1/seq_data/210614_Raja/megalodon/barcode01/per_read_modified_base_calls.txt --db /Data2/ce11_ACmeth.db
#Note need to confirm that motifsize = 2
# Note this command currently does not work.

methylartist segmeth --bam /Data1/seq_data/210614_Raja/megalodon/barcode01/mod_mappings.bam -i /Data1/reference/ce11-10kb_windows_labels_X.bed -p 96
methylartist segmeth --bam /Data1/seq_data/210614_Raja/megalodon/barcode02_CpG/mod_mappings.sorted.bam -i /Data1/reference/sacCer3-10kb_windows_labels.bed -p 96


methylartist segplot -s ce11-10kb_windows_labels_X.mod_mappings.01.sorted.segmeth.tsv -v
methylartist segplot -s saccer3_10kb_windows_labels.mod_mappings.sorted.segmeth.tsv -v

methylartist locus --bam /Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam -i CHROMOSOME_X:2996577-2997578 -p 1,6,1,3,4
methylartist locus --bam /Data1/seq_data/210614_Raja/megalodon/barcode02_CpG/mod_mappings.sorted.bam -i chrI:100000-101000 -p 1,6,1,3,4


methylartist region -i chrX --bam /Data1/seq_data/210614_Raja/megalodon/barcode01/mod_mappings.01.sorted.bam -p 96 -n CG -r /Data1/reference/c_elegans.WS235.genomic.fa --skip_align_plot --panelratios 1,0,1,4 --height 4.5 --genepalette viridis --samplepalette viridis

# To determine  modified base alphabet
megalodon_extras modified_bases describe_alphabet --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --fast5s_dir /Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5/split/barcode01 --guppy-params "-d /opt/ont/guppy/data" --guppy-server-path /usr/bin/guppy_basecall_server
#output Loaded model calls canonical alphabet ACGT and modified bases Y=6mA (alt to A); Z=5mC (alt to C)

#sync output files back over to s3 using a command like:
#cd /Data1
#rm -rf workspace
#aws s3 sync . s3://bucket/output_directory
  aws s3 sync /seq_data/ "s3://nanopore-1/nanopore first run/"

# For viewing in samtools, all context models should have their modification tags converted using ChEBI standards
# See section 1.7 https://samtools.github.io/hts-specs/SAMtags.pdf
samtools view -h -@ 8 mod_mappings.01.sorted.bam | sed 's/A+Y/A+a/g' | sed 's/C+Z/C+m/g' | samtools view -@ 8 -bh -o mod_mappings.01.sorted.m6Aonly.bam

#This is required for viewing in IGV
#For details on setting up IGV with AWS S3 see: https://umccr.org/blog/igv-amazon-backend-setup/