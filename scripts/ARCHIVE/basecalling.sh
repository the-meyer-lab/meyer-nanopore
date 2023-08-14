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
guppy_basecaller_supervisor --num_clients 15 --input_path /Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5 --save_path /Data1/seq_data/210614_Raja/fastq_3_23 --fast5_out --bam_out --port 5555 --trim_barcodes --config ./Data1/software/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg --progress_stats_frequency 30 --barcode_kits EXP-NBD104

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