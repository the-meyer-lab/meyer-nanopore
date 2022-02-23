#run guppy with demultiplexing to get barcode assignments for all reads (here it also does reference-free modification calling)
#change --barcode_kits parameter to match kit used
#1) initiate server in one screen

guppy_basecall_server -d ./Data1/software/rerio/basecall_models/ -c ./Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg --log_path ~/Data1/log --port 5555 --ipc_threads 2 --device cuda:all --num_callers 2 --gpu_runners_per_device 24 --chunks_per_runner 1024 --chunk_size 2000 --barcode_kits "SQK-LSK109" --trim_barcodes
#2) call supervisor in another screen once server is running
guppy_basecaller_supervisor --num_clients 16 --input_path ./Data1/seq_data/210614_Raja/20210614_2006_MN36964_FAQ21673_44f55e74/fast5 --save_path ./Data1/seq_data/210614_Raja/fastq --fast5_out --bam_out --bam_methylation_threshold 0 --config ~/Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg --port 5555 --barcode_kits "SQK-LSK109" --trim_barcodes

: << 'END'
#merge sequencing summary files
cd Data1/seq_dat/210614_Raja
awk 'FNR==1 && NR!=1{next;}{print}' sequencing_summary_*.txt >sequencing_summary_ALL.txt


#split fast5 files by barcode (for this sample barcode indices 13-24 were used)
#you could choose to split the raw input fast5s or the fast5s output by guppy (which contain base/mod call info)
#this submits the commands to run in the background to parallelize by barcode (by using & at the end of the line)
cd /Data1
for i in 13 14 15 16 17 18 19 20 21 22 23 24; do
	grep "barcode$i" /Data1/sequencing_summary_ALL.txt | cut -f 2 > /Data1/barcode$i.reads.txt
done
mkdir /Data1/fast5/split
for i in 13 14 15 16 17 18 19 20 21 22 23 24; do
	mkdir /Data1/fast5/split/barcode$i
	fast5_subset -i /Data1/fast5 -s /Data1/fast5/split/barcode$i -l /Data1/barcode$i.reads.txt -n 4000 -t 7 &
done


#now run megalodon on one or more individual split fast5 files
mkdir /Data1/megalodon
for i in 13 14; do
	megalodon /Data1/fast5/split/barcode$i --output-directory /Data1/megalodon/barcode$i --overwrite --guppy-params "-d /opt/ont/guppy/data" --guppy-server-path /usr/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg  --outputs mappings mod_mappings --reference /Data1/reference/chm13.draft_v1.1.fasta --device cuda:all --processes 92 --mod-min-prob 0
done
#sort megalodons output bam
for i in 13 14;do
	samtools sort -m 7000M -@ 4 -o /Data1/megalodon/barcode$i/mod_mappings.$i.sorted.bam /Data1/megalodon/barcode$i/mod_mappings.bam
done
END

#sync output files back over to s3 using a command like:
#cd /Data1
#rm -rf workspace
#aws s3 sync . s3://bucket/output_directory