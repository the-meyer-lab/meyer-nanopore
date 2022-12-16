### Commands for visualization nanopore data using methylartist

# Create bed files for methylartist viz

# Create window labels for violin plot
bedtools makewindows -g ce11_size.txt -w 200000 > ce11_200kb_windows.bed
sed 's/chr/CHROMOSOME_/' ce11_200kb_windows.bed > ce11_200kb_windows_c.bed

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
awk -v f2=ce11_200kb_windows_c.bed ' { c = $1; getline < f2; print $0, c; } ' ce11_200kb_windows_c.bed  > ce11-200kb_windows_labels.bed

awk -v f2=sacCer3_10k_windows.txt ' { c = $1; getline < f2; print $0, c; } ' sacCer3_10kb_windows.txt  > sacCer3-10kb_windows_labels.bed

methylartist segmeth --bam /Data2/seq_data/Tube4_b2_2uM-Hia5_fiber-seq_11_21_22/basecalls/mod_mappings.sorted.m6Aonly.bam -i /Data2/reference/ce11_200kb_windows_c.bed -p 96

methylartist segmeth --bam /Data1/seq_data/210614_Raja/megalodon/barcode01/mod_mappings.bam -i /Data1/reference/ce11-10kb_windows_labels_X.bed -p 96
methylartist segmeth --bam /Data1/seq_data/210614_Raja/megalodon/barcode02_CpG/mod_mappings.sorted.bam -i /Data1/reference/sacCer3-10kb_windows_labels.bed -p 96
methylartist segmeth --bam /Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam -i /Data1/reference/ce11-chrm-regions_Xonly.bed -p 24

methylartist segmeth --bam /Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam -i /Data1/reference/ce11-100kb_windows_labels.bed -p 48
methylartist segmeth --bam /Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam -i /Data1/reference/ce11-100kb_windows_labels.bed -p 48

methylartist segmeth --bam /Data1/seq_data/20221003_TUBE27Vc_DPY27/basecalls/mod_mappings.sorted.bam -i /Data1/reference/ce11-100kb_windows_labels.bed -p 48

methylartist segplot -s ce11-10kb_windows_labels_X.mod_mappings.01.sorted.segmeth.tsv -v
methylartist segplot -s saccer3_10kb_windows_labels.mod_mappings.sorted.segmeth.tsv -v

methylartist locus --bams /Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam -i CHROMOSOME_X:2996577-2997578 -p 1,6,1,3,4
methylartist locus --bams /Data1/seq_data/210614_Raja/megalodon/barcode02_CpG/mod_mappings.sorted.bam -i chrI:100000-101000 -p 1,6,1,3,4

methylartist locus --bams /Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam,/Data1/seq_data/210614_Raja/megalodon/barcode01_m6A/mod_mappings.01.sorted.m6Aonly.bam -i CHROMOSOME_X:2992077-3002077 -l CHROMOSOME_X:2997071-2997083 --mods "a,m" --motifsize 1 -p 1,6,1,3,4


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

# For alfred QC
/Data1/sofware/alfred/src/alfred qc -r /Data1/reference/c_elegans.WS235.genomic.fa -j /Data1/viz/qc.json.gz -o /Data1/viz/qc.tsv.gz /Data1/seq_data/210614_Raja/megalodon/barcode01_CpG/mod_mappings.sorted.bam
