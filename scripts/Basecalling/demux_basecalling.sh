### Originally Authored by Nicolas Altemose 2022

#1) run megalodon on gpu instance & output basecalls + mod-mappings (& mappings just in case)
megalodon /fast5/path --output-directory /output/path --overwrite --guppy-params "-d /opt/ont/guppy/data" --guppy-server-path /usr/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg  --outputs basecalls mappings mod_mappings --reference /reference/path --device cuda:all --processes 10 --mod-min-prob 0

#2) run guppy_barcoder on the fastq outputs (use of GPU is optional, specified by --device)
mkdir /output/path/fastq
mv /output/path/basecalls.fastq /output/path/fastq
mkdir /output/path/demux
guppy_barcoder -i /output/path/fastq -s /output/path/demux --barcode_kits SQK-NBD110-24 --device cuda:all

#3) sort megalodon mod_mappings.bam file
cd /output/path
samtools sort -m 3800M -@ 8 -o mod_mappings.sorted.bam mod_mappings.bam

#4) pipe sorted mod_mappings and barcoding_summary info into perl script (e.g. "perl AddBarcodeToBam.pl /path/to/barcoding_summary.txt 04 05 06" would only bother to label barcodes 4,5,6 in the bam file, while putting anything else into a "0" bin)
samtools view -h mod_mappings.sorted.bam | perl AddBarcodeToBam.pl /output/path/demux/barcoding_summary.txt 04 05 06 | samtools view -b >mod_mappings_barcode.bam

#this can be trivially parellized by chromosome, for example:
for i in chr1 chr2 chr3;do
    samtools view -h mod_mappings.sorted.bam $i | perl AddBarcodeToBam.pl /output/path/demux/barcoding_summary.txt 04 05 06 | samtools view -b >mod_mappings_barcode.$i.bam &
done

#5) use samtools to split bam by barcode
samtools split -@ 10 -f '%*_%!.%.' mod_mappings_barcode.bam

#6) if using IGV for visualization, you need to adjust the mod tags a little bit

#this converts both the mA and mC tags to their proper form
samtools view -h -@ 8 mod_mappings.sorted.bam | sed 's/A+Y/A+a/g' | sed 's/C+Z/C+m/g' | samtools view -@ 8 -bh -o mod_mappings.sorted.m6AplusmC.bam

#7) if you want to only visualize mA or mC, or if you want to only visualize methyl calls above a certain threshold, this script can be used to apply different thresholds to mA and mC for IGV browsing (and if you set a threshold to 256 it will set all scores for that base to 0, allowing you to create a mA only or mC only file).
samtools view -h mod_mappings.sorted.m6AplusmC.bam | perl ThresholdBAMModCalls_separateAC_v0.pl <threshold mA, e.g. 155> <threshold mC, e.g. 256 to remove> | samtools view -b >mod_mappings.sorted.m6AplusmC.thresholded.bam

