# Initial configuration for AWS instance for nanopore sequence data analysis.

######
# install/update guppy, samtools, bedtools, python & pip
cd /Data1/software
sudo apt-get update
sudo apt-get install wget lsb-release
export PLATFORM=$(lsb_release -cs)
wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get update
sudo apt update
printf "Y" | sudo apt install samtools
sudo apt-get install bedtools
printf "Y" | sudo apt install ont-guppy --no-install-recommends

# obtain the rerio all-contexts model file
git clone https://github.com/nanoporetech/rerio
rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001
sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 /opt/ont/guppy/data/
sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg /opt/ont/guppy/data/
sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.jsn /opt/ont/guppy/data/
sudo cp -r /Data1/software/rerio/basecall_models/barcoding/* /opt/ont/guppy/data/barcoding/

#install/update megalodon
pip install megalodon
pip install ont_pyguppy_client_lib

#install latest versions of python & pip
sudo apt-get install python3.7
curl -O https://bootstrap.pypa.io/get-pip.py
python3 get-pip.py --user
######

######
# Download reference genome
#download T2T human reference sequence
cd /Data1/reference
# NOTE UPDATE BELOW FOR C ELEGANS REF GENOME
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
gunzip chm13.draft_v1.1.fasta.gz
######
