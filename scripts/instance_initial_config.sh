# Initial configuration for AWS instance for nanopore sequence data analysis.

###### Initiate file structure
# Navigate to home directory
cd ~

### Make folder tree (Uncomment on AWS server)
mkdir ./Data1/software
#mkdir ./Data1/seq_data
######

###### Install required software packages
# install/update guppy, samtools, bedtools, python & pip
cd ./Data1/software
# For remote server: ../Data1/software
sudo apt-get update
printf "Y" | sudo apt-get install wget lsb-release
export PLATFORM=$(lsb_release -cs)

wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get update
sudo apt update
printf "Y" | sudo apt upgrade
printf "Y" | sudo apt install samtools
sudo apt-get install bedtools
printf "Y" | sudo apt install ont-guppy --no-install-recommends

if [ ! -f /home/ubuntu/Data1/software/v2.03.tar.gz ]; then
  wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz
fi

tar zxvf v2.03.tar.gz
cd Winnowmap-2.03
make -j8
# cd ..
# mv Winnowmap-2.03 Winnowmap

#install/update megalodon
sudo pip install megalodon

#Install minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2 | tar -jxvf -

cd ../
######

###### Download required files
# obtain the rerio all-contexts model file
if [ ! -d /home/ubuntu/Data1/software/rerio ]; then
  git clone https://github.com/nanoporetech/rerio
  rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001
  sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 /opt/ont/guppy/data/
  sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg /opt/ont/guppy/data/
  sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.jsn /opt/ont/guppy/data/
  sudo cp -r /Data1/software/rerio/basecall_models/barcoding/* /opt/ont/guppy/data/barcoding/
else
  echo "rerio folder aready exists"
fi

# Download reference genome
# For remote server: ../Data1/reference
# Update below for your specific organism.
if [ ! -f /home/ubuntu/Data1/reference/c_elegans.WS235.genomic.fa ]; then
  mkdir ./Data1/reference
  wget https://downloads.wormbase.org/species/c_elegans/sequence/genomic/c_elegans.WS235.genomic.fa.gz -O ./reference/c_elegans.WS235.genomic.fa.gz
  gunzip ./reference/c_elegans.WS235.genomic.fa.gz
else
  echo "reference genome "reference/c_elegans.WS235.genomic.fa" already exists"
fi

#index using minimap
if [ ! -f /home/ubuntu/Data1/reference/ws235.mmi ]; then
  ./software/minimap2-2.20_x64-linux/minimap2 -d ./reference/ws235.mmi ./reference/c_elegans.WS235.genomic.fa
  else
  echo "indexed reference already exists"
fi

#download nanopore fast5 files
#note: replace s3 bucket with where your sequencing data is stored.
if [ ! -d /home/ubuntu/Data1/seq_data/210614_Raja ]; then
  aws s3 sync "s3://nanopore-1/nanopore first run/" /seq_data/
else
  echo "210614_Raja nanopore data already synced from S3."
fi
######




