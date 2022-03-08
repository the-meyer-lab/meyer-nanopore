# Initial configuration for AWS instance for nanopore sequence data analysis.

###### Initiate file structure
# Navigate to home directory
cd ~

###### Install required software packages
cd ./Data1/software

# Add Oxford Nanopore's deb repository to your system (this is to install Oxford Nanopore Technologies-specific dependency packages):
# Get operating system name
sudo apt-get update
sudo apt-get update
printf "Y" | sudo apt-get install wget lsb-release
export PLATFORM=$(lsb_release -cs)
echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
sudo apt-get update

sudo apt update

# display ont-guppy requirements
apt-cache show ont-guppy
# Guppy requires Cuda 11.1.1
wget https://developer.download.nvidia.com/compute/cuda/11.1.1/local_installers/cuda_11.1.1_455.32.00_linux.run
sudo sh cuda_11.1.1_455.32.00_linux.run
#accept and install driver 455 & CUDA toolkit 11.1.1

sudo apt install nvidia-cuda-toolkit
export PATH="/usr/local/cuda-11.1/bin"

 -   LD_LIBRARY_PATH includes /usr/local/cuda-11.1/lib64, or, add /usr/local/cuda-11.1/lib64 to /etc/ld.so.conf and run ldconfig as root


#install samtools
printf "Y" | sudo apt install samtools

#install bedtools
sudo apt-get install bedtools

sudo apt-get install ont_guppy --no-install-recommends



#printf "Y" | sudo apt install ont-guppy --no-install-recommends
#if [ ! -f /home/ubuntu/Data1/software/ont_guppy_6.0.1-1~bionic_amd64.deb ]; then
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_5.0.16_linux64.tar.gz
#fi
#tar -xf ont-guppy_6.0.1_linuxaarch64_cuda10.tar.gz

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

### Install nvidia driver and CUDA
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.6.1/local_installers/cuda-repo-ubuntu1804-11-6-local_11.6.1-510.47.03-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu1804-11-6-local_11.6.1-510.47.03-1_amd64.deb
sudo apt-key add /var/cuda-repo-ubuntu1804-11-6-local/7fa2af80.pub
sudo apt-get update
sudo apt-get -y install cuda

export PATH=/usr/local/cuda-11.6/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-11.6/lib64\ ${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

# Install methylartist for analysis ad visualization
sudo pip install methylartist
sudo pip install methplotlib
sudo pip install methplotlib --upgrade

#install mysql needed to query chromosome sizes for methylartist violin plots
sudo apt install mysql-client-core-5.7
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




