# Initial configuration for AWS instance for nanopore sequence data analysis.

###### Initiate file structure
# Create folder for software and navigate to it
if [ ! -f /Data1/software ]; then
  mkdir /Data1/software
fi
cd /Data1/software

###### Install required software packages

### Install Ont-Guppy
# DO NOT USE THESE COMMAND FOR TIME BEING. THEY BREAK GRAPHICS CARD DRIVERS. USE NEXT SECTION TO INSTALL FROM TAR
# Add Oxford Nanopore's deb repository to your system (this is to install Oxford Nanopore Technologies-specific dependency packages):
# Taken from installation manual here: https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revac_14dec2018/linux-guppy
#sudo apt-get update
#printf "Y" | sudo apt-get install wget lsb-release
#export PLATFORM=$(lsb_release -cs)
#echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list

#sudo apt-get update
#sudo apt update

# If error saying signatures couldn't be verified because public key is not available, copy public key and replace in below command:
# sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys AD9D09AF2EBBA5A8
# sudo apt-get update
# sudo apt update
#sudo apt install ont-guppy

### Install Ont-Guppt from Tar

if [ ! -f /home/ubuntu/Data1/software/ont-guppy_5.0.16_linux64.tar.gz ]; then
  wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_5.0.16_linux64.tar.gz
  tar -xf ont-guppy_5.0.16_linux64.tar.gz
  pip install ont-pyguppy-client-lib==5.0.16
  #add this folder to PATH
  export PATH="/Data1/software/ont-guppy/bin:$PATH"

# Note: current incompatibility exists between guppy 6.0 and Megalodon. Recommend using version 5.0
# wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_6.0.6_linux64.tar.gz
# tar -xf ont-guppy_6.0.1_linuxaarch64_cuda10.tar.gz
# pip install ont-pyguppy-client-lib==6.0.1
fi

#install samtools
printf "Y" | sudo apt install samtools

#install bedtools
sudo apt-get install bedtools

#install winnowmap
if [ ! -f /home/ubuntu/Data1/software/v2.03.tar.gz ]; then
  wget https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz
  tar zxvf v2.03.tar.gz
  cd Winnowmap-2.03
  make -j8
  # cd ..
  # mv Winnowmap-2.03 Winnowmap
fi

#install/update megalodon
sudo pip install megalodon

#Install minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.20/minimap2-2.20_x64-linux.tar.bz2 | tar -jxvf -

### IF HAVING TROUBLE WITH NVIDIA DRIVERS
# display ont-guppy requirements
#apt-cache show ont-guppy

# Guppy requires Cuda 11.1.1
#wget https://developer.download.nvidia.com/compute/cuda/11.1.1/local_installers/cuda_11.1.1_455.32.00_linux.run
#sudo sh cuda_11.1.1_455.32.00_linux.run
#accept and install driver 455 & CUDA toolkit 11.1.1
#sudo apt install nvidia-cuda-toolkit
#export PATH="/usr/local/cuda-11.1/bin"

# -   LD_LIBRARY_PATH includes /usr/local/cuda-11.1/lib64, or, add /usr/local/cuda-11.1/lib64 to /etc/ld.so.conf and run ldconfig as root

### Install nvidia driver and CUDA
#wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
#sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
#wget https://developer.download.nvidia.com/compute/cuda/11.6.1/local_installers/cuda-repo-ubuntu1804-11-6-local_11.6.1-510.47.03-1_amd64.deb
#sudo dpkg -i cuda-repo-ubuntu1804-11-6-local_11.6.1-510.47.03-1_amd64.deb
#sudo apt-key add /var/cuda-repo-ubuntu1804-11-6-local/7fa2af80.pub
#sudo apt-get update
#sudo apt-get -y install cuda

#export PATH=/usr/local/cuda-11.6/bin${PATH:+:${PATH}}
#export LD_LIBRARY_PATH=/usr/local/cuda-11.6/lib64\ ${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}

# Install methylartist for analysis ad visualization
pip install methylartist
pip install methplotlib
pip install methplotlib --upgrade

#install mysql needed to query chromosome sizes for methylartist violin plots
sudo apt install mysql-client-core-5.7
######

###### Download required files
# obtain the rerio model file
if [ ! -d /Data1/software/rerio ]; then
  git clone https://github.com/nanoporetech/rerio /Data1/software/rerio
#  /Data1/software/rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001
#  sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001 /opt/ont/guppy/data/
#  sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.cfg /opt/ont/guppy/data/
#  sudo cp /Data1/software/rerio/basecall_models/res_dna_r941_min_modbases-all-context_v001.jsn /opt/ont/guppy/data/
#  sudo cp -r /Data1/software/rerio/basecall_models/barcoding/* /opt/ont/guppy/data/barcoding/
  /Data1/software/rerio/download_model.py /Data1/software/rerio/basecall_models/*
else
  echo "rerio folder aready exists"
fi

# Download reference genome
# Update below for your specific organism.
if [ ! -f /Data1/reference/c_elegans.WS235.genomic.fa ]; then
  mkdir /Data1/reference
  wget https://downloads.wormbase.org/species/c_elegans/sequence/genomic/c_elegans.WS235.genomic.fa.gz -O /Data1/reference/c_elegans.WS235.genomic.fa.gz
  gunzip /Data1/reference/c_elegans.WS235.genomic.fa.gz
else
  echo "reference genome "reference/c_elegans.WS235.genomic.fa" already exists"
fi

#index using minimap
if [ ! -f /Data1/reference/ws235.mmi ]; then
  /Data1/software/minimap2-2.20_x64-linux/minimap2 -d /Data1/reference/ws235.mmi /Data1/reference/c_elegans.WS235.genomic.fa
  else
  echo "indexed reference already exists"
fi

#download nanopore fast5 files
#note: replace s3 bucket with where your sequencing data is stored.
if [ ! -d /Data1/seq_data/210614_Raja ]; then
  aws s3 sync "s3://nanopore-1/nanopore first run/" /Data1/seq_data/
else
  echo "210614_Raja nanopore data already synced from S3."
fi
######
