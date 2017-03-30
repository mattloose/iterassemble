# Virtual Genome Walking
VGW is an iterative assembler that will locally assemble genome reads based on existing transcripts. This is designed for hard to assemble genomes where WGS assembly approaches fail.

## Installation

To run on all platforms using docker (may limit RAM or CPU use):

```
docker run -v <host drive>:/data <name> assemble.py -h
docker run -v <host drive>:/data <name> parse_output.py -h
docker run -v <host drive>:/data <name> summarise_vgw.py -h
```

To run on UNIX platforms without any limitations:

```
git clone https://github.com/mattloose/iterassemble
cd iterassemble/
make
```

Prior to running VGW there are several dependencies required, which can either be installed manually or using conda:


#### 1. Conda
  
  Download and install fermi-lite (https://github.com/lh3/fermi-lite), ensure program is in the PATH. Then simply use conda to install the remaining packages, either to the root:
  ```
  conda install -c bioconda -c auto -c biobuilds --file conda_install.txt
  ```
  or within a virtual environment:
  ```
  conda create -n vgwEnv
  source activate vgwEnv
  conda install -c bioconda -c auto -c biobuilds --file conda_install.txt
  ```
  
#### 2. Manual

  Spend a long time downloading and installing the following:
  * BWA
  * Tophat
  * SOAPdenovo2
  * cap3
  * fermi-lite 
  * gawk
  * BLAST+
  * BioPython
  * multiprocessing
  * GNU parallel
  * mafft
  * samtools
  * bedtools
  * gmap

