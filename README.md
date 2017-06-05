# Virtual Genome Walking
VGW is an iterative assembler that will locally assemble genome reads based on existing transcripts. This is designed for hard to assemble genomes where WGS assembly approaches fail.

VGW can either be used in Docker, or can be installed and maintained locally. Note that this is still under active development so frequent updates may be required.

## Docker

The docker image contains all dependencies and can be run on all operating systems. Docker itself may be limited on CPU/MEM usage. To install and run using docker:

```
docker run terievans/docker-vgw 
```

To mount a local drive that contains your input files and run the three main programs with docker use the following. In this example we are mounting the local drive into `/data` but this could be anything so long as it is consistant. All output files will appear in your local drive. For more information on what these programs do, see below. Note: Your local drive must be accessible by docker! 

```
docker run -v <local drive>:/data terievans/docker-vgw assemble.py /data/<transcripts.fa> /data/<genome1.fq> /data/<genome2.fq> --docker_vol /data
docker run -v <local drive>:/data terievans/docker-vgw parse_output.py /data/<VGW output> /data/<transcripts.fa> --docker_vol /data
docker run -v <local drive>:/data terievans/docker-vgw summarise_vgw.py /data/<parsed VGW output> /data/<transcripts.fa> --docker_vol /data
```

## UNIX Installation

To download onto UNIX platforms:

```
git clone https://github.com/mattloose/iterassemble
cd iterassemble/
make
sudo make install
```

This will install the C/C++ programs, you will also need to add the `iterassemble/` directory to the PATH.

### Dependencies 

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

## Running VGW

There are three key programs within VGW, the first (`assemble.py`) is the base program that will map, extend and assemble genome fragments. It requires a fasta file of transcript sequences and two paired-end genome read files. These genome read files must contain identical IDs either prior to a blankspace or '/1' or '/2'. For extremely large genomes we recommend first repeat-depleting the read set by k-mer analysis.

The longest process will likely be the inital parsing and indexing of the genome reads, the results of which will be output to a folder. So long as the vgw is run from the same location with the same sequences this index will not need to be remade.

We include axolotl genome reads that map to the BAC JF490016 as well as the axolotl nanog transcript in the example folder. To run VGW for 5 iterations on these data see below, for a complete list of running variables run `assemble.py -h`. 

```
assemble.py example/Transcript.fa example/NanogBAC_mappedpairs.1.fastq.gz example/NanogBAC_mappedpairs.2.fastq.gz -m 5
```

The second script (`parse_output.py`) will compare the genome walked output (`Final_sequences.fasta`) against the original transcripts and remove any contigs that are incorrectly incorporated. Again, for a full list of options run `parse_output.py -h`.  

```
parse_output.py example/Final_sequences.fasta example/Transcript.fa
```

The output from this program (`Final_sequences.fasta_gmapparsed`) can then be processed by `summarise_vgw.py` to generate a gff file and assembly metrics. 

```
summarise_vgw.py example/Final_sequences.fasta_gmapparsed example/Transcript.fa
```
