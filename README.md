# Virtual Genome Walking
VGW is an iterative assembler that will locally assemble genome reads based on existing transcripts. This is designed for hard to assemble genomes where WGS assembly approaches fail.

## Installation

#### 1. Docker

The docker image contains all dependencies and can be run on all operating systems. Docker itself may be limited on CPU/MEM usage. To install and run using docker:

```
docker run terievans/docker-vgw 
```

#### 2. UNIX

To download onto UNIX platforms:

```
git clone https://github.com/mattloose/iterassemble
cd iterassemble/
make
sudo make install
```
This will install the C/C++ programs, we also recommend you add the `iterassemble/` directory to the PATH.

## Dependencies 

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
```
assemble.py -h 
usage: assemble.py [options] cDNA.fa in1.fq in2.fq

Run an iterative "genome walking" assembly

positional arguments:
  cDNA.fa               cDNA/mRNA sequence file
  in1.fq                FASTQ file of first paired-end reads
  in2.fq                FASTQ file of second paired-end reads

optional arguments:
  -h, --help            show this help message and exit
  -t [INT]              Number of processors to use (default: 10)
  -m [INT]              Maximum number of iterations to make (default: 50)
  -d [dir/]             Directory for fastq indexes, will be overwritten
                        (default: Split_files/)
  -i [INT], --insert [INT]
                        Average insert size (default: 200)
  -l [INT], --length [INT]
                        Maximum read length (default: 100)
  -e [INT], --endsize [INT]
                        Number of bases from each end of the contigs to map
                        (default: 600)
  -r [INT], --remove [INT]
                        After 10 iterations, remove contigs shorter than INT
                        (default: 200)
  -F [INT], --fastmap1 [INT]
                        Minimum SMEM length permited in bwa fastmap for first
                        iteration (default: 40)
  -f [INT], --fastmap [INT]
                        Minimum SMEM length permited in bwa fastmap for later
                        iterations (default: 60)
  -c [INT], --culling [INT]
                        After each iteration blast step uses this
                        culling_limit (default: 2)
  -n [INT], --nreads [INT]
                        Contigs with n reads mapped across are kept after each
                        iteration (default: 5)
  -M [INT], --maxcontigs [INT]
                        Maximum number of contigs per gene after each
                        iteration (default: 500)
  -s [INT], --split [INT]
                        FASTQ files will be split on this line number, must be
                        divisible by 4 (default: 4000000)
  --continue_from [INT]
                        Continue iterating from this number (default: 0)
  --soap [SOAP]         SOAP denovo program to use, must be in path (default:
                        SOAPdenovo-63mer)
  --end_process_only    No iterative assembly will be performed, just the end
                        process based on existing files (default: False)
  ```

The second script (`parse_output.py`) will compare the genome walked output (`Final_sequences.fasta`) against the original transcripts and remove any contigs that are incorrectly incorporated. 

```
parse_output.py -h
usage: parse_output.py [-h] [-o] [-g [GFF]] genome transcripts

Parse vgw output using gmap

positional arguments:
  genome                vgw output file
  transcripts           transcript sequences

optional arguments:
  -h, --help            show this help message and exit
  -o, --overwrite       overwrite gmap database
  -g [GFF], --gff [GFF]
                        intermediate gff filename, default (parsing.gff)
```

The output from this program (`Final_sequences.fasta_gmapparsed`) can then be processed by `summarise_vgw.py` to generate a gff file and assembly metrics. 
```
summarise_vgw.py -h
usage: summarise_vgw.py [-h] [-o] [-g [GFF3]] [-s [SUMMARY]]
                        genome transcripts

Summarise vgw results using gmap

positional arguments:
  genome                vgw output file
  transcripts           transcript sequences

optional arguments:
  -h, --help            show this help message and exit
  -o, --overwrite       overwrite gmap database
  -g [GFF3], --gff3 [GFF3]
                        output gff3 filename, default (vgw.gff)
  -s [SUMMARY], --summary [SUMMARY]
                        output file name, default (vgw_summary.txt)
```
