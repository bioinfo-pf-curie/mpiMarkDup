# Documentation

* [Installation](#installation)
* [Prerequisites](#prerequisites)
* [Usage](#usage)
    * [Input](#input)
    * [Options](#options)
    * [Output](#output)
* [Informatic resources](#informatic-resources)
    * [Memory](#memory)
    * [Cpu](#cpu)
* [Examples](#examples)
    * [Standard](#standard)
    * [Slurm](#slurm)
    * [PBS/Torque](#pbstorque)
* [Performance](#performance)
* [Parallel filesystems](#parallel-filesystems)
* [Algorithm](#algorithm)
* [References](#references)


## Installation

Follow the [Installation](INSTALL.md) guidelines to compile and install `mpiMarkDup`.

## Prerequisites

As `mpiMarkDup` relies on the Message Passing Interface (MPI) standard, `mpirun` must be available to run the program. Several MPI implementations exist
such as [mpich](https://www.mpich.org/), [open-mpi](https://www.open-mpi.org/) or [Intel® MPI Library](https://software.intel.com/en-us/mpi-library). The `mpirun` program must be available in your PATH.

### Install mpirun on CentOS

* for [open-mpi](https://www.open-mpi.org/): `sudo yum install openmpi openmpi-devel`
* for [mpich](https://www.mpich.org/): `sudo yum install mpich`

### Install mpirun on Ubuntu

* for [open-mpi](https://www.open-mpi.org/): `sudo apt-get install openmpi-bin`
* for [mpich](https://www.mpich.org/): `sudo apt-get install mpich`



## Usage

`mpiMarkDup` is executed with the `mpirun` program, for example:

`mpirun -n 4 mpiMarkDup examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiSORTExample -q 0 -d 1000 -v 4`

The `-n` options passed to `mpirun` indicates the number of processes to run in parallel (this is basically the number of cores that will be used). For for details on how to choose the number processes, see the [Informatic resources](#informatic-resources) section.

`mpiMarkDup` requires two mandatory arguments:

* [Input](#input): the path to the SAM file to be sorted and marked
* [Output](#output): the directory in which the results will be written

### Input

A SAM file produced by an aligner (such as [BWA](https://github.com/lh3/bwa)) with paired reads and compliant with the SAM format. Only paired reads are accepted so far.


### Options

* `-q INTEGER` filters the reads according to their quality. Reads quality under the threshold are ignored in the sorting results. Default is 0 (all reads are kept).
* `-d INTEGER` is the optical distance for duplicates.
* `-v INTEGER` is the level of the verbose:

	- 0 is `LOG_OFF`
	- 1 is `LOG_ERROR`
	- 2 is `LOG_WARNING`
	- 3 is `LOG_INFO (default)`
	- 4 is `LOG_DEBUG`
	- 5 is `LOG_TRACE` 

### Output

The output consists of gz files:
* one per chromosome (e.g. chr11.gz)
* one for discordant reads (discordant.gz): discordants reads are reads where one pair aligns on a chromosome and the other pair aligns on another chromosome
* one for unmapped reads (unmapped.gz): unmapped reads are reads without coordinates on any chromosome

Only discordant and chromosome files are marked.

To index the SAM:
`tabix -p sam chr11.gz`

`tabix` is available from [Samtools](http://www.htslib.org/doc/tabix.html)

To uncompress:
`bgzip -d chr11.gz -c > chr11.sam`

`bgzip` is available from [Samtools](http://www.htslib.org/doc/bgzip.html)

## Informatic resources

### Memory

Like for the sorting the total memory used during the marking is around one and a half the size of the SAM file.
For example, to sort a 1.3TB SAM file (such as the [NA24631](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/) from [GIAB](https://github.com/genome-in-a-bottle/about_GIAB) which is a 300X Whole Genome (2x150-base) paired reads that we aligned with [mpiBWA](https://github.com/bioinfo-pf-curie/mpiBWA)), 1.7 TB of memory are required and splitted over 512 MPI workers (i.e. cores) that corresponds to makes 3.3 Gb of memory per core.

NA24631 sample is available here: ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams

### Cpu

Due to the bitonic sorting, the algorithm is optimized for power of 2 number of CPU. Therefore, it is recommended to set the `-n` parameter of `mpirun` to 2, 4, 8, 16, 32, etc. in order to ensure for optimal performance. For example, `mpirun -n 4 mpiMarkDup examples/data/HCC1187C_70K_READS.sam -q 0 -d 0 -v 4 ${HOME}/mpiMarkDupExample`


However, the `-n` parameter can be set to any other value but extra  MPI communications will be added to fit power of 2 required by the bitonic algorithm. In this case, additonal memory is needed for the rank 0 worker. This rank is responsible for collecting and dispatching the data before and after bitonic.


## Examples

There are many ways to distribute and bind MPI jobs according to your architecture. We provide below several examples to launch MPI jobs in a standard manner or with a job scheduling system such as [Slurm](https://slurm.schedmd.com/sbatch.html) and [PBS/Torque](https://support.adaptivecomputing.com/support/documentation-index/torque-resource-manager-documentation/).

A toy dataset (SAM file) is provided in the [examples/data](../examples/data) folder for testing the program. We test it with from 1 to 8 jobs and 2 nodes with a normal network and file system. You can use this sample to test the program.

### Standard

`mpirun` can be launched in a standard manner without using any job scheduling systems. For example:

`mpirun -n 4 mpiMarkDup examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiMarkDupExample -q 0 -d 1000 -v 4`

If needed, a file with the server name in `-host` option can be provided to `mpirun`. We invite you to read the `mpirun` documentation for more details.


### Slurm

```shell
#! /bin/bash
#SBATCH -J MPIMARKDUP_MYSAM_4_JOBS
#SBATCH -N 2                       	# Ask 2 nodes
#SBATCH -n 4                       	# Total number of cores
#SBATCH -c 1			   	# use 1 core per mpi job
#SBATCH --tasks-per-node=2         	# Ask 2 cores per node
#SBATCH --mem-per-cpu=${MEM}	   	# See Memory ressources
#SBATCH -t 01:00:00
#SBATCH -o STDOUT_FILE.%j.o
#SBATCH -e STDERR_FILE.%j.e

mpirun mpiMarkDup examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiMarkDupExample -q 0 -d 1000 -v 4

```

### PBS/Torque

```shell
#! /bin/bash
#PBS -N MPIMARKDUP_MYSAM_4_JOBS
#PBS -l	nodes=2:ppn=2:mem=${MEM}     	# Ask 2 nodes and 2 jobs per node
#PBS -l walltime=24:00:00
#PBS -o STDOUT_FILE.%j.o
#PBS -e STDERR_FILE.%j.e

mpirun mpiMarkDup examples/data/HCC1187C_70K_READS.sam ${HOME}/mpiMarkDupExample -q 0 -d 1000 -v 4

```

## Parallel filesystems

As the `mpiMarkDup` program uses MPI functions for reading and writing you can take advantage of a parallel file system to tackle the IOs bottleneck and speed-up the sorting and the marking of a SAM file. Using a parallel filesystem such as [Lustre](http://lustre.org/) or [BeeGFS](https://www.beegfs.io/) will be mandatory as long as the SAM file to sort is bigger and bigger.


Note that for advanced users, it is even possible to fine tune the source code such that `mpiMarkDup` can optimally used the settings of the parallel filesystem on which you will run the program. Refer to the [Advanced](ADVANCED.md) guidelines to do so.


We don't recommend to use MPI with [NFS](https://en.wikipedia.org/wiki/Network_File_System) (it works but it does not scale very well).

## Performance

(In construction)


## Algorithm

Firstly, the program sorts the reads by genome position and extract discordant and unmapped (end and mate) reads with the same method described in [mpiSORT](https://github.com/bioinfo-pf-curie/mpiSORT).

Secondly, the program marks the duplicated reads for each chromosome and discordant reads according to Picard Markduplicate method. The unmapped and unmapped mates are not marked. To limit memory overhead, we build a distributed perfect hash table (see [perfectHash.c](../src/perfectHash.c) for details) for fragment list and end list. This way the memory usage is under the memory usage of [mpiSORT](https://github.com/bioinfo-pf-curie/mpiSORT).

Finally, each chromosome is marked and compressed with bgzf and written down in the output folder.

We test the reproducibility by comparing both pipelines: mpiMarkDup and [mpiSORT](https://github.com/bioinfo-pf-curie/mpiSORT) + Picard (MarkDuplicate).
We use the same number of cpu for each pipeline and we obtained the same results.

If the number of cpu differs, then the reproducibility is not garanted. Indeed, tie cases are solved using the index of the read in the sorted file. This index can differ with the number of cpu.

This problem does not impact the results in the downstream analysis.

In conclusion when you test reproducibility always take the same number of cpu.


