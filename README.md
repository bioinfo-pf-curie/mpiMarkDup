# mpiMarkDup

This software allows the sorting and the marking of duplicates reads for high-throughput sequencing data after alignment from [SAM files](https://samtools.github.io/hts-specs/). `mpiMarkDup` relies on the [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface) standard to perform the parallelisation of the sorting processing over multiple cores and nodes of high performance computing clusters.

* [Installation](docs/INSTALL.md)
* [Documentation](docs/README.md)
* [Credits](#credits)
* [Acknowledgements](#acknowledgements)
* [Release notes](CHANGELOG.md)
* [Citation](#citation)

# Quick introduction

The programm sorts and markduplicates SAM files.
The SAM files must contain read group.
This programm has only been tested for short reads and upon whole genome.
The total memory needed is approximately 2 times the SAM file size.
It produces BGZF files for each chromosome present in the header and are
compatible with samtools.
The programm markduplicates discordant reads and optical duplicates.
Discordant reads are also marked in a separates BGZF file.
example of command lines:
mpirun numactl --interleave=all  $MPIMD $SAM $OUTPUT -q 0 -d 1000 -v 4
OUTPUT is the output directory \
d is for optical distance \
v is for log level \
q if for quality \

## Contacts

Contacts: [frederic.jarlier@curie.fr](mailto:frederic.jarlier@curie.fr])

