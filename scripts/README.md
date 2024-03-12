# Scripts

## compareSAMWithExternalTools.sh

### Description

Compare two SAM files given in input. First the SAM files are converted into BAM
files by samtools view function. Then, BAM files are compared with
the function `cmpbams` of jvarkit. We assume that the SAM files given
in input are coherent. By default the script assumes that `samtools` and
`cmpbams.jar` are located in script directory. If there is no `-d` option,
it assumes that the script directory is the current directory.

It creates a directory containing samtools and
jvarkit.jar for temporary files.

### Options

* -g indicates that input are gzipped
* -d directory where the script are. It is where
     compareSAMWithExternalTools.sh creates/uses a temporary directory
* -s indicates where is `samtools`. By default the script assumes that `samtools` is located in the script directory.
* -j indicates where is `cmpbams.jar`. By default the script assumes that
     `cmpbams.jar` is located in the script directory.
* -h prints help


### How to used it

`compareSAMWithExternalTools.sh [-g] [-s samtools] [-j cmpbams.jar] -d DIRECTORY INPUT1 INPUT2`

where:
* DIRECTORY is the directory where the executable `samtools` and jvarkit.jar are.
* INPUT1 and INPUT2 are the inputs.

## compareMpiMarkDupVsMarkDuplicates.sh

### Description

Compare the output of mpiMarkDup and MarkDuplicates. It takes a directory
containing directories of class of files to test and to compare the
files. The script can create a very simple log file which contains the
output and the records that don’t match. It creates also a temporary
directory in the directory where the script is launched.

### Options

* -h for help
* -s indicates the path of samtools
* -p indicates the path of picard.jar
* -j indicates the path of jvarkit.jar
* -S indicates the path of compareSAMWithExternalTools.sh
* -m indicates the path of mpiMarkDup
* -l indicates that a log file is generated. You must indicate the name of the log file. Information are adding at the end of the file

### How to used it

```shell
compareMpiMarkDupVsMarkDuplicates.sh [-h] [-s samtools] [-p picard.jar] [-j
jvarkit.jar] [-S compareSAMWithExternalTools.sh] [-m mpiMD]
[-l logfile] DIRECTORY
```

where DIRECTORY is the class test directory path.

## cmpWithMD.sh

### Description

Given a sam, compare the results between `markDuplicates` and `mpiMarkDup`.

### Synopsis

`./cmpWithMD.sh SAM`

### User variables

User needs to set *user variables* at first use. They are

- PICARDTOOLS : picard tools jar path
- MPI_MD : `mpiMarkDup` executable path
- NUM_PROCS : number of processes for `mpiMD`
- JVARKIT : `cmpbams` path from jvarkit
- SAMTOOLS : `samtools` executable paths

### Details

- This script compares `markDuplicates` and `mpiMarkDup` by chromosome which appears in the header. Unmapped (only the case where read and its mate are unmapped) and discordants are removed from `markDuplicates` result. 
- This script givse some stats in the standard output.
- As output, a directory is created, containing `markDuplicates` output, `mpiMarkDup` output  and `jvarkit` output.

#### What this script doesn't do

- This script do not compare discordant.gz and unmapped.gz since they aren't handle yet.

## runCmpOnDir.sh

Given a directory of sam, compare for each sam the result of `markDuplicates` and `mpiMarkDup`. There is a logfile as output and directories created by cmpWithMD.

### Synopsis

`./runCmpOnDir.sh DIR`

