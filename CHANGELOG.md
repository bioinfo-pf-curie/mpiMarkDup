


Release 1.0 from the 25/10/2019 <br />

1) Not a good idea to play with magic number bgzip complain. We come back to previous version

2) TODO: now the major bottle neck is the computation of the collision in perfect hash table. <br />
		 we shall compute collision only one time according to the biggest chromosom to mark <br /> 

3) TODO: problem with the option -q > 0 (it works with -q 0)

Release 1.0 from the 23/10/2019 <br />

1) Add the magic number at end of BGZF file. This way samtools doesn't complain <br />

example of usage with samtools
- samtools flagstat chr1.gz 
- samtools view -f 1024 chr1.gz 

Other tools shall complain (not tested). <br />

Release 1.0 from the 21/10/2019 <br />

1) Review of the search of mates after the bruck exchange of inter-process reads. <br />
     We now use the hash table to look for the corresponding pair of a read.  <br />


Release 1.0 from the 02/07/2019 <br />

1) Fix a bug when pairs overlap 2 ranks <br />

2) Fix a memory leak <br />

Release 1.0 from the 01/07/2019 <br />

1) Fix memory leaks <br />

2) Remove global variable COMM_WORLD <br />

Release 1.0 from the 27/06/2019 <br />

1) Fix a corner case. We check intra-node overlap before extra-node overlap in mpiSort.c. <br />

2) Fix a bug in sort_any_dim.c  <br />

Release 1.0 from the 19/06/2019 <br />

1) Fix a corner case <br />

Release 1.0 from the 13/06/2019 <br />

1) Replace All2all with a Bruck when exchange external frags <br />

2) Change some types in readInfo and mateInfo structures <br />

Release 1.0 from the 23/05/2019 <br />

1) Fix issues with openssl 1.0.2k-fips. <br />

Release 1.0 from the 15/05/2019 <br />

1) Cleaning up the code. <br />

Release 1.0 from the 07/05/2019 <br />

1) Fix reproducibility issue with 1 cpu. <br />

Release 1.0 from the 06/05/2019 <br />

1) Fix tie case (tested only with power of 2 cpu). <br />

Release 1.0 from the 08/04/2019 <br />

1) fix reproducibility issue but still have a corner case.  <br />

Release 1.0 from the 31/03/2019 <br />

1) fix corner case when reads distribution is unbalanced. <br />

Release 1.0 from the 28/03/2019 <br />

1) fix some integer conversion and prototype. <br />
2) fix a corner case when a rank recieve no mate from other rank. <br />

Release 1.0 from the 22/03/2019 <br />

1) fix the coordinates overlapped between 2 ranks in case we are not power of 2 dimension (again it's better if you are in power of 2 dimension) <br />

