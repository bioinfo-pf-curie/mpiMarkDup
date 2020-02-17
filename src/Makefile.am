
bin_PROGRAMS = mpiMarkDup
noinst_HEADERS = $(mpiMarkDupHDR) $(mpiSortHDR)

mpiMarkDup_SOURCES = $(mpiMarkDupSRC) $(mpiSortSRC)

mpiMarkDupSRC = markDuplicates.c markDuplicatesDiscordant.c createLBList.c perfectHash.c log.c mpiMarkDuplicatesUtils.c reads.c llist.c

mpiMarkDupHDR = markDuplicates.h markDuplicatesDiscordant.h createLBList.h perfectHash.h log.h mpiMarkDuplicatesUtils.h common.h reads.h llist.h

mpiSortSRC = mpiSort.c mpiSortUtils.c mergeSort.c parallelBitonicSort.c \
	parallelBitonicSort2.c parallelBitonicSort3.c qksort.c parser.c \
	tokenizer.c preWrite.c write.c  writeUtils.c sortAnyDim.c  

mpiSortHDR = mpiSortUtils.h mergeSort.h parallelBitonicSort.h \
	parallelBitonicSort2.h parallelBitonicSort3.h qksort.h \
	parser.h tokenizer.h preWrite.h write.h writeUtils.h \
	bgzf.c bgzf.h khash.h sortAnyDim.h \
	compat.h 

AM_CFLAGS= -g -O2 -Wall -Wextra -std=gnu99 
#AM_CFLAGS= -g -O2
#AM_LDFLAGS=--coverage
