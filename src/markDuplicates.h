/*
   This file is part of mpiMarkDup
   
   Copyright Institut Curie 2020
   
   This software is a computer program whose purpose is to sort SAM file and mark duplicates.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/


/*
   Module:
     markDuplicates.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#ifndef MARKDUPLICATES_H
#define MARKDUPLICATES_H

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include <mpi.h>

#include "mpiMarkDuplicatesUtils.h"
#include "perfectHash.h"
#include "reads.h"
#include "parser.h"
#include "llist.h"
#include "createLBList.h"

typedef struct {
    unsigned int firstCoord;
    unsigned int lastCoord;
} Interval;


int areComparableForDuplicates(readInfo *pair1, readInfo *pair2, const int isPaired);
void markDuplicateFragments(llist_t *cluster, int *totalDuplica, const int containsPairs);
void markDuplicatePairs(llist_t *cluster, hashTable *htbl, int *totalDuplica, int *totalOpticalDuplicate);
int countExternalMateInList(llist_t *list, size_t *externalMate);

readInfo *readParsing (char *sam_buff, 
                        Interval *intervalByProc, 
                        size_t readIndex,
                        size_t counter, 
                        chrInfo *chr, 
                        lbInfo *lb,
                        size_t *buffer_offset_source,  
                        MPI_Comm comm);

void insertReadInList(llist_t *l, readInfo *read);
int ComputePercentageIncrement(size_t readNum, int readSecRatio, int second);
char* addDuplicateFlag(char *oldReadLine, readInfo *read);
void checkMateBuffer(readInfo *read, Interval *intervalByProc, MPI_Comm comm);
//readInfo* cloneRead(readInfo *read);
int cmpMateRank(const void *a, const void *b);
int computeCountAndDispl(int *scounts, int **sdispls, int **rcounts, int **rdispls, MPI_Comm comm);
coord computePhysicalLocation(readInfo *read);
int countExternalMateInArray(readInfo **readArr, size_t readNum, unsigned int checkWithBruck);
void createCMType(MPI_Datatype *CMtype);
void createMateType(MPI_Datatype *mate_type);
void createMateType2(MPI_Datatype *mate_type);
void ensureMateRank(readInfo **readArr, Interval *intervalByProc, size_t readNum, MPI_Comm comm);
int exchangeAndFillMate(readInfo ***matesByProc, mateInfo *mates, size_t numberOfReadsToSend, MPI_Comm comm);
int fillReadAndFictitiousMate(readInfo **readArr, readInfo ***readArrWithExternal, size_t readNum);
int fillReadLBValue(char *token, readInfo *read);
void fillUnclippedCoord(readInfo *read);
char* flag2string(int flagValue);
Interval* gatherIntervalByProc(Interval interval, MPI_Comm comm);
Interval getIntervalFromBuffer(char *bufferReads, size_t readNum);
size_t getMateRankReadSizeBeforeBruck(llist_t *list, mateInfo **mates) ;
int getPosFromLine(char *samLine);
char* getReadTagValue(char *token, const char *tag);

int parseLibraries(char *bufferReads, 
                    Interval *intervalByProc, 
                    llist_t *fragList, 
                    llist_t *readEndsList, 
                    readInfo ***readArr, 
                    char ***samTokenLines, 
                    size_t readNum, 
                    size_t readIndex, 
                    chrInfo *chr, 
                    lbInfo *lb, 
                    size_t *reads_buffer_offset_source,
                    MPI_Comm comm);


char *markDuplicate (char *bufferReads, 
                    size_t readNum, 
                    char *header, 
                    MPI_Comm comm, 
                    size_t *buffer_read_offset_source,
                    size_t *all_discord_dup_offset, 
                    size_t *total_disc_dup);

int markMateDuplicateFlag(hashTable *htbl, readInfo *read, int d);
char* writeBuff(char **samTokenLines, readInfo **readArr, size_t readNum);
void zeroCopyBruck(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
CMInfo read2CM(readInfo *read) ;
size_t exchangeAndFillMate_with_Bruck(readInfo ***matesByProc, mateInfo *mates, size_t numberOfReadsToSend, MPI_Comm comm);
size_t exchangeAndFillMate_with_Bruck_v2(readInfo ***matesByProc, mateInfo *mates, size_t numberOfReadsToSend, MPI_Comm comm);
readInfo *buildReadEnds(readInfo *read1, readInfo *read2, llist_t *readEndsList);
readInfo *buildReadEnds_v2(readInfo *read1, readInfo *read2, llist_t *readEndsList);

#endif
