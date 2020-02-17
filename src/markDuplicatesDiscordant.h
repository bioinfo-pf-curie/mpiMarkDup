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
     markDuplicatesDiscordant.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#ifndef MARKDUPLICATESDISCORDANT_H
#define MARKDUPLICATESDISCORDANT_H

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
#include "markDuplicates.h"
 
readInfo *readParsingDiscordant (char *sam_buff, 
                                Interval *intervalByProc, 
                                size_t readIndex,
                                size_t counter, 
                                chrInfo *chr, 
                                lbInfo *lb,
                                size_t *disc_offsets_source,  
                                MPI_Comm comm); 

void exchangeExternFragDiscordant(llist_t *fragList, llist_t *readEndsList, hashTable *htbl, MPI_Comm comm);

size_t parseLibrariesDiscordant(char *bufferReads, 
                                Interval *intervalByProc, 
                                llist_t *fragList, 
                                llist_t *readEndsList, 
                                readInfo ***readArr, 
                                char ***samTokenLines, 
                                size_t readNum, 
                                size_t readIndex, 
                                chrInfo *chr, 
                                lbInfo *lb, 
                                size_t *disc_offsets_source,
                                MPI_Comm comm);

readInfo *buildReadEndsDiscordant(readInfo *read1, readInfo *read2, llist_t *readEndsList, int case_insert);

char *markDuplicateDiscordant (char *bufferReads, 
                                size_t readNum, 
                                char *header, 
                                MPI_Comm comm,
                                size_t *disc_offsets_source, 
                                size_t **disc_dup_offset_source, 
                                size_t *disc_dup_number);

void findDuplicaDiscordant(llist_t *fragList, 
                            llist_t *readEndsList, 
                            hashTable *htbl, 
                            int *totalDuplica, 
                            int *totalOpticalDuplicate, 
                            MPI_Comm comm);

#endif
