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
     perfectHash.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmin Martin,      Paris Descartes University
*/

#ifndef PERFECT_HASH_H
#define PERFECT_HASH_H


/**
 *   @date 2018 Feb 26
 *   @file perfectHash.h
 *   @author Firmin Martin
 *   @brief Perfect hashing implementation adapted to multiple processes context.
 *   @details Each process implement a perfect hashing table with a common main universal hash function. 
 *            Hence, processes can exchange read which have a locally given fingerprint.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

// use -lcrypto to compile
// TODO:should be replaced by another open project
#include <openssl/md5.h>
#include <mpi.h>

#include "reads.h"
#include "mpiMarkDuplicatesUtils.h"
#include "log.h"

typedef struct {
    size_t m;
    int a;
    int b;
} hashParam;

typedef struct {
    hashParam h;
    readInfo **table;
} secTable;

typedef struct {
    hashParam h;
    size_t prime;
    size_t size;
    secTable **secTable;
} hashTable;

// Utils
size_t computeSpace(hashParam hp, int p, readInfo **arr, size_t size);
int haveCollision(hashParam hp, int prime, size_t *arr, size_t size);
void printPerfectHashTable(hashTable *htbl);

// Conversion
size_t read2mateFP(readInfo* read) ;
size_t read2Fingerprint(readInfo *read);
size_t string2MD5(const char *str);

// Construction 
void constructMainUnivHash(hashParam *hp, int *p, readInfo **arr, size_t size);
void constructMainUnivHashWithHp(hashParam hp, int *p, readInfo **arr, size_t size);
void constructSecTable(hashTable *htbl, readInfo **arr, size_t size);
void hashTableInit(hashTable *htbl, readInfo **arr, size_t size);
void hashTableInitWithHp(hashTable *htbl, hashParam hp, readInfo **arr, size_t size);
void shareHpAndConstructHtbl(hashTable *htbl, readInfo **arr, size_t size, MPI_Comm comm) ;

// get element
readInfo *getReadFromFingerprint(hashTable *htbl, size_t fingerprint);
readInfo *getReadFromQnameAndPairNum(hashTable *htbl, char *QName, int PairNum) ;
readInfo *getMateFromRead(hashTable *htbl, readInfo *read) ;

// Methods
void hashTableDestroy(hashTable *htbl);
void hashTableInsert(hashTable *htbl, readInfo *read) ;

// Hash function
int univHash(hashParam hp, int p, size_t k);

#endif
