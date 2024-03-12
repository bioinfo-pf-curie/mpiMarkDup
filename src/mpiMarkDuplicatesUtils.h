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
     mpiMarkDuplicatesUtils.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#ifndef MPIMD_UTILS_H
#define MPIMD_UTILS_H 

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "log.h"
#include "reads.h"

char digit2char(int n) ;
int fastqToPhred (char ch) ;
int getLine(char** offset, char** tokenLine);
int getTokenTab(char **offset, char **tokenTab) ;
char* strapp(char* dest, char* src ) ;
unsigned int readBits(unsigned int x, int k) ;
void writeBits(unsigned int *x, int k) ;
void toggleBits(unsigned int *x, int k) ;
void switchBits(unsigned int *x, int k, int j);
#endif /* ifndef MPIMD_UTILS_H */
