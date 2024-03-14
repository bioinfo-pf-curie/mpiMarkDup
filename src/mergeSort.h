/*
   This file is part of mpiSORT
   
   Copyright Institut Curie 2022
   
   This software is a computer program whose purpose is to sort SAM file.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     bufferized_read.c

  Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/


#ifndef MERGESORT_H
#define MERGESORT_H


#include <stdlib.h>
#include <stdbool.h>
#include "parser.h"

#define SORT_TYPE_MS size_t
#define SORT_CMP_MS compare_size_MS

extern size_t *base_arr2;

/* structure to represent ranges within the array */
typedef struct {
        SORT_TYPE_MS start;
        SORT_TYPE_MS end;
} Range;

SORT_TYPE_MS Range_length(Range range);
Range Range_new(const SORT_TYPE_MS start, const SORT_TYPE_MS end);
Read* mergeSort(Read* c, size_t n);
Read* structMerge(Read* c, size_t p, Read* d, size_t q);

int compare_size_MS(const SORT_TYPE_MS *a, const SORT_TYPE_MS *b);
void InsertionSort(SORT_TYPE_MS *array, const Range range);
/* find the index of the last value within the range that is equal to array[index], plus 1 */
size_t BinaryLast(const SORT_TYPE_MS *array, const SORT_TYPE_MS value, const Range range);
void MergeSortRec(SORT_TYPE_MS *array, const Range range, SORT_TYPE_MS *buffer);
void MergeSortMain(SORT_TYPE_MS *array, const SORT_TYPE_MS array_count);

#endif


