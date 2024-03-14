/*
   This file is part of mpiSORT
   
   Copyright Institut Curie 2020
   
   This software is a computer program whose purpose is to sort SAM file.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     mergeSort.c

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mergeSort.h"

Read* mergeSort(Read* c, size_t n){
	size_t q, p;
	Read* d;

	q = n / 2;
	p = n - q;

	if(p > 1){
		d = mergeSort(c, p);
		if(q > 1)
			mergeSort(d, q);
	}
	else
		d = c->next;
	d = structMerge(c, p, d, q);

	return d;
}

Read* structMerge(Read* c, size_t p, Read* d, size_t q){

	Read* t;

	while(1){

		if (c->next->coord > d->next->coord){
			t = d->next;
			d->next = t->next;
			t->next = c->next;
			c->next = t;
			if (q == 1)
				break;
			--q;
		}
		else {
			if(p == 1){
				while (q > 0){
					d = d->next;
					--q;
				}
				break;
			}
			--p;
		}
		c = c->next;
	}
	return d;
}

SORT_TYPE_MS Range_length(Range range) { return range.end - range.start; }

Range Range_new(const SORT_TYPE_MS start, const SORT_TYPE_MS end) {
        Range range;
        range.start = start;
        range.end = end;
        return range;
}

int compare_size_MS(const SORT_TYPE_MS *a, const SORT_TYPE_MS *b){
         size_t aa = *(size_t *)a, bb = *(size_t *)b;
         return (base_arr2[aa] < base_arr2[bb]);
}

/* Use during merge sort to stop the recursion */
void InsertionSort(SORT_TYPE_MS *array, const Range range) {
        size_t i, j;
        for (i = range.start + 1; i < range.end; i++) {
                const SORT_TYPE_MS temp = array[i];
                for (j = i; j > range.start && SORT_CMP_MS(&temp, &array[j - 1]); j--)
                        array[j] = array[j - 1];
                array[j] = temp;
        }
}

/* find the index of the last value within the range that is equal to array[index], plus 1 */
size_t BinaryLast(const SORT_TYPE_MS *array, const SORT_TYPE_MS value, const Range range) {
        size_t start = range.start, end = range.end - 1;
        if (range.start >= range.end) return range.end;
        while (start < end) {
                size_t mid = start + (end - start)/2;
                if (!SORT_CMP_MS(&value, &array[mid]))
                        start = mid + 1;
                else
                        end = mid;
        }
        if (start == range.end - 1 && !SORT_CMP_MS(&value, &array[start])) start++;
        return start;
}

void MergeSortRec(SORT_TYPE_MS *array, const Range range, SORT_TYPE_MS *buffer) {
        size_t mid, A_count = 0, B_count = 0, insert = 0;
        Range A, B;

	//use to stop recursion	
	if (Range_length(range) < 32){
		InsertionSort(array, range);
		return;
	}

	//chose a pivot
        mid = range.start + (range.end - range.start)/2;
        A = Range_new(range.start, mid);
        B = Range_new(mid, range.end);


	//start the recursion
        MergeSortRec(array, A, buffer);
        MergeSortRec(array, B, buffer);

        //merge part only A is copied 
        A = Range_new(BinaryLast(array, array[B.start], A), A.end);
        memcpy(&buffer[0], &array[A.start], Range_length(A) * sizeof(SORT_TYPE_MS));
        while (A_count < Range_length(A) && B_count < Range_length(B)) {
                if (!SORT_CMP_MS(&array[A.end + B_count], &buffer[A_count])) {
                        array[A.start + insert] = buffer[A_count];
                        A_count++;
                } else {
                        array[A.start + insert] = array[A.end + B_count];
                        B_count++;
                }
                insert++;
        }

        memcpy(&array[A.start + insert], &buffer[A_count], (Range_length(A) - A_count) * sizeof(SORT_TYPE_MS));
}



void MergeSortMain(SORT_TYPE_MS *array, const SORT_TYPE_MS array_count) {
        SORT_TYPE_MS *buffer = malloc(((array_count +1)/2) * sizeof(SORT_TYPE_MS));
	MergeSortRec(array, Range_new(0, array_count), buffer);
        free(buffer);
}

