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
     mpiSortUtils.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifndef MPISORTUTILS_H
#define MPISORTUTILS_H


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include<time.h>
#include<limits.h>
/**************************************************
 *
 *      list interface to handle chromosoms name
 *      collision in the hash table
 *  
 **************************************************/
typedef struct Read Read;

struct Read
{
        /*
         * task: REMOVE STRING
         */
        char* string;
        size_t coord;
        size_t coord_mate;
        size_t offset_source_file;
        size_t offset;
        char quality;
        struct Read* next;
        struct Read* link;
}; //attribute packed slows by x2 times





typedef struct Elmt_ {
                int index; //hold the chromosom index
} Elmt;

/**************************************************
 *
 *  Hash function prototype for chromosoms string name
 *  
 **************************************************/
//for hg19
//#define            PRIME_TBLSIZ       65531

//for hg38
#define              PRIME_TBLSIZ       1049897

unsigned int hashpjw(const void *key);
/**************************************************
 *
 *  Hash function prototype for chromosoms string name
 *  
 **************************************************/


typedef struct CHTbl_ {
        unsigned int (*h)(const void *key);
        int buckets;
        int size;
        Elmt *elements;

} CHTbl;

void print_table(CHTbl *htbl);

int chtbl_init(CHTbl *htbl, int buckets, unsigned int (*h)(const void *key));

void chtbl_destroy(CHTbl *htbl);

int chtbl_insert(CHTbl *htbl, const void *data, int i);

int chtbl_remove(CHTbl *htbl, char *data);

int chtbl_lookup(const CHTbl *htbl, char *data);

void chtbl_destroy(CHTbl *htbl);

#define chtbl_size(htbl) ((htbl)->size)

void get_coordinates_and_offset_source_and_size_and_free_reads(
		int rank,
		int *local_read_rank,
		size_t *coordinates,
		size_t* offset,
		int* size,
		Read* data_chr,
		size_t local_readNum
		);

size_t init_coordinates_and_size(
		int rank,
		int *local_reads_rank,
		size_t *local_reads_index,
		size_t* coordinates,
		int* size,
		Read* data_chr,
		size_t local_readNum
		);


void chosen_split_rank_gather_size_t(
		MPI_Comm split_comm,
		int rank,
		int num_proc,
		int master,
		size_t size,
		size_t *size_per_jobs,
		size_t *start_size_per_job,
		size_t *all_data,
		size_t *data,
		size_t start_index
		);

#define SUCCESS 0
#define ERROR_CREATION 1 
#define ERROR_DESTROY 2
#define NOT_DELETED_LIST 3
#define EMPTY_LIST 4
#define ELEMENT_NOT_FOUND 6
#define ERROR_WALK 7
#define ERROR_CREATION_CAPA_SUPERIOR 8

#define INFINITE_POSITIVE INT_MAX
#define INFINITE_NEGATIVE INT_MIN

// Used for showing elements in test
#define INDETERMINATE (INT_MAX-1)

typedef struct skip_node_t {
	int key, level; 
	struct skip_node_t *prev, *next; 
	struct skip_node_t *down, *up; 
} skip_node; 

typedef struct skip_list_t {
	skip_node *head, *tail; 
	int nElems;
} skip_list;

int Create(skip_list *sl);
int Destroy(skip_list *sl);
int Insert(skip_list *sl, int elem);
int Delete(skip_list *sl, int elem);
int Length(skip_list sl, int *length);
int Search(skip_list sl, int elem, bool *found);
int Cost_Search(skip_list sl, int elem, int *cost);

// Funcions auxiliars
bool Exist(skip_list sl); 
int createNewUpperLayer(skip_list *sl); 
int deleteUpperLayer(skip_list *sl); 
int deleteOneElement(skip_list *sl, int elem); 


#endif
