/*
   This file is part of mpiMarkDup
   
   Copyright Institut Curie 2020
   
   This software is a computer program whose purpose is to sort SAM file and mark duplicates.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
i*/

/*
   Module:
     mpiSortUtils.c

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

/*
 * function used in mpiSort.c
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "log.h"
#include "parser.h"

/**************************************************
 *
 *  Has function for chromosoms string name
 *  
 **************************************************/
unsigned int hashpjw(const void *key){

    const char *ptr;
    unsigned int val = 0;
    ptr = key;
    while (*ptr != '\0'){

        unsigned int tmp;

        val = (val << 4) + (*ptr);
        tmp = val & 0xf0000000;

        if (tmp){
                val = val ^ (tmp >> 24);
                val = val ^ tmp;
            }
        ptr++;
    }

    //replace PRIME_TBLSIZ with the actual table size
    return val%PRIME_TBLSIZ;

}

/*
 *      Hash table implementation
 *
 */

int chtbl_init(CHTbl *htbl, int buckets, unsigned int (*h)(const void *key)) {

    Elmt *element;
    int i;
    if ((htbl->elements = (Elmt *)malloc(buckets * sizeof(Elmt))) == NULL)
           return -1;

    htbl->buckets = buckets;

    for (i = 0; i < htbl->buckets; i++){
           element = &htbl->elements[i];
           element->index = -1;
    }

    htbl->h = h;
    htbl->size = 0;
    return 0;
}


void chtbl_destroy(CHTbl *htbl) {

    free(htbl->elements);
    return;
}

int chtbl_insert(CHTbl *htbl, const void *data, int index) {

    int bucket, retval = 0;
    Elmt    *element;
    bucket = htbl->h(data) % htbl->buckets;
    element = &htbl->elements[bucket];
    if (element->index != -1){
        fprintf(stderr,"!!!!problem with %s during htable creation change the PRIME_TBLSIZ size in mpiSortUtils.h. We exit the program! \n",(char *)data );
        int res = MPI_Finalize();
        assert(res == MPI_SUCCESS);
        exit(2);

    }
    else{

        element->index=index;
    }
    htbl->size++;

    return retval;
}

int chtbl_lookup(const CHTbl *htbl, char *data) {

    Elmt    *element;
    int bucket;
    bucket = htbl->h(data) % htbl->buckets;
    element = &htbl->elements[bucket];
    if ( element->index != -1 ){
        return element->index;
    }
    else {
        return -1;
    }
}


void get_coordinates_and_offset_source_and_size_and_free_reads(int rank, int *local_read_rank, size_t *coordinates, size_t* offset, int* size, Read* data_chr, size_t local_readNum){

   	size_t j;
   	Read* chr = data_chr;
   	Read* tmp;

   	//we initialize offset source and size_source
   	for(j = 0; j < local_readNum; j++){
   		coordinates[j] = 0;
   		size[j] = 0;
   		offset[j] = 0;
   		local_read_rank[j]=rank;
   	}

   	//first we are going to read in the source file
   	for(j = 0; j < local_readNum; j++){
   		//offset is the read size
   		coordinates[j] = chr->coord;
   		offset[j] = chr->offset_source_file;
   		size[j] = (int)chr->offset; //read length
        tmp = chr;
        chr = chr->next;
        free(tmp);
   	}

}


size_t init_coordinates_and_size(int rank, int *local_reads_rank, size_t *local_reads_index,
		size_t* coordinates, int* size, Read* data_chr, size_t local_readNum)
{
	size_t dataSize = 0;
	size_t j;
	Read* chr = data_chr;

	//we initialize offset source and size_source
	for(j = 0; j < local_readNum; j++){
		size[j] = 0;
		coordinates[j] = 0;
	}

	//first we are going to read in the source file
	for(j = 0; j < local_readNum; j++){
		//offset is the read size
		coordinates[j] = chr->coord;
		size[j] = (int)chr->offset; //read length
		local_reads_rank[j] = rank;
		local_reads_index[j] = j;
		dataSize += chr->offset;

		chr = chr->next;
	}

	return dataSize;
}

size_t init_coordinates_and_size2(int rank, int *local_reads_rank,
		size_t* coordinates, int* size, Read* data_chr, size_t local_readNum)
{
	size_t dataSize = 0;
	size_t j;
	Read* chr = data_chr;

	//we initialize offset source and size_source
	for(j = 0; j < local_readNum; j++){
		size[j] = 0;
		coordinates[j] = 0;
	}

	//first we are going to read in the source file
	for(j = 0; j < local_readNum; j++){
		//offset is the read size
		coordinates[j] = chr->coord;
		size[j] = (int)chr->offset; //read length
		local_reads_rank[j] = rank;
		dataSize += chr->offset;

		chr = chr->next;
	}

	return dataSize;
}


void chosen_split_rank_gather_size_t(MPI_Comm split_comm, int rank, int num_proc, int master, size_t size, size_t *size_per_jobs,
		size_t *start_size_per_job, size_t *all_data, size_t *data, size_t start_index)
{
	MPI_Status status;
	size_t k;
	int j2;
	//size_t *temp_buf = malloc(sizeof(size_t));


	if (rank == master){
		// we copy element for rank master_2
		// we eliminate zeros from the
		// beginning of the vector
		size_t st = start_size_per_job[master];
		for (k = 0; k < size_per_jobs[master]; k++){
			all_data[st] = data[k+start_index];
			st++;
		}

		for(j2 = 0; j2 < num_proc; j2++){

			if (j2 != master && size_per_jobs[j2] != 0){

				/*
				temp_buf  = realloc(temp_buf, size_per_jobs[j]*sizeof(size_t));
				assert( temp_buf !=0);
				*/

				size_t temp_buf[size_per_jobs[j2]];
				temp_buf[size_per_jobs[j2]] = 0;
				assert(temp_buf !=0 );
				MPI_Recv(temp_buf, size_per_jobs[j2], MPI_LONG_LONG_INT, j2, 0, split_comm, &status);

				size_t st = start_size_per_job[j2];
				for (k = 0; k < size_per_jobs[j2]; k++){
					all_data[st] = temp_buf[k];
					st++;
				}

			}
		}
	}
	else{
		MPI_Send(data, size, MPI_LONG_LONG_INT, master,  0, split_comm);
	}

}

//FROM https://github.com/yago-123/skip-list-c
int Create(skip_list *sl) {
	srand(time(0)); 
	sl->nElems = 0; 

	sl->head = (skip_node*)malloc(sizeof(skip_node)); 
	sl->head->level = -1; 

	sl->tail = (skip_node*)malloc(sizeof(skip_node)); 
	sl->tail->level = -1; 

	if(createNewUpperLayer(&(*sl)) == SUCCESS) {
		return SUCCESS; 
	} else {
		return ERROR_CREATION; 
	}
}

int Destroy(skip_list *sl) {
	bool end; 
	int i, j, maxLevel;
	skip_node *tmp, *tmp2, *tmp3; 

	if(Exist(*sl)) {
		maxLevel = sl->head->level+1; 
		tmp = sl->head; 
		while(tmp->level != 0) {
			tmp = tmp->down; 
		}

		for(j = 0; j < (sl->nElems+2); j++) {
			tmp2 = tmp->next; 
			end = false; 
			i = 0; 
			while(!end && i < maxLevel) {
				if(tmp != NULL) {
					tmp3 = tmp; 
					tmp = tmp->up;
					free(tmp3);  
				} else 	end = true; 
				i++; 
			} 
			tmp = tmp2; 
		}

		sl->nElems = -1;
		return SUCCESS;
	} else {
		return ERROR_DESTROY;  
	}
}

int Insert(skip_list *sl, int elem) {
	int i; 
	bool newLevel; 
	skip_node *node, *tmp, *tmp2; 

	node = sl->head; 
	while(node->level != 0) {
		if(node->key <= elem && node->next->key <= elem) {
			node = node->next;
		} else {
			node = node->down; 	
		}
	} 
	while(elem > node->key) {
		node = node->next; 
	}

	tmp = (skip_node*)malloc(sizeof(skip_node)); 
	tmp->key = elem; 
	tmp->level = 0;
	tmp->up = NULL;  
	tmp->down = NULL; 
	tmp->next = node; 
	tmp->prev = node->prev;
	node->prev->next = tmp;
	node->prev = tmp;
	node = tmp->prev;  
	i = 1;  
	while((rand()%2) == 0) {
		newLevel = false; 
		while(node->up == NULL && !newLevel) {
			if((node->key == INFINITE_NEGATIVE) && (node->level == sl->head->level)) newLevel = true; 
			else node = node->prev;
		}
		if(newLevel) {
			createNewUpperLayer(&(*sl)); 
			node = sl->head; 
		} else node = node->up;

		tmp->up = (skip_node*)malloc(sizeof(skip_node)); 
		tmp->up->key = elem; 
		tmp->up->level = i;
		tmp2 = tmp;  
		tmp = tmp->up;
		tmp->down = tmp2;  
		tmp->prev = node; 
		tmp->next = node->next;
		node->next->prev = tmp; 
		node->next = tmp; 
		i++;   
	}

	tmp->up = NULL;  
	sl->nElems++; 

	return 0;
}

int Delete(skip_list *sl, int elem) {
	bool present; 

	if(Exist(*sl)) {
		Search((*sl), elem, &present);
		while(present) { 
			deleteOneElement(&(*sl), elem);
			Search((*sl), elem, &present); 
		}
		return SUCCESS; 
	} else {
		return NOT_DELETED_LIST; 
	}
}

int Length(skip_list sl, int *length) {
	if(Exist(sl)) {
		*length = sl.nElems; 
		return SUCCESS; 
	} else {
		return NOT_DELETED_LIST; 
	} 
}

int Search(skip_list sl, int elem, bool *found) {
	skip_node *tmp; 
	tmp = sl.head; 
	*found = false; 
	if(Exist(sl)) {
		while(!(*found) && elem >= tmp->key) {
			if(elem == tmp->key) *found = true; 
			else {
				if(tmp->next->key > elem && tmp->level != 0) tmp = tmp->down; 	
				else tmp = tmp->next; 
				} 
			
		}
		return SUCCESS; 
	} else return NOT_DELETED_LIST; 
	
}

int Cost_Search(skip_list sl, int elem, int *cost) {
	skip_node *tmp; 
	bool found; 

	tmp = sl.head; 
	found = false; 

	*cost = 0; 
	if(Exist(sl)) {
		while(!found && elem >= tmp->key) {
			if(elem == tmp->key) found = true; 
			else {
				if(tmp->next->key > elem && tmp->level != 0) tmp = tmp->down;
				else tmp = tmp->next; 
				(*cost)++; 
			}
		}
		if(found) {
			return SUCCESS;
		} else {
			return ELEMENT_NOT_FOUND; 
		} 
	} else {
		return NOT_DELETED_LIST; 
	}
}

bool Exist(skip_list sl) {
	if(sl.nElems != -1) {
		return true; 
	} else {
		return false; 
	}
} 

int createNewUpperLayer(skip_list *sl) {
	skip_node *tmp, *tmp2; 

	tmp = sl->head; 
	tmp2 = sl->tail; 

	sl->head->up = (skip_node*)malloc(sizeof(skip_node)); 
	sl->tail->up = (skip_node*)malloc(sizeof(skip_node)); 

	if((sl->head->up != NULL) && (sl->tail->up != NULL)) {	
		sl->head->up->level = sl->head->level + 1; 
		sl->tail->up->level = sl->tail->level + 1; 

		sl->head = sl->head->up;  
		sl->tail = sl->tail->up; 

		sl->head->key = INFINITE_NEGATIVE; 
		sl->tail->key = INFINITE_POSITIVE;

		sl->head->next = sl->tail; 
		sl->tail->prev = sl->head; 

		sl->head->down = tmp; 
		sl->tail->down = tmp2; 

		sl->head->up = NULL; 
		sl->tail->up = NULL; 

		return SUCCESS; 
	} else return ERROR_CREATION_CAPA_SUPERIOR; 
	
}

int deleteUpperLayer(skip_list *sl) {  
	sl->head = sl->head->down;  
	sl->tail = sl->tail->down; 

	free(sl->head->up); 
	free(sl->tail->up); 

	sl->head->up = NULL; 
	sl->tail->up = NULL;

	return SUCCESS; 
}
 
int deleteOneElement(skip_list *sl, int elem) {
	skip_node *tmp, *tmp2; 
	bool found, end; 

	tmp = sl->head; 
	found = false; 
	while(!found && elem >= tmp->key) {
		if(elem == tmp->key) found = true; 
		else {
			if(tmp->next->key > elem && tmp->level != 0) tmp = tmp->down;
			else tmp = tmp->next; 
		}
	}
	if(found) {
		end = false; 
		do {
			if((tmp->prev->key == INFINITE_NEGATIVE) && (tmp->next->key == INFINITE_POSITIVE) && (tmp->level > 0))
				deleteUpperLayer(&(*sl)); 
			else {
				tmp->prev->next = tmp->next;
				tmp->next->prev = tmp->prev;   
			}

			tmp2 = tmp; 
			if(tmp->level == 0) end = true; 
			else tmp = tmp->down; 
			free(tmp2);   
		} while(!end && tmp != NULL); 
		sl->nElems--;
		return SUCCESS;
	} else {
		return ELEMENT_NOT_FOUND; 
	} 	
}



