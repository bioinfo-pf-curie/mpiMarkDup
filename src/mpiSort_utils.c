/*
   mpiSORT
   Copyright (C) 2016-2019 Institut Curie / Institut Pasteur

   mpiSORT is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   mpiSORT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser Public License
   along with mpiSORT.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
   Module:
     mpiSort_utils.c

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

#include "log.h"
#include "parser.h"

void get_coordinates_and_offset_source_and_size_and_free_reads(int rank, int *local_read_rank, size_t *coordinates, size_t* offset, int* size, Read* data_chr, size_t local_readNum){

   	size_t j;
   	Read* chr = data_chr;

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
        chr = chr->next;
        //free(tmp);
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
