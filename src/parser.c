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
     parser.c

   Authors:
    Frederic Jarlier,   Institut Curie
    Nicolas Joly,       Institut Pasteur
    Nicolas Fedy,       Institut Curie
    Leonor Sirotti,     Institut Curie
    Thomas Magalhaes,   Institut Curie
    Paul Paganiban,     Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "parser.h"
#include "log.h"

char parse_mode;
/*
size_t hash_name(char *line, size_t max) {
    
    size_t i, j, result = 0;

    for (i = 0, j = 0; i < strlen(line); i++) {
        if ( (line[i] >= 'A' && line[i] <= 'Z') || (line[i] >= 'a' && line[i] <= 'z') ) {
            j = i;
        }
    }

    i = j + 1;

    for (j = 0; i < strlen(line) && j < max; i++) {
        if (line[i] >= '0' && line[i] <= '9') {
            char tmp;
            tmp = line[i] - '0';
            result <<= 4;
            result += tmp;
            j++;
        }
    }

    return result;
}
*/

size_t hash_name(char *begin_name, char *end_name, int max)
 {
        size_t i = 0;
        int j = 0;
        size_t result = 0;
        size_t name_length = end_name - begin_name;
        char *line = begin_name;

        for(i = 0, j = 0; i < name_length; i++) {
                if( (line[i] >= 'A' && line[i] <= 'Z') || (line[i] >= 'a' && line[i] <= 'z') )
                        j = i;
        }

        i = j + 1;
        for(j = 0; i < name_length && j < max; i++) {
                if(line[i] >= '0' && line[i] <= '9') {
                        char tmp;
                        tmp = line[i] - '0';
                        result <<= 4;
                        result += tmp;
                        j++;
                }
        }
        return result;
 }

void init_goff(MPI_File mpi_filed, unsigned int headerSize, size_t fsize, int numproc, size_t *goff) {


    char *current_line = NULL;
    MPI_Status status;
    int i = 0;
    size_t j = 0;

    size_t lsize = fsize / numproc;
    goff[0] = headerSize;

    for (i = 1; i < numproc; i++) {
        goff[i] = lsize * i + headerSize;
    }

    goff[numproc] = fsize;

    for (i = 1; i < numproc; i++) {
        current_line = (char *)calloc(1000, sizeof(char));
        MPI_File_read_at(mpi_filed, (MPI_Offset)goff[i], current_line, 1000, MPI_CHAR, &status);
        assert(strlen(current_line) != 0);
        j = 0;

        while (j < fsize && current_line[j] != '\n') {
            j++;
        }

        goff[i] += (j + 1);
        free(current_line);
    }


}


int parser_paired(char *localData, size_t size2compute, int rank, size_t start_offset, unsigned char threshold,
                int nbchrom, size_t **preadNumberByChr, Read ***preads, CHTbl *ref_htbl){

                #define _read_token_tab(q) (q); do { char *tab = strchr((q), '\t'); if (!tab) goto err_ret;  (q) = tab + 1; } while (0)
                #define _read_token_tab2(q) do { char *tab = strchr((q), '\t'); if (!tab) goto err_ret;  (q) = tab + 1; } while (0)
                #define _read_token_ret(p) (p); do { char *tab = strchr((p), '\n'); if (!tab) goto err_ret;  (p) = tab + 1; } while (0)

                char *p = localData;
                char *q;
                char *q1;
                char *q2;
                char *name;
                char *chr_txt;
                char *mchr_txt;
                char *coord_txt;
		char *qual_txt;
                char chr_name[100];
                char mchr_name[100];
                unsigned char quality;
                unsigned int i, chr, nbchr = 0, mchr;
                size_t lineSize, offset_read_in_source_file;
                size_t coord = 0;
		//size_t mate_coord = 0;
                size_t *readNumberByChr;
                size_t counter = 0;
                size_t total_computed = 0;
                Read **reads = *preads;

                offset_read_in_source_file = start_offset;
                nbchr = nbchrom;
                readNumberByChr = (size_t*)calloc(nbchr, sizeof(size_t));

                while( total_computed < size2compute ){

                        for (i = 0; i<100; i++){chr_name[i]=0;mchr_name[i]=0;}
                        q = _read_token_ret(p);

                        assert(*p != '\n');
                        total_computed += p - q;
                        lineSize = p - q;
                        //get the name
                        name = _read_token_tab(q);
                        if (parse_mode == MODE_NAME) { q1 = name; q2 = q; coord = hash_name(q1, q2, 16);}

                        //get the flag
                        _read_token_tab2(q);
                        //get the chr name)
                        chr_txt = _read_token_tab(q);
                        if ( *chr_txt == '*' ) chr = (nbchr-1);
                        else{
                               memcpy(chr_name, chr_txt , q - chr_txt - 1);
                                //look the index of the chromosom in hash table
                                chr = chtbl_lookup(ref_htbl, chr_name);
                        }
//GO TO COORD
                        //coord_txt = _read_token_tab(q);
                        if (parse_mode != MODE_NAME) {
                                coord_txt = _read_token_tab(q);
                                coord = strtoull(coord_txt, &q, 10);
                        }
                        //got to qual 
                        qual_txt = _read_token_tab(q);
                        quality = strtoull(q, &qual_txt, 10);
                        //goto cigar
                        if (parse_mode != MODE_NAME) { _read_token_tab2(q); }
                        _read_token_tab2(q);
                        //goto mate mchr
                        mchr_txt = _read_token_tab(q);

                        if( *mchr_txt == '=') mchr = chr;
                        else if( *mchr_txt == '*') mchr = (nbchr-1);
                        else {
                                memcpy(mchr_name, mchr_txt , q - mchr_txt - 1);
                                mchr = chtbl_lookup(ref_htbl, mchr_name);
                        }
                        //we take the second coordinate

                        /*      
                        fprintf(stderr, " :::::[MPISORT][PARSER] found a duplicate = READS = %zu ::: MATE= %zu \n",
                                                coord, mate_coord);
                        */

                        //first we check if reads mapped on the same chromosome
                        if ((chr < nbchr-2) && (chr == mchr)){
                                        //then we found concordant reads
                                        if(quality >= threshold){

                                                reads[chr]->next = malloc(sizeof(Read));
                                                reads[chr]->next->coord = coord;
                                                reads[chr]->next->quality = quality;
                                                reads[chr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[chr]->next->offset = lineSize;
                                                reads[chr] = reads[chr]->next;
                                                readNumberByChr[chr]++;
                                        }
                        }
                        else if ((chr < (nbchr-2)) && ( mchr < (nbchr -2)) && ( chr != mchr)){

                                        //we found discordant reads
                                        if(quality >= threshold){
						
						reads[nbchr - 2]->next = malloc(sizeof(Read));
                				reads[nbchr - 2]->next->coord = coord;
                				reads[nbchr - 2]->next->quality = quality;
                				reads[nbchr - 2]->next->offset_source_file = offset_read_in_source_file;
                				reads[nbchr - 2]->next->offset = lineSize;
                				reads[nbchr - 2] = reads[nbchr - 2]->next;
                				readNumberByChr[nbchr - 2]++;

                                                reads[chr]->next = malloc(sizeof(Read));
                                                reads[chr]->next->coord = coord;
                                                reads[chr]->next->quality = quality;
                                                reads[chr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[chr]->next->offset = lineSize;
						reads[chr] = reads[chr]->next;
                                                readNumberByChr[chr]++;
                                        }

                       }
                        else if ((chr == (nbchr-1)) && ( mchr < (nbchr -2))){

                                        //we found discordant reads with one pair unmapped
                                        if(quality >= threshold){
						/*	
                                                reads[mchr]->next = malloc(sizeof(Read));
                                                reads[mchr]->next->coord = coord;
                                                reads[mchr]->next->quality = quality;
                                                reads[mchr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[mchr]->next->offset = lineSize;
                                                reads[mchr] = reads[mchr]->next;
                                                readNumberByChr[mchr]++;
						*/

						reads[nbchr - 1]->next = malloc(sizeof(Read));
            					reads[nbchr - 1]->next->offset_source_file = offset_read_in_source_file;
            					reads[nbchr - 1]->next->offset = lineSize;
            					reads[nbchr - 1] = reads[nbchr - 1]->next;
            					readNumberByChr[nbchr - 1]++;

                                        }
                        }
                        else if ((mchr == (nbchr-1)) && ( chr < (nbchr -2))){

                                        //we found discordant reads with one pair unmapped

                                        if(quality >= threshold){
						/*
                                                reads[chr]->next = malloc(sizeof(Read));
                                                reads[chr]->next->coord = coord;
                                                reads[chr]->next->quality = quality;
                                                reads[chr]->next->offset_source_file=offset_read_in_source_file;
                                                reads[chr]->next->offset = lineSize;
                                                reads[chr] = reads[chr]->next;
                                                readNumberByChr[chr]++;
						*/

						reads[nbchr - 1]->next = malloc(sizeof(Read));
            					reads[nbchr - 1]->next->offset_source_file = offset_read_in_source_file;
            					reads[nbchr - 1]->next->offset = lineSize;
            					reads[nbchr - 1] = reads[nbchr - 1]->next;
            					readNumberByChr[nbchr - 1]++;
                                        }
                        }
                        else{
                                        //we found unmapped pairs reads
                                        //fprintf(stderr, "\n \n");
                                        //fprintf(stderr, "rank %d unmapped currentLine =  %s \n", rank, currentLine);
                                        //fprintf(stderr, "rank %d chr =  %d \n", rank, chr);
                                        //fprintf(stderr, "rank %d mchr =  %d \n", rank, mchr);
                                        //fprintf(stderr, "rank %d quality =  %s \n", rank, quality);
                                        reads[nbchr-1]->next = malloc(sizeof(Read));
                                        reads[nbchr-1]->next->offset_source_file=offset_read_in_source_file;
                                        reads[nbchr-1]->next->offset = lineSize;
                                        reads[nbchr-1] = reads[nbchr-1]->next;
                                        readNumberByChr[nbchr-1]++;
                        }
			offset_read_in_source_file += lineSize;
                        q = p;
                        counter++;
                        //we update the offset_read_in_source_file
                        /*for(i=0;i<MAX_LINE_SIZE;i++){
                                currentLine[i]=0;
                        }
                        next = tokenizer(NULL, '\n', currentLine);
                        */
                }
                fprintf(stderr, "rank %d ::: counter = %zu \n", rank, counter);
                assert(total_computed == size2compute);
                for(i=0;i<nbchr;i++){
                        preadNumberByChr[0][i] += readNumberByChr[i];
                }
                free(readNumberByChr);
                return 0;

                err_ret:
                        return -2;

}





/*
void parser_paired(char *localData, int rank, size_t start_offset, unsigned char threshold,
                   int nbchrom, size_t **preadNumberByChr, char **chrNames, Read ***preads,  CHTbl *ref_htbl) {

    char *currentCarac;
    char currentLine[MAX_LINE_SIZE];
    unsigned char quality;
    int i, chr, mchr, nbchr = 0;
    int lastChr = -1;
    int next;
    size_t lineSize, offset_read_in_source_file;
    size_t coord;
    size_t *readNumberByChr;
    size_t counter = 0;
    Read **reads = *preads;

    for (i = 0; i < MAX_LINE_SIZE; i++) {
        currentLine[i] = 0;
    }

    //we take the first line *
    //before calling parsepaired, we know that localdata is at the begining of a read
    next = tokenizer(localData, '\n', currentLine);
    offset_read_in_source_file = start_offset;

    nbchr = nbchrom;
    readNumberByChr = (size_t *)calloc(nbchr, sizeof(size_t));

    while (next) {

        lineSize = strlen(currentLine) + 1;

        //we update the offset in the
        //source file
        currentLine[lineSize - 1] = '\n';
        currentLine[lineSize] = '\0';

        //GO TO FLAG
        currentCarac = strstr(currentLine, "\t");

        *currentCarac = '\0';
        currentCarac++;

        //GO TO RNAME (Chr name)
        currentCarac = strstr(currentCarac + 1, "\t");

        if (lastChr == (nbchr - 1)) {
            chr = (nbchr - 1);

        } else {
            chr = getChr(currentCarac, chrNames, nbchr);
        }

        //GO TO COORD
        currentCarac = strstr(currentCarac + 1, "\t");


        if (parse_mode == MODE_NAME) {
            coord = strtoull(currentLine, NULL, strlen(currentLine));
            //coord = strtoull(currentLine, NULL, strlen(currentLine));
            coord = hash_name(currentLine, 16);
            
            //if (!rank) {
            //    printf("%s => %zu\n", currentLine, coord);
            //}
            

            strtoull(currentCarac, &currentCarac, 10);

        } else {
            //TAKE COORD AND GO TO MAPQ
            coord = strtoull(currentCarac, &currentCarac, 10);
        }

        //TAKE MAPQ AND GO TO CIGAR
        quality = strtoull(currentCarac, &currentCarac, 10);

        //GO TO RNEXT
        currentCarac = strstr(currentCarac + 1, "\t");

        if (currentCarac[1] == '=') {
            mchr = chr;

        } else if (currentCarac[1] == '*') {
            mchr = (nbchr - 1);

        } else {
            mchr = getChr(currentCarac, chrNames, nbchr);
        }

        //first we check if reads mapped on the same chromosome

        if ((chr < nbchr - 2) && (chr == mchr) ) {
            //then we found concordant reads
            if (quality >= threshold) {

                reads[chr]->next = malloc(sizeof(Read));
                reads[chr]->next->coord = coord;
                reads[chr]->next->quality = quality;
                reads[chr]->next->offset_source_file = offset_read_in_source_file;
                reads[chr]->next->offset = lineSize;
                reads[chr] = reads[chr]->next;
                readNumberByChr[chr]++;
            }

        } else if ((chr < (nbchr - 2)) && ( mchr < (nbchr - 2))) {

            
            if (quality >= threshold) {
                
                //we found discordant reads
                reads[nbchr - 2]->next = malloc(sizeof(Read));
                reads[nbchr - 2]->next->coord = coord;
                reads[nbchr - 2]->next->quality = quality;
                reads[nbchr - 2]->next->offset_source_file = offset_read_in_source_file;
                reads[nbchr - 2]->next->offset = lineSize;
                reads[nbchr - 2] = reads[nbchr - 2]->next;
                readNumberByChr[nbchr - 2]++;

                //we add it in reads[chr] too
                //as we want to keep it in the bam
                reads[chr]->next = malloc(sizeof(Read));
                reads[chr]->next->coord = coord;
                reads[chr]->next->quality = quality;
                reads[chr]->next->offset_source_file = offset_read_in_source_file;
                reads[chr]->next->offset = lineSize;
                reads[chr] = reads[chr]->next;
                readNumberByChr[chr]++;
      
            }

        } else if ((chr == '*') && ( mchr < (nbchr - 2))) {

            //we found discordant reads with one pair unmapped
            reads[nbchr - 1]->next = malloc(sizeof(Read));
            reads[nbchr - 1]->next->offset_source_file = offset_read_in_source_file;
            reads[nbchr - 1]->next->offset = lineSize;
            reads[nbchr - 1] = reads[nbchr - 1]->next;
            readNumberByChr[nbchr - 1]++;

        } else if ((mchr == '*') && ( chr < (nbchr - 2))) {

            //we found discordant reads with one pair unmapped
            reads[nbchr - 1]->next = malloc(sizeof(Read));
            reads[nbchr - 1]->next->offset_source_file = offset_read_in_source_file;
            reads[nbchr - 1]->next->offset = lineSize;
            reads[nbchr - 1] = reads[nbchr - 1]->next;
            readNumberByChr[nbchr - 1]++;
        }

        else {
            //we found unmapped pairs reads
            reads[nbchr - 1]->next = malloc(sizeof(Read));
            reads[nbchr - 1]->next->offset_source_file = offset_read_in_source_file;
            reads[nbchr - 1]->next->offset = lineSize;
            reads[nbchr - 1] = reads[nbchr - 1]->next;
            readNumberByChr[nbchr - 1]++;
        }



        //we update the offset_read_in_source_file
        offset_read_in_source_file += lineSize;
        //we read the next line

        for (i = 0; i < MAX_LINE_SIZE; i++) {
            currentLine[i] = 0;
        }

        next = tokenizer(NULL, '\n', currentLine);

        counter++;
    }

    //fprintf(stderr, "rank %d ::: counter = %zu \n", rank, counter);
    md_log_rank_trace(rank, "counter = %zu\n", counter);

    for (i = 0; i < nbchr; i++) {
        preadNumberByChr[0][i] += readNumberByChr[i];
    }

    free(readNumberByChr);
}

int getChr(char *str, char **chrNames, int nbchr) {
    int i = 0, found = 0, size;
    char *str1 = str, *str2;

    str2 = str1 + 1;

    for (; *str2 != '\t'; str2++);

    size = strlen(str1) - strlen(str2);

    char *tmp_chr = (char *)malloc(sizeof(char) * size + 1);
    tmp_chr[0] = 0;

    for (i = 0; i < size; i++) {
        tmp_chr[i] = str1[i + 1];
    }

    tmp_chr[size - 1] = 0;

    assert(strlen(tmp_chr) != 0);

    for (i = 0, found = 0; i < nbchr && !found; i++) {
        found = !strcmp(tmp_chr, chrNames[i]);
    }

    free(tmp_chr);
    return i - 1;
}
*/
