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
     mpiMD_utils.c
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#include "mpiMarkDuplicatesUtils.h"


/**
 * @date 2018 Feb 26
 * @brief convert a digit (0-9) to char
 * @param[in] n {0, ..., 9}
 * @return the corresponding char
 */

inline char digit2char(int n) {
    char c = n + '0';
    return c;
}


/*
 * This fonction return the Phred value of a char
 */

inline int fastqToPhred (char ch) {
    int ich = ch;
    assert(ich > 33 || ich < 126);
    return (ich - 33);
}


void printRead(readInfo *read) {
    printf("Qname=%s, fingerprint=%zu, mateRank=%d, matePos=%zu, coord=%zu, check_with_bruck=%u\n", read->Qname, read->fingerprint, read->mate_rank, read->coordMatePos, read->coordPos, read->check_with_bruck);
}


/**
 *   @date 2018 Feb 26
 *   @brief allocate a line of sam and step forward offset
 *   @param[in, out] offset offset in text
 *   @param[in] text text we parse
 *   @param[out] tokenLine string allocated
 *   @return -1 if failed, 0 otherwise
 *   @todo merge getLine and getTokenTab
 */

int getLine(char **offset, char **tokenLine) {

    //TODO:Do this without strlen(text)
    //if (*offset < text || *offset > text + strlen(text)) {
    //    *offset = NULL;
    //    return -1;
    //}

    char *lineHead = *offset, *next;
    next = strchr(*offset, '\n');

    if (next) {
        *offset = next;
        int counter = (*offset - lineHead) / sizeof(char);
        (*offset)++;
        *tokenLine = strndup(lineHead, counter + 1);
        return counter;

    } else {
        *tokenLine = NULL;
        return 0;
    }
}

/**
 *   @date 2018 Feb 26
 *   @brief allocate a token delimited by tab and step forward offset
 *   @param[in, out] offset offset in text
 *   @param[in] text text we parse
 *   @param[out] tokenLine token allocated
 *   @return -1 if failed, remains text size otherwise
 */

int getTokenTab(char **offset, char **tokenTab) {
    //if (*offset < text || *offset > text + strlen(text)) {
    //    *offset = NULL;
    //    return -1;
    //}

    char *lineHead = *offset, *next;
    next = strchr(*offset, '\t');

    if (next) {
        *offset = next;
        int counter = (*offset - lineHead) / sizeof(char);
        (*offset)++;
        *tokenTab = strndup(lineHead, counter);
        return counter;

    } else {
        char *lineEnd = strchr(*offset, '\n');
        size_t remainsTextSize = lineEnd ? (lineEnd - *offset) / sizeof(char) : strlen(*offset);
        *offset += remainsTextSize;
        *tokenTab = remainsTextSize ? strndup(lineHead, remainsTextSize) : NULL;
        return remainsTextSize;
    }
}


/*
 * Function that reads n bits of x at position p
 */

inline unsigned int readBits(unsigned int x, int k) {
    // from K&R page 49
    return ((x & ( 1 << k )) >> k);
}

inline void writeBits(unsigned int *x, int k) {
    *x |= (1 << k);
}

inline void toggleBits(unsigned int *x, int k) {
    *x ^= readBits(*x, k);
}

inline void switchBits(unsigned int *x, int k, int j){
    toggleBits(x, k);
    toggleBits(x, j);
}



/**
 *   @date 2018 Mar 30
 *   @brief string appending, return concatenated string end
 *   @param[in] dest destination pointer
 *   @param[in] src source pointer
 *   @return end of the concatenated string
 */

inline char *strapp(char *dest, char *src ) {
    while (*dest) {
        dest++;
    }

    while ((*dest++ = *src++));

    return --dest;
}

