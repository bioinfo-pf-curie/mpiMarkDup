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
     TestList.c
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmin Martin,      Paris Descartes   
*/

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>
#include <cmocka.h>

#include "perfectHash.h"
#include "reads.h"

/**
 * @date 2018 Apr 1
 */

//note: must be even for commodity. (each read has its mate in the table).
#define TABLE_SIZE 300

typedef struct {
    hashTable *htbl;
    readInfo **reads;
} test_t;

static readInfo *readsInit(int i, int pairNumber) {
    readInfo* r = calloc(1, sizeof(readInfo));
    char *str = malloc(sizeof(char) * 100);
    sprintf(str, "QName-%04d", i);
    r->cigar = NULL;
    r->Qname = str;
    r->valueFlag = 64 * pairNumber;
    r->pair_num = pairNumber;
    r->fingerprint = read2Fingerprint(r);
    return r;
}

static int PH_setup(void **state) {
    test_t *test = malloc(sizeof(test_t));
    test->htbl = malloc(sizeof(hashTable));
    readInfo **reads = malloc(sizeof(readInfo *) * TABLE_SIZE);

    for (int i = 0; i < TABLE_SIZE; ++i) {
        // reads Qname : QName-0000, QName-0002, QName-0004 ...
        reads[i] = readsInit(i - i % 2, i % 2 + 1);
    }
    test->reads = reads;
    hashTableInit(test->htbl, reads, TABLE_SIZE);

    *state = test;
    return 0;
}

static int PH_teardown(void **state) {
    test_t *test = *state;
    hashTable *htbl = test->htbl;
    readInfo** reads = test->reads;

    hashTableDestroy(htbl);
    for (int i = 0; i < TABLE_SIZE; ++i) {
        freeRead(reads[i]);   
    }
    free(reads);
    free(test);
    return 0;
}

static void PH_getElement(void **state) {
    test_t *test = *state;
    readInfo** reads = test->reads;
    hashTable* htbl = test->htbl;

    unsigned long long fingerprint;
    readInfo* r;
    char* str = malloc(100 * sizeof(char));

    
    /** Test getReadFromFingerprint **/

    // test the first read in array
    fingerprint = reads[0]->fingerprint;
    r = getReadFromFingerprint(htbl, fingerprint);
    assert_ptr_equal(reads[0], r);

    // test the last read in array
    fingerprint = reads[TABLE_SIZE - 1]->fingerprint;
    r = getReadFromFingerprint(htbl, fingerprint);
    assert_ptr_equal(reads[TABLE_SIZE - 1], r);

    // test non-existent read
    r = getReadFromFingerprint(htbl, 0);
    assert_null(r);

    /** Test getMateFromRead **/

    // reads[0] mate is reads[1]
    r = getMateFromRead(htbl, reads[0]);
    assert_ptr_equal(reads[1], r);
    
    // reads[1] mate is reads[0]
    r = getMateFromRead(htbl, reads[1]);
    assert_ptr_equal(reads[0], r);
     
    // reads[TABLE_SIZE - 1] mate is reads[TABLE_SIZE - 2]
    r = getMateFromRead(htbl, reads[TABLE_SIZE - 1]);
    assert_ptr_equal(reads[TABLE_SIZE - 2], r);

    // getMateFromRead is idempotent...
    r = getMateFromRead(htbl, getMateFromRead(htbl, reads[0]));
    assert_ptr_equal(reads[0], r);

    /** Test getReadFromQnameAndPairNum **/

    // reads[0] Qname is QName-0000, its pair Number is 1
    r = getReadFromQnameAndPairNum(htbl, "QName-0000", 1);
    assert_ptr_equal(reads[0], r);
    
    // reads[TABLE_SIZE - 1] 
    sprintf(str, "QName-%04d", (TABLE_SIZE - 1) - (TABLE_SIZE - 1)%2);
    r = getReadFromQnameAndPairNum(htbl, str, (TABLE_SIZE - 1) % 2 + 1);
    assert_ptr_equal(reads[TABLE_SIZE - 1], r);

    // test non-existent read
    r = getReadFromQnameAndPairNum(htbl, "QName-48932793487", 1);
    assert_null(r);

    free(str);
}

static void PH_conversion(void **state) {
    test_t *test = *state;
    readInfo** reads = test->reads;
    hashTable* htbl = test->htbl;

    unsigned long long fingerprint;
    readInfo* r;

    /** Test read2mateFP **/

    r = getMateFromRead(htbl, reads[0]);
    fingerprint = read2mateFP(reads[0]); 
    assert_int_equal(r->fingerprint, fingerprint);
        
}

static void PH_insert(void **state) {
    test_t *test = *state;
    readInfo** reads = test->reads;
    hashTable* htbl = test->htbl;

    readInfo* r1, *r2;

    /** Test hashTableInsert **/
    // insert in the same place, so free first
    freeRead(reads[0]);
    r1 = readsInit(0, 1);
    reads[0] = r1;
    hashTableInsert(htbl, r1);
    r2 = getReadFromQnameAndPairNum(htbl, "QName-0000", 1);
    assert_ptr_equal(r1, r2);
        
}

static void PH_known_digests(void **state) {
    (void)state;
    static const struct {
        const char *input;
        unsigned char digest[16];
    } vectors[] = {
        {"", {0xd4, 0x1d, 0x8c, 0xd9, 0x8f, 0x00, 0xb2, 0x04,
              0xe9, 0x80, 0x09, 0x98, 0xec, 0xf8, 0x42, 0x7e}},
        {"abc", {0x90, 0x01, 0x50, 0x98, 0x3c, 0xd2, 0x4f, 0xb0,
                 0xd6, 0x96, 0x3f, 0x7d, 0x28, 0xe1, 0x7f, 0x72}},
    };

    for (size_t i = 0; i < sizeof(vectors) / sizeof(vectors[0]); i++) {
        size_t a, b;
        memcpy(&a, vectors[i].digest, sizeof(a));
        memcpy(&b, vectors[i].digest + 8, sizeof(b));
        assert_int_equal(string2MD5(vectors[i].input), a ^ b);
    }
}

int main(void) {

    const struct CMUnitTest PHPtr[] = {
        cmocka_unit_test_setup_teardown(PH_getElement, PH_setup, PH_teardown),
        cmocka_unit_test_setup_teardown(PH_conversion, PH_setup, PH_teardown),
        cmocka_unit_test_setup_teardown(PH_insert, PH_setup, PH_teardown),
        cmocka_unit_test(PH_known_digests),
    };
    return cmocka_run_group_tests(PHPtr, NULL, NULL);
}
