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
     TestUtils.c
     
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

#include <stdio.h>

#include "mpiMarkDuplicatesUtils.h"

/**
 * @date 2018/02/19
 * @author Firmin Martin
 */

#define TEXT \
"Mademoiselle Baptistine's ambition had been to be able to\n"\
"purchase a set of drawing-room furniture in yellow Utrecht vel-\n"\
"vet, stamped with a rose pattern, and with mahogany in swan's-\n"\
"neck style, with a sofa.\n"\


#define TEXT_TAB \
    "HISEQ:001:001\t99\tchr1\t10000"


static int getLine_setup(void **state) {
    char *str = strdup(TEXT);
    *state = str;
    return 0;
}

static int getLine_teardown(void **state) {
    char *str = *state;
    free(str);
    return 0;
}

/* Test token and offset */
static void getLine_ptr(void **state) {
    char *str  = *state, *tok, *offset = str;
    assert_string_equal (str, TEXT);
    getLine(&offset, &tok);
    assert_string_equal(tok, "Mademoiselle Baptistine's ambition had been to be able to\n");
    free(tok);
    getLine(&offset, &tok);
    assert_string_equal(tok, "purchase a set of drawing-room furniture in yellow Utrecht vel-\n");
    free(tok);
    getLine(&offset, &tok);
    assert_string_equal(tok, "vet, stamped with a rose pattern, and with mahogany in swan's-\n");
    free(tok);
    getLine(&offset, &tok);
    assert_string_equal(tok, "neck style, with a sofa.\n");
    free(tok);
    getLine(&offset, &tok);
    assert_null(tok);
}

/* Test return value (length of token) */
static void getLine_size(void **state) {
    char *str  = *state, *tok, *offset = str;
    assert_int_equal(getLine(&offset, &tok), strlen("Mademoiselle Baptistine's ambition had been to be able to"));
    free(tok);
    assert_int_equal(getLine(&offset, &tok), strlen("purchase a set of drawing-room furniture in yellow Utrecht vel-"));
    free(tok);
    assert_int_equal(getLine(&offset, &tok), strlen("vet, stamped with a rose pattern, and with mahogany in swan's-"));
    free(tok);
    assert_int_equal(getLine(&offset, &tok), strlen("neck style, with a sofa."));
    free(tok);
    assert_int_equal(getLine(&offset, &tok), 0);
}

static int getTokenTab_setup(void **state) {
    char *str = strdup(TEXT_TAB);
    *state = str;
    return 0;
}

static int getTokenTab_teardown(void **state) {
    char *str = *state;
    free(str);
    return 0;
}

static void getTokenTab_ptr(void **state) {
    char *str  = *state, *tok, *offset = str;
    assert_string_equal (str, TEXT_TAB);
    getTokenTab(&offset, &tok);
    assert_string_equal(tok, "HISEQ:001:001");
    free(tok);
    getTokenTab(&offset, &tok);
    assert_string_equal(tok, "99");
    free(tok);
    getTokenTab(&offset, &tok);
    assert_string_equal(tok, "chr1");
    free(tok);
    getTokenTab(&offset, &tok);
    assert_string_equal(tok, "10000");
    free(tok);
    getTokenTab(&offset, &tok);
    assert_null(tok);
}

/* Test return value (length of token) */
static void getTokenTab_size(void **state) {
    char *str  = *state, *tok, *offset = str;
    assert_int_equal(getTokenTab(&offset, &tok), strlen("HISEQ:001:001"));
    free(tok);
    assert_int_equal(getTokenTab(&offset, &tok), strlen("99"));
    free(tok);
    assert_int_equal(getTokenTab(&offset, &tok), strlen("chr1"));
    free(tok);
    assert_int_equal(getTokenTab(&offset, &tok), strlen("10000"));
    free(tok);
    assert_int_equal(getTokenTab(&offset, &tok), 0);
}

static void readBits_correctness(void** state){
    (void)state;
    assert_int_equal(readBits(99, 5), 1);
    assert_int_equal(readBits(99, 6), 1);
    assert_int_equal(readBits(99, 7), 0);
    assert_int_equal(readBits(99, 8), 0);
}

static void fastQToPhred_correctness(void** state){
    (void)state;
    assert_int_equal(fastqToPhred('D'),35);
}


int main(void) {

    const struct CMUnitTest TestgetLine[] = {
        cmocka_unit_test_setup_teardown(getLine_ptr, getLine_setup, getLine_teardown),
        cmocka_unit_test_setup_teardown(getLine_size, getLine_setup, getLine_teardown),
    };
    int failures = cmocka_run_group_tests(TestgetLine, NULL, NULL);

    const struct CMUnitTest TestgetTokenTab[] = {
        cmocka_unit_test_setup_teardown(getTokenTab_ptr, getTokenTab_setup, getTokenTab_teardown),
        cmocka_unit_test_setup_teardown(getTokenTab_size, getTokenTab_setup, getTokenTab_teardown),
    };
    failures += cmocka_run_group_tests(TestgetTokenTab, NULL, NULL);

    const struct CMUnitTest TestreadBits[] = {
        cmocka_unit_test_setup_teardown(readBits_correctness, NULL, NULL),
    };

    failures += cmocka_run_group_tests(TestreadBits, NULL, NULL);

    const struct CMUnitTest TestFastQToPhred[] = {
        cmocka_unit_test_setup_teardown(fastQToPhred_correctness, NULL, NULL),
    };
    failures += cmocka_run_group_tests(TestFastQToPhred, NULL, NULL);
    return failures ? 1 : 0;
}
