#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>
#include <mpi.h>
#include "mpiSortUtils.h"

static void empty_list(void **state) {
    (void)state;
    skip_list list;
    assert_int_equal(Create(&list), SUCCESS);
    assert_int_equal(Destroy(&list), SUCCESS);
    assert_false(Exist(list));
    assert_null(list.head);
    assert_null(list.tail);
    assert_int_equal(Destroy(&list), ERROR_DESTROY);
}

static void populated_list(void **state) {
    (void)state;
    skip_list list;
    bool found;
    assert_int_equal(Create(&list), SUCCESS);
    for (int i = 0; i < 100; i++) {
        assert_int_equal(Insert(&list, i), SUCCESS);
    }
    assert_int_equal(Search(list, 42, &found), SUCCESS);
    assert_true(found);
    assert_int_equal(Destroy(&list), SUCCESS);
    assert_false(Exist(list));
}

int main(void) {
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(empty_list),
        cmocka_unit_test(populated_list),
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
