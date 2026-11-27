CC = mpicc
SRC = ../src
CPPFLAGS += -I$(SRC) -D_GNU_SOURCE
CFLAGS ?= -g -O2 -Wall -Wextra -std=gnu99
LDLIBS += -lcmocka -lcrypto -lm

BIN = TestUtils TestPerfectHash TestLlist TestSkipList
COMMON = $(SRC)/log.c $(SRC)/mpiMarkDuplicatesUtils.c

.PHONY: all check clean
all: $(BIN)

check: all
	@set -e; for test in $(BIN); do ./$$test; done

TestUtils: TestUtils.c $(COMMON)
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

TestPerfectHash: TestPerfectHash.c $(SRC)/perfectHash.c $(SRC)/reads.c $(COMMON)
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

TestLlist: TestLlist.c $(SRC)/llist.c $(SRC)/reads.c $(COMMON)
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

TestSkipList: TestSkipList.c $(SRC)/mpiSortUtils.c $(COMMON)
	$(CC) $(CPPFLAGS) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f $(BIN) *.gcda *.gcno
