# Unit tests

The runnable tests require an MPI C compiler, CMocka, OpenSSL, and make.
For example:

```sh
conda activate mpimarkdup-mpi
make -C unit-tests check
```

`make` builds the tests; `make check` builds and runs them, returning a
nonzero status on failure. `make clean` removes generated test binaries and
coverage files. Override `CC`, `CPPFLAGS`, `CFLAGS`, `LDFLAGS`, or `LDLIBS`
as needed for a local installation.

The current suite covers utilities, perfect hashing, and the linked list.
`TestList.c`, `TestLBList.c`, and `TestClusterFunctions.c` are retained as
legacy tests but are not built: they depend on the removed `readList` API
and need to be rewritten for the current implementation. This suite does
not test MPI communication or end-to-end duplicate marking.
