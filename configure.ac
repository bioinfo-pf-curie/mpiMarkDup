# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([mpimd], [1.0])
#AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([1.13 foreign -Wall])

# Checks for programs.
AC_PROG_CC([mpicc])

# Check for functions
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([memmove]) 
AC_CHECK_FUNCS([memset]) 
AC_CHECK_FUNCS([munmap]) 
AC_CHECK_FUNCS([sqrt]) 
AC_CHECK_FUNCS([strchr]) 
AC_CHECK_FUNCS([strdup]) 
AC_CHECK_FUNCS([strndup]) 
AC_CHECK_FUNCS([strstr]) 
AC_CHECK_FUNCS([strtoull])

AC_CHECK_HEADERS([fcntl.h]) 
AC_CHECK_HEADERS([limits.h]) 
AC_CHECK_HEADERS([malloc.h]) 
AC_CHECK_HEADERS([stddef.h]) 
AC_CHECK_HEADERS([stdint.h]) 
AC_CHECK_HEADERS([stdlib.h]) 
AC_CHECK_HEADERS([string.h]) 
AC_CHECK_HEADERS([unistd.h]) 
AC_CHECK_HEADER_STDBOOL

AC_C_INLINE
AC_FUNC_FSEEKO 
AC_FUNC_MALLOC 
AC_FUNC_MMAP 
AC_FUNC_REALLOC 
AC_TYPE_INT16_T 
AC_TYPE_INT64_T 
AC_TYPE_INT8_T 
AC_TYPE_OFF_T 
AC_TYPE_SIZE_T 
AC_TYPE_UINT16_T 
AC_TYPE_UINT32_T 
AC_TYPE_UINT8_T 


# Checks for libraries.
AC_CHECK_LIB([m], [log2])
AC_CHECK_LIB([m], [max])
AC_CHECK_LIB([z], [inflate])
AC_CHECK_LIB([crypto], [MD5_Init])
AC_CHECK_LIB([crypto], [MD5_Update])
AC_CHECK_LIB([crypto], [MD5_Final])

# Checks for library functions.
AC_CHECK_FUNCS([malloc_trim])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
