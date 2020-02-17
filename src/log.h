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
     log.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#ifndef MD_LOG_H
#define MD_LOG_H

/**
 * @file log.h
 */

#include <mpi.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#define LOG_USE_COLOR

#define LOG_BUF_SIZ 4096

#define MPI_ASSERT(fx) do { \
  assert(fx == MPI_SUCCESS); \
} while(0)

/// Log levels.
enum mdLogLevel {
    MD_LOG_OFF = -1,  ///< All logging disabled.
    MD_LOG_ERROR,     ///< Logging of errors only.
    MD_LOG_WARNING,   ///< Logging of errors and warnings.
    MD_LOG_INFO,      ///< Logging of errors, warnings, and normal but significant events.
    MD_LOG_DEBUG,     ///< Logging of all except the most detailed debug events.
    MD_LOG_TRACE      ///< All logging enabled.
};

/// Sets the selected log level.
void md_set_log_level(enum mdLogLevel level);
void md_set_log_comm(MPI_Comm comm) ;

/// Gets the selected log level.
int md_get_log_level();
MPI_Comm md_get_log_comm() ;

/// Selected log level.
/**
 * One of the MD_LOG_* values. The default is MD_LOG_WARNING.
 * @note Avoid direct use of this variable. Use md_set_log_level and md_get_log_level instead.
 */
extern int md_verbose;

void md_log_rank(int rank, enum mdLogLevel severity, const char *file,  const char *context, const int line,  const char *format, ...) ;
void md_log_all(enum mdLogLevel severity, const char *file,  const char *context, const int line, const char *format, ...) ;

#define md_log_error(format, ...)   \
 if ( MD_LOG_ERROR  <= md_get_log_level()) {   \
    md_log_rank(0, MD_LOG_ERROR, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_warning(format, ...) \
 if ( MD_LOG_WARNING  <= md_get_log_level()) { \
    md_log_rank(0, MD_LOG_WARNING, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_info(format, ...)    \
 if ( MD_LOG_INFO  <= md_get_log_level()) {   \
    md_log_rank(0, MD_LOG_INFO, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_debug(format, ...)   \
 if ( MD_LOG_DEBUG <= md_get_log_level())  {  \
    md_log_rank(0, MD_LOG_DEBUG, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_trace(format, ...)   \
 if ( MD_LOG_TRACE  <= md_get_log_level()) {  \
    md_log_rank(0, MD_LOG_TRACE, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }

#define md_log_rank_error(rank, format, ...)   md_log_rank(rank, MD_LOG_ERROR, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_warning(rank, format, ...) md_log_rank(rank, MD_LOG_WARNING, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_info(rank, format, ...)    md_log_rank(rank, MD_LOG_INFO, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_debug(rank, format, ...)   md_log_rank(rank, MD_LOG_DEBUG, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_trace(rank, format, ...)   md_log_rank(rank, MD_LOG_TRACE, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)

#define md_log_all_error(format, ...)          md_log_all(MD_LOG_ERROR, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_warning(format, ...)        md_log_all(MD_LOG_WARNING, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_info(format, ...)           md_log_all(MD_LOG_INFO, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_debug(format, ...)          md_log_all(MD_LOG_DEBUG, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_trace(format, ...)          md_log_all(MD_LOG_TRACE, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)

#endif /* ifndef MD_LOG_H */
