/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */
/*-------------------------------------------------------------------------*/
/**  @file trios_logger.h
 *
 *   @brief Method prototypes for the logger API.
 *
 *   The logger API is a simple API for logging events to a file
 *   (or to stderr or stdout)
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 1014 $.
 *   $Date: 2006-10-09 15:59:10 -0600 (Mon, 09 Oct 2006) $.
 *
 */

#ifndef _LOGGER_H_
#define _LOGGER_H_

/* removes the warning when compiled with C++ */
#ifdef __GNUC__
#pragma GCC system_header
#endif

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef FALSE
#	define FALSE (0)
#endif

#ifndef TRUE
#	define TRUE (1)
#endif

//#define DISABLE_DEBUG_LOGGING

#if defined(DISABLE_DEBUG_LOGGING)

#define logging_debug(level) (0)
#define log_debug(level, ...)

#else

/**
 * @brief Boolean function that returns TRUE if we are logging
 *        debug statements.
 */
#define logging_debug(level) \
        (((level == LOG_UNDEFINED) && (default_log_level >= LOG_DEBUG)) \
        || ((level != LOG_UNDEFINED) && (level >= LOG_DEBUG)))

/**
 * @brief Inline function that outputs a DEBUG
 * message to the log file.
 *
 * @param level   The log level to use.
 * @param args    A formatted message (like printf).
 */
#define log_debug(level, ...) if (logging_debug(level)) \
                log_output("DEBUG",__FUNCTION__,__FILE__,__LINE__, ## __VA_ARGS__)

#endif

/**
 * @brief Boolean function that returns TRUE if we are logging
 *        info statements.
 */
#define logging_info(level) \
        (((level == LOG_UNDEFINED) && (default_log_level >= LOG_INFO)) \
        || ((level != LOG_UNDEFINED) && (level >= LOG_INFO)))

/**
 * @brief Inline function that outputs an INFO
 * message to the log file.
 *
 * @param level   The log level to use.
 * @param args    A formatted message (like printf).
 */
#define log_info(level, ...) if (logging_info(level)) \
                log_output("INFO",__FUNCTION__,__FILE__,__LINE__, ## __VA_ARGS__)

/**
 * @brief Boolean function that returns TRUE if we are logging
 *        warning statements.
 */
#define logging_warn(level) \
        (((level == LOG_UNDEFINED) && (default_log_level >= LOG_WARN)) \
        || ((level != LOG_UNDEFINED) && (level >= LOG_WARN)))

/**
 * @brief Inline function that outputs a WARN
 * message to the log file.
 *
 * @param level   The log level to use.
 * @param args    A formatted message (like printf).
 */
#define log_warn(level, ...) if (logging_warn(level)) \
                log_output("WARN",__FUNCTION__,__FILE__,__LINE__, ## __VA_ARGS__)

/**
 * @brief Boolean function that returns TRUE if we are logging
 *        error statements.
 */
#define logging_error(level) \
        (((level == LOG_UNDEFINED) && (default_log_level >= LOG_ERROR)) \
        || ((level != LOG_UNDEFINED) && (level >= LOG_ERROR)))

/**
 * @brief Inline function that outputs an ERROR
 * message to the log file.
 *
 * @param level   The log level to use.
 * @param args    A formatted message (like printf).
 */
#define log_error(level, ...) if (logging_error(level)) \
                log_output("ERROR",__FUNCTION__,__FILE__,__LINE__, ## __VA_ARGS__)

/**
 * @brief Boolean function that returns TRUE if we are logging
 *        error statements.
 */
#define logging_fatal(level) \
        (((level == LOG_UNDEFINED) && (default_log_level >= LOG_FATAL)) \
        || ((level != LOG_UNDEFINED) && (level >= LOG_FATAL)))


/**
 * @brief Inline function that outputs a FATAL
 * message to the log file.
 *
 * @param level   The log level to use.
 * @param args    A formatted message (like printf).
 */
#define log_fatal(level, ...) if (logging_fatal(level)) \
                log_output("FATAL",__FUNCTION__,__FILE__,__LINE__, ## __VA_ARGS__ )

enum log_level {
        LOG_UNDEFINED = -1,
        LOG_OFF = 0,
        LOG_FATAL = 1,
        LOG_ERROR = 2,
        LOG_WARN = 3,
        LOG_INFO = 4,
        LOG_DEBUG = 5,
        LOG_ALL = 6
};
typedef enum log_level log_level;

extern log_level default_log_level;

/* the functions */

#if defined(__STDC__) || defined(__cplusplus)

extern void logger_mutex_lock(void);
extern void logger_mutex_unlock(void);
extern int logger_init(const log_level debug_level, const char *file);
extern int logger_not_initialized(void);
extern void logger_set_file(FILE *);
extern FILE *logger_get_file(void);
extern void logger_set_default_level(const log_level);
extern log_level logger_get_default_level(void);

void log_output(const char *prefix,
                const char *func_name,
                const char *file_name,
                const int line_no,
                const char *msg, ...);

extern int logger_fini(void);

#endif


#ifdef __cplusplus
}
#endif

#endif /* !LOGGER_H_ */
