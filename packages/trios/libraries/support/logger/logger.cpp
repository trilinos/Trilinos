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
/**  @file logger.c
 *
 *   @brief This file contains method defintions for the logger API.
 *
 *   The logger API is a simple API for logging events to a file
 *   (or to stderr or stdio)
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 1576 $.
 *   $Date: 2007-10-01 17:31:49 -0600 (Mon, 01 Oct 2007) $.
 *
 */

#include "Trios_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <iostream>
#include <sstream>

#ifdef HAVE_TRIOS_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_TRIOS_SYSCALL_H
#include <syscall.h>
#endif

#include "Trios_logger.h"
#include "Trios_threads.h"

static bool mutex_initialized = false;
static nthread_lock_t logger_mutex;

void logger_mutex_lock()
{
    if (!mutex_initialized) {
        if (nthread_lock_init(&logger_mutex) == -1) {
            fprintf(stderr, "nthread_lock_init failed.\n");
            fflush(stderr);
            return;
        }

        mutex_initialized = true;
    }

    nthread_lock(&logger_mutex);
}
void logger_mutex_unlock()
{
    if (!mutex_initialized) {
        fprintf(stderr, "logger_mutex_unlock: mutex not intialized.\n");
        fflush(stderr);
        return;
    }

    nthread_unlock(&logger_mutex);
}

log_level default_log_level = LOG_WARN;
static FILE *log_file = NULL;



/**
 * @brief Initialize the logging capabilities.
 *
 * @param level   Log level (LOG_ALL, LOG_DEBUG, LOG_INFO, LOG_WARN, LOG_ERROR, LOG_FATAL, or LOG_OFF)
 * @param logfile File name of the logging output.
 *
 */
int logger_init(const log_level debug_level,  const char *logfile)
{
    int rc = 0;

    if (!mutex_initialized) {
        if (nthread_lock_init(&logger_mutex) == -1) {
            fprintf(stderr, "nthread_lock_init failed.\n");
            fflush(stderr);
            return(-1);
        }

        mutex_initialized = true;
    }

    /* initialize the default debug level */
    if (debug_level == 0)
        logger_set_default_level(LOG_OFF);
    else if (debug_level > 5)
        logger_set_default_level(LOG_ALL);
    else {
        int new_level = (int)debug_level - LOG_OFF;
        logger_set_default_level((log_level)new_level);
    }

    /* initialize the logfile */
    if ((logfile == NULL) || (logfile[0] == '\0')) {
        logger_set_file(stdout);
    }
    else if (strcasecmp("stdout", logfile) == 0) {
        logger_set_file(stdout);
    }

    else if (strcasecmp("stderr", logfile) == 0) {
        logger_set_file(stderr);
    }

    else {
        FILE *fp = fopen(logfile, "w+");
        if (fp == NULL) {
            fprintf(stderr, "could not create log file \"%s\"\n",logfile);
            return(-1);
        }
        else {
            logger_set_file(fp);
        }
    }

    return(rc);
}

int logger_not_initialized()
{
    return log_file == NULL;
}

/**
 * @brief Set the file for the log information.
 */
void logger_set_file(FILE *newfile)
{
    log_file = newfile;
}

FILE *logger_get_file()
{
    if (!log_file)
        return stdout;
    else
        return log_file;
}


/**
 * @brief Set the default log level.
 *
 * The different log levels are LOG_ALL, LOG_DEBUG, LOG_INFO, LOG_WARN,
 * LOG_ERROR, LOG_FATAL, and LOG_OFF.
 */
void logger_set_default_level(const log_level newlevel)
{
    default_log_level = newlevel;
}

/**
 * @brief Return the default log level.
 */
log_level logger_get_default_level(void)
{
    return default_log_level;
}

/**
 * @brief Output a log message.
 *
 * This method should be called by one of the inline
 * methods log_debug, log_info, log_warn, log_error, or
 * log_fatal.
 */
void log_output(const char *prefix,
        const char *func_name,
        const char *file_name,
        const int line_num,
        const char *msg,
        ...)
{
    va_list ap;
    const char *file;
    char buf1[256];
    char buf2[256];

    if (logger_not_initialized()) {
        logger_init(LOG_ERROR, NULL);
    }

    /* path from last '/' */
    file = strrchr(file_name, '/');

    va_start(ap, msg);

#ifdef HAVE_TRIOS_GETTID
    sprintf(buf1, "%s [%s:%s:%d:t%lu]: ",
            prefix,
            func_name,
            (file == NULL) ? file_name : &(file[1]),
            line_num,
            syscall(SYS_gettid));
#else
    sprintf(buf1, "%s [%s:%s:%d]: ",
            prefix,
            func_name,
            (file == NULL) ? file_name : &(file[1]),
            line_num);
#endif

    vsprintf(buf2, msg, ap);
    logger_mutex_lock();
    fprintf(log_file, "%s %s\n", buf1, buf2);
    logger_mutex_unlock();
    va_end(ap);
    fflush(log_file);
}

/**
 * @brief Finalize the logging capabilities.
 *
 */
int logger_fini(void)
{
    int rc = 0;

    if (log_file) {
        fclose(log_file);
    }

    if (mutex_initialized) {
        nthread_lock_fini(&logger_mutex);
        mutex_initialized = false;
    }

    return rc;
}
