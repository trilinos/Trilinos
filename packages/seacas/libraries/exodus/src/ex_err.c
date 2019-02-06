/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE2 USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "exodusII.h" // for exoptval, MAX_ERR_LENGTH, etc
#include "exodusII_int.h"
#include <stdio.h>  // for fprintf, stderr, fflush
#include <stdlib.h> // for exit
#include <string.h> // for strcpy

/*!
\fn{void ex_err_fn(exoid, const char *module_name, const char *message, int err_num)}

The function ex_err_fn(exoid, ) logs an error to stderr. It is intended
to provide explanatory messages for error codes returned from other
exodus routines.

The passed in error codes and corresponding messages are listed in
???. The programmer may supplement the error message printed
for standard errors by providing an error message. If the error code
is provided with no error message, the predefined message will be
used. The error code EX_MSG is available to log application
specific messages.

\param[in]  module_name  This is a string containing the name of the calling
function.
\param[in]  message      This is a string containing a message explaining the
error
                         or problem. If EX_VERBOSE (see ex_opts()) is true,
                         this message will be printed to stderr. Otherwise,
                         nothing will be printed. Maximum length is \c
MAX_ERR_LENGTH.

\param[in] err_num       This is an integer code identifying the error. exodus C
functions
                         place an error code value in exerrval, an external
int. Negative
                         values are considered fatal errors while positive
values are
                         warnings. There is a set of predefined values defined
in
                         \file{exodusII.h}. The predefined constant \c
EX_PRTLASTMSG will
                         cause the last error message to be output, regardless
of the setting
                         of the error reporting level (see ex_opts()).

The following is an example of the use of this function:

~~~{.c}
int exoid, CPU_word_size, IO_word_size, errval;
float version;
char errmsg[MAX_ERR_LENGTH];

CPU_word_size = sizeof(float);
IO_word_size = 0;

\comment{open exodus file}
if (exoid = ex_open ("test.exo", EX_READ, &CPU_word_size,
                     &IO_word_size, &version)) {
   errval = 999;
   snprintf(errmsg, MAX_ERR_LENGTH,"ERROR: cannot open file test.exo");
   ex_err_fn(exoid, __func__, errmsg, errval);
}
~~~

*/

#if defined(EXODUS_THREADSAFE)
EX_errval_t *ex_errval = NULL;
#define EX_PNAME ex_errval->last_pname
#define EX_ERRMSG ex_errval->last_errmsg
#define EX_ERR_NUM ex_errval->last_err_num
#else
int exerrval = 0; /* clear initial global error code value */

static char last_pname[MAX_ERR_LENGTH + 1];
static char last_errmsg[MAX_ERR_LENGTH + 1];
static int  last_err_num;

#define EX_PNAME last_pname
#define EX_ERRMSG last_errmsg
#define EX_ERR_NUM last_err_num
#endif

void ex_reset_error_status()
{
#if !defined(EXODUS_THREADSAFE)
  exerrval   = 0;
  EX_ERR_NUM = 0;
#endif
}

void ex_err(const char *module_name, const char *message, int err_num)
{
  EX_FUNC_ENTER_INT();
  if (err_num == 0) { /* zero is no error, ignore and return */
    exerrval = err_num;
    EX_FUNC_VOID();
  }

  /* save the error message for replays */
  if (message != NULL) {
    strncpy(EX_ERRMSG, message, MAX_ERR_LENGTH);
    EX_ERRMSG[MAX_ERR_LENGTH - 1] = '\0';
  }
  if (module_name != NULL) {
    strncpy(EX_PNAME, module_name, MAX_ERR_LENGTH);
    EX_PNAME[MAX_ERR_LENGTH - 1] = '\0';
  }

  if (err_num == EX_PRTLASTMSG) {
    fprintf(stderr, "\n[%s] %s\n", EX_PNAME, EX_ERRMSG);
    fprintf(stderr, "    exerrval = %d\n", EX_ERR_NUM);
    if (EX_ERR_NUM < 0) {
      fprintf(stderr, "\t%s\n", ex_strerror(EX_ERR_NUM));
    }
    EX_FUNC_VOID();
  }

  if (err_num == EX_LASTERR) {
    err_num = EX_ERR_NUM;
  }
  else {
    exerrval   = err_num;
    EX_ERR_NUM = err_num;
  }

  if (err_num == EX_NULLENTITY) {
    if (exoptval & EX_NULLVERBOSE) {
      fprintf(stderr, "\nExodus Library Warning: [%s]\n\t%s\n", module_name, message);
    }
  }

  else if (exoptval & EX_VERBOSE) { /* check see if we really want to hear this */
    fprintf(stderr, "\nExodus Library Warning/Error: [%s]\n\t%s\n", module_name, message);
    if (err_num < 0) {
      fprintf(stderr, "\t%s\n", ex_strerror(err_num));
    }
  }
  fflush(stderr);

  /* with netCDF 3.4, (fatal) system error codes are > 0;
     so all EXODUS fatal error codes are > 0    */
  if ((err_num > 0) && (exoptval & EX_ABORT)) {
    exit(err_num);
  }
  EX_FUNC_VOID();
}

void ex_err_fn(int exoid, const char *module_name, const char *message, int err_num)
{
  EX_FUNC_ENTER_INT();
  if (err_num == 0) { /* zero is no error, ignore and return */
    exerrval = err_num;
    EX_FUNC_VOID();
  }

  /* save the error message for replays */
  if (message != NULL) {
    strncpy(EX_ERRMSG, message, MAX_ERR_LENGTH);
    EX_ERRMSG[MAX_ERR_LENGTH - 1] = '\0';
  }
  if (module_name != NULL) {
    strncpy(EX_PNAME, module_name, MAX_ERR_LENGTH);
    EX_PNAME[MAX_ERR_LENGTH - 1] = '\0';
  }

  if (err_num == EX_PRTLASTMSG) {
    fprintf(stderr, "\n[%s] %s\n", EX_PNAME, EX_ERRMSG);

    struct ex_file_item *file = ex_find_file_item(exoid);
    if (file) {
      size_t pathlen = 0;
      nc_inq_path(exoid, &pathlen, NULL);
      char *path = NULL;
      if (pathlen > 0) {
        path = malloc(pathlen + 1);
        if (path != NULL) {
          nc_inq_path(exoid, NULL, path);
          fprintf(stderr, "    in file '%s'", path);
          free(path);
        }
      }
    }

    fprintf(stderr, "    exerrval = %d\n", EX_ERR_NUM);

    if (EX_ERR_NUM < 0) {
      fprintf(stderr, "\t%s\n", ex_strerror(EX_ERR_NUM));
    }
    EX_FUNC_VOID();
  }

  if (err_num == EX_LASTERR) {
    err_num = EX_ERR_NUM;
  }
  else {
    exerrval   = err_num;
    EX_ERR_NUM = err_num;
  }

  if (err_num == EX_NULLENTITY) {
    if (exoptval & EX_NULLVERBOSE) {
      fprintf(stderr, "\nExodus Library Warning: [%s]\n\t%s\n", module_name, message);
    }
  }

  else if (exoptval & EX_VERBOSE) { /* check see if we really want to hear this */
    char *               path = NULL;
    struct ex_file_item *file = ex_find_file_item(exoid);
    if (file) {
      size_t pathlen = 0;
      nc_inq_path(exoid, &pathlen, NULL);
      if (pathlen > 0) {
        path = malloc(pathlen + 1);
        nc_inq_path(exoid, NULL, path);
      }
    }
    if (path) {
      fprintf(stderr, "\nExodus Library Warning/Error: [%s] in file '%s'\n\t%s\n", module_name,
              path, message);
      free(path);
    }
    else {
      fprintf(stderr, "\nExodus Library Warning/Error: [%s]\n\t%s\n", module_name, message);
    }
    if (err_num < 0) {
      fprintf(stderr, "\t%s\n", ex_strerror(err_num));
    }
  }
  fflush(stderr);

  /* with netCDF 3.4, (fatal) system error codes are > 0;
     so all EXODUS fatal error codes are > 0    */
  if ((err_num > 0) && (exoptval & EX_ABORT)) {
    exit(err_num);
  }
  EX_FUNC_VOID();
}

void ex_set_err(const char *module_name, const char *message, int err_num)
{
  EX_FUNC_ENTER_INT();
  /* save the error message for replays */
  strncpy(EX_ERRMSG, message, MAX_ERR_LENGTH);
  strncpy(EX_PNAME, module_name, MAX_ERR_LENGTH);
  EX_ERRMSG[MAX_ERR_LENGTH - 1] = '\0';
  EX_PNAME[MAX_ERR_LENGTH - 1]  = '\0';
  if (err_num != EX_LASTERR) {
    /* Use last set error number, but add new function and message */
    EX_ERR_NUM = err_num;
  }
  EX_FUNC_VOID();
}

void ex_get_err(const char **msg, const char **func, int *err_num)
{
  EX_FUNC_ENTER_INT();
  if (msg) {
    (*msg) = EX_ERRMSG;
  }

  if (func) {
    (*func) = EX_PNAME;
  }

  if (err_num) {
    (*err_num) = EX_ERR_NUM;
  }
  EX_FUNC_VOID();
}

const char *ex_strerror(int err_num)
{
  switch (err_num) {
  case EX_MEMFAIL: return "Memory allocation failure";
  case EX_BADFILEMODE: return "Bad file mode -- cannot specify both EX_READ and EX_WRITE";
  case EX_BADFILEID: return "Bad file id. Could not find exodus file associated with file id.";
  case EX_WRONGFILETYPE: return "Integer sizes must match for input and output file in ex_copy.";
  case EX_LOOKUPFAIL:
    return "Id lookup failed for specified entity type. Could not find entity with specified id.";
  case EX_BADPARAM: return "Bad parameter.";
  case EX_INTERNAL: return "Internal logic error in exodus library.";
  case EX_NOTROOTID: return "File id is not the root id; it is a subgroup id.";
  case EX_NULLENTITY: return "Null entity found.";
  case EX_DUPLICATEID: return "Duplicate entity id found.";
  default: return nc_strerror(err_num);
  }
}
