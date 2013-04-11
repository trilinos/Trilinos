/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifndef __ZOLTAN_UTIL_H
#define __ZOLTAN_UTIL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* This block should be executed for Autotools and CMake builds. */
/* The Zoltan classic build defines TRILINOS_NO_CONFIG_H.        */

#ifndef TRILINOS_NO_CONFIG_H

#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

/* This file passes values from configure to the source code. */
#include "Zoltan_config.h"

#ifdef HAVE_PARMETIS
#define ZOLTAN_PARMETIS
#endif

#ifdef HAVE_METIS
#define ZOLTAN_METIS
#endif

#ifdef HAVE_SCOTCH
#define ZOLTAN_SCOTCH
#  ifdef HAVE_MPI
#  define ZOLTAN_PTSCOTCH
#  endif
#endif

#ifdef HAVE_OVIS
#define ZOLTAN_OVIS
#endif

#ifdef HAVE_PATOH
#define ZOLTAN_PATOH
#endif

#ifdef HAVE_CCOLAMD
#define ZOLTAN_CCOLAMD
#endif

#ifdef HAVE_ZOLTAN_HUND
#define ZOLTAN_HUND
#endif

#ifdef HAVE_PARKWAY
#define ZOLTAN_PARKWAY
#endif

#ifdef HAVE_PURIFY
#define ZOLTAN_PURIFY
#define strcasecmp Zoltan_strcasecmp
#define strncasecmp Zoltan_strncasecmp
#endif

#else /* !AUTOTOOLS_BUILD */

  /* With the manual build system we support only Parallel Version of Scotch */

#ifdef ZOLTAN_SCOTCH
#define ZOLTAN_PTSCOTCH
#endif

#endif /* !AUTOTOOLS_BUILD */

#define ZOLTAN_HIER
/*****************************************************************************/
/* 
 *  Macros and definitions that are common to all Zoltan modules and 
 *  utilities.
 */
/*****************************************************************************/

/*****************************************************************************/
/*
 *  Macros for consistently printing error and warning messages.
 */
/*****************************************************************************/

void Zoltan_add_back_trace(char *yo);
void Zoltan_remove_back_trace();
void Zoltan_print_trace(int rank);

#define ZOLTAN_PRINT_ERROR(proc,yo,str) { \
  fprintf(stderr, "[%d] Zoltan ERROR in %s (line %d of %s):  %s\n", \
          proc, yo, __LINE__, __FILE__, str);                       \
  Zoltan_print_trace(proc); }

#define ZOLTAN_PRINT_WARN(proc,yo,str) \
  fprintf(stderr, "[%d] Zoltan WARNING in %s (line %d of %s):  %s\n", \
          proc, yo, __LINE__, __FILE__, str);

#define ZOLTAN_TRACE(proc,where,yo,str) \
  printf("ZOLTAN (Processor %d) %s %s  %s\n", (proc), (where), (yo), \
         ((str) != NULL ? (char *)(str) : " "));

#define ZOLTAN_TRACE_IN(proc,yo,str) \
  ZOLTAN_TRACE((proc),"Entering",(yo),(str));

#define ZOLTAN_TRACE_OUT(proc,yo,str) \
  ZOLTAN_TRACE((proc),"Exiting",(yo),(str));

#define ZOLTAN_PRINT_INFO(proc,yo,str) \
  printf("ZOLTAN (Processor %d) %s: %s\n", (proc), (yo), \
         ((str) != NULL ? (char *)(str) : " "));

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* !__ZOLTAN_UTIL_H */
