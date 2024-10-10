/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/


#ifndef EUCLID_CONF_DH
#define EUCLID_CONF_DH

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#define MAX_MPI_TASKS 50000

/* for use in printTriples functions */
#define TRIPLES_FORMAT    "%i %i %1.8e\n"
/* #define TRIPLES_FORMAT    "%i %i %1.19e\n" */

#undef PRIVATE_TIMING_DH
  /* primarily for experimental purposes; if defined, TimeLog_dh
     marks are entered in Mat_dh and Euclid_apply modules.
   */


  /* top-level error handlers. redefine to do what you want, or don't
     use it at all.  Intended usage for calling Euclid functions from
     main is:

     Euclid_dhPhoo(); ERRCHKA;
   */

#ifdef USING_MPI
#define EUCLID_EXIT MPI_Abort(comm_dh, -1)
#else
#define EUCLID_EXIT exit(-1);
#endif

#define EXIT_NOW(msg) \
      { setError_dh(msg, __FUNC__, __FILE__, __LINE__); \
        ERRCHKA; \
      }

#define ERRCHKA   \
    if (errFlag_dh) {  \
      setError_dh("", __FUNC__, __FILE__, __LINE__); \
      if (logFile != NULL) {  \
        printErrorMsg(logFile);  \
        closeLogfile_dh();  \
      } \
      printErrorMsg(stderr);  \
      if (myid_dh == 0) { \
        Mem_dhPrint(mem_dh, stderr, false); \
      } \
      EUCLID_EXIT; \
    }

#define ERRCHKA_CHKERRA(ierr)   \
    if (errFlag_dh) {  \
      setError_dh("", __FUNC__, __FILE__, __LINE__); \
      if (logFile != NULL) {  \
        printErrorMsg(logFile);  \
        fprintf(logFile, "\n[%i] ierr = %i, errFlag_dh = %i\n", myid_dh, ierr, errFlag_dh); \
        closeLogfile_dh();  \
      } \
      printErrorMsg(stderr);  \
      fprintf(stderr, "\n[%i] ierr = %i, errFlag_dh = %i\n", myid_dh, ierr, errFlag_dh); \
      CHKERRA(ierr); \
    }


#define MAX_SUBDOMAINS  20
  /* The maximum number of subdomains into which
     the matrix may be partitioned.  Rule of thumb:
     MAX_SUBDOMAINS >= number of threads.

     Note: this is only for shared-memory.
   */


#define PIVOT_FIX_DEFAULT  1e-3

/*---------------------------------------------------------------------
 * Memory management.  These macros work with functions in Mem_dh.c;
 * Change if you want to use some memory management and reporting schemes 
 * other than that supplied with Euclid.   These depend on the global
 * object "Mem_dh mem_dh" which is defined in globalObjects.c (yuck!)
 ---------------------------------------------------------------------*/

#define MALLOC_DH(s)  Mem_dhMalloc(mem_dh, (s))
#define FREE_DH(p)    Mem_dhFree(mem_dh, p)


  /* The actual calls used by Mem_dh objects to allocate/free memory 
   * from the heap.
   */
#define PRIVATE_MALLOC  malloc
#define PRIVATE_FREE    free

/*------------------ Memory management end -----------------------------*/

/*

This is currently accomplished in the makefile system;
If you're building an interface to a solver package,
you need to write EUCLID_GET_ROW() functions: see src/getRow.c
*/

#endif
