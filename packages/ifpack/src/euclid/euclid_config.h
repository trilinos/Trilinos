/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/


#ifndef EUCLID_CONF_DH
#define EUCLID_CONF_DH

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
