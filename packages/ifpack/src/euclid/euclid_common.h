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

#ifndef COMMON_DH
#define COMMON_DH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

#define REAL_DH double

/*-----------------------------------------------------------------------
 * compile-time dependent includes from other libraries.
 * maintainer's note: this is the only place where non-Euclid
 * files are included.
 *-----------------------------------------------------------------------*/


#include <mpi.h>


/*-----------------------------------------------------------------------
 * Euclid includes
 *-----------------------------------------------------------------------*/

#include "euclid_config.h"	/* contains various user-configurable settings;
				   edit this when building an interface with
				   other libraries.  
				 */

#include "macros_dh.h"		/* macros for error checking, etc */

/*----------------------------------------------------------- 
 *  Euclid classes 
 *-----------------------------------------------------------*/
typedef struct _matgenfd *MatGenFD;
typedef struct _subdomain_dh *SubdomainGraph_dh;
typedef struct _timer_dh *Timer_dh;
typedef struct _parser_dh *Parser_dh;
typedef struct _timeLog_dh *TimeLog_dh;
typedef struct _mem_dh *Mem_dh;
typedef struct _mat_dh *Mat_dh;
typedef struct _factor_dh *Factor_dh;
typedef struct _vec_dh *Vec_dh;
typedef struct _numbering_dh *Numbering_dh;
typedef struct _hash_dh *Hash_dh;
typedef struct _hash_i_dh *Hash_i_dh;
typedef struct _mpi_interface_dh *Euclid_dh;
typedef struct _sortedList_dh *SortedList_dh;
typedef struct _extrows_dh *ExternalRows_dh;
typedef struct _stack_dh *Stack_dh;
typedef struct _queue_dh *Queue_dh;
typedef struct _sortedset_dh *SortedSet_dh;

/*
typedef struct _localPerm_dh*       LocalPerm_dh;
typedef struct _procGrid_dh*        ProcGrid_dh;
typedef struct _globalPerm_dh*      GlobalPerm_dh;
typedef struct _apply_dh*           Apply_dh;
typedef struct _externalRows_dh*    ExternalRows_dh;
*/

/*---------------------------------------------------------------------
 * misc.
 *---------------------------------------------------------------------*/


#if defined(__cplusplus)
#else
typedef int bool;
#define true   1
#define false  0
#endif

/* ------------------------------------------------------------------
 * Globally scoped variables, error handling functions, etc.
 * These are all defined in /src/globalObjects.c 
 * ------------------------------------------------------------------*/
extern Parser_dh parser_dh;	/* for setting/getting runtime options */
extern TimeLog_dh tlog_dh;	/* internal timing  functionality */
extern Mem_dh mem_dh;		/* memory management */
extern FILE *logFile;
extern int np_dh;		/* number of processors and subdomains */
extern int myid_dh;		/* rank of this processor (and subdomain) */
extern MPI_Comm comm_dh;


extern bool ignoreMe;		/* used to stop compiler complaints */
extern int ref_counter;		/* for internal use only!  Reference counter
				   to ensure that global objects are not
				   destroyed when Euclid's destructor is called,
				   and more than one instance of Euclid has been
				   instantiated.
				 */


/* Error and message handling.  These are accessed through
 * macros defined in "macros_dh.h"
 */
extern bool errFlag_dh;

#ifdef __cplusplus
extern "C"
{
#endif

  extern void setInfo_dh (char *msg, char *function, char *file, int line);
  extern void setError_dh (char *msg, char *function, char *file, int line);
  extern void printErrorMsg (FILE * fp);

#ifndef MPI_MAX_ERROR_STRING
#define MPI_MAX_ERROR_STRING 256
#endif

#define MSG_BUF_SIZE_DH MAX(1024, MPI_MAX_ERROR_STRING)
  extern char msgBuf_dh[MSG_BUF_SIZE_DH];

/* Each processor (may) open a logfile.
 * The bools are switches for controlling the amount of informational 
 * output, and where it gets written to.  Function trace logging is only 
 * enabled when compiled with the debugging (-g) option.
 */
  extern void openLogfile_dh (int argc, char *argv[]);
  extern void closeLogfile_dh ();
  extern bool logInfoToStderr;
  extern bool logInfoToFile;
  extern bool logFuncsToStderr;
  extern bool logFuncsToFile;
  extern void Error_dhStartFunc (char *function, char *file, int line);
  extern void Error_dhEndFunc (char *function);
  extern void dh_StartFunc (char *function, char *file, int line,
			    int priority);
  extern void dh_EndFunc (char *function, int priority);
  extern void printFunctionStack (FILE * fp);

  extern void EuclidInitialize (int argc, char *argv[], char *help);	/* instantiates global objects */
  extern void EuclidFinalize ();	/* deletes global objects */
  extern bool EuclidIsInitialized ();
  extern void printf_dh (char *fmt, ...);
  extern void fprintf_dh (FILE * fp, char *fmt, ...);

  /* echo command line invocation to stdout.
     The "prefix" string is for grepping; it may be NULL.
   */
  extern void echoInvocation_dh (MPI_Comm comm, char *prefix, int argc,
				 char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
