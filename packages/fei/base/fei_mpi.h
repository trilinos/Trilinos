#ifndef _fei_mpi_h_
#define _fei_mpi_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"

#ifdef FEI_SER
/**
  If FEI_SER is defined, the user wants to build/run in purely serial mode,
  without linking against MPI.
  To minimize #ifdefs in FEI code, we do a few #defines for
  some common MPI symbols that appear in the code.
*/
#define MPI_Comm int
#define MPI_Request int
#define MPI_COMM_WORLD 0
#define MPI_Abort(a, b) abort()
#define MPI_Wtime() 0.0
#define MPI_Barrier( a ) (void)a
#define MPI_SUCCESS 0

#else
#include <mpi.h>
#endif

#endif // _fei_mpi_h_

