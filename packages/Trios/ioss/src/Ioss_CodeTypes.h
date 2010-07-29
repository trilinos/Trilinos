/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_code_types_h
#define IOSS_code_types_h

#if !defined(HAVE_MPI)
#if defined(SIERRA_PARALLEL_MPI) || defined(STK_BUILT_IN_SIERRA)
#define HAVE_MPI
#else
#include <Trios_config.h>
#endif
#endif

#if defined(HAVE_MPI)
#include <mpi.h>
#else
#ifndef MPI_COMM_WORLD
#define MPI_COMM_WORLD 0
typedef int MPI_Comm;
#endif
#endif

#include <complex>
typedef std::complex<double> Complex;
#endif 
