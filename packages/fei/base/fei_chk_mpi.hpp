
#ifndef _fei_chk_mpi_hpp_
#define _fei_chk_mpi_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_mpi.h>

#ifdef CHK_MPI
#undef CHK_MPI
#endif

static const char fei_mpi_file[] = __FILE__;

#define CHK_MPI(a) { int snl_fei_mpiErrorCode = a; \
                     if (snl_fei_mpiErrorCode != MPI_SUCCESS) {\
                      fei::console_out() << fei_mpi_file << ", line " << __LINE__  \
                           <<" MPI ERROR " << FEI_ENDL; \
                      return(snl_fei_mpiErrorCode); \
                    } }
#endif

