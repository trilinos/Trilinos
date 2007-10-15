#ifndef _cfei_dash_aztec_h_
#define _cfei_dash_aztec_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#ifndef FEI_SER

#define AZ_MPI AZ_MPI
#ifndef MPI
#define MPI MPI
#endif
#ifndef PLL
#define PLL PLL
#endif

#include <mpi.h>

#endif

#include <cfei_aztec.h>

#endif /*_cfei_dash_aztec_h_*/

