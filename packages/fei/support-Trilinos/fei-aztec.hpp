#ifndef _fei_aztec_hpp_
#define _fei_aztec_hpp_
/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef FEI_SER

#ifndef AZTEC_MPI
#define AZTEC_MPI AZTEC_MPI
#endif
#ifndef AZ_MPI
#define AZ_MPI AZ_MPI
#endif
#ifndef MPI
#define MPI MPI
#endif
#ifndef PLL
#define PLL
#endif

#endif

//az_aztec.h redefines max and min, and sometimes Warnings occur if you
//redefine macros without undef'ing them first...
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

#ifdef F77_FUNC
#undef F77_FUNC
#endif
#ifdef F77_FUNC_
#undef F77_FUNC_
#endif

#include <az_aztec.h>

#include <fei_Aztec_Map.hpp>
#include <fei_Aztec_Vector.hpp>
#include <fei_AztecDMSR_Matrix.hpp>
#include <fei_Aztec_BlockMap.hpp>
#include <fei_AztecDVBR_Matrix.hpp>

#include <fei_Aztec_LinSysCore.hpp>

#endif /* _fei_aztec_h_ */

