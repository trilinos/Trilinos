/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef SIERRA_shared_fortran_h
#define SIERRA_shared_fortran_h

#include <SIERRA_code_types.h>
#include <stk_util/util/Fortran.hpp>


extern "C" void
SIERRA_FORTRAN(asol_ug3d_so_ratedef)( const Int& nelem,
			       const Int& nint,
			       const Real& dt,
			       const Real* def_grad_inv,
				     Real* lcg_inv,
				     Real* eigval,
				     Real* eigvec,
				     Real* rate_def );

#endif
