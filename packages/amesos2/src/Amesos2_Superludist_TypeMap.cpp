// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_Superludist_TypeMap.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Fri Jul 22 11:23:47 2011
 * 
 * \brief  Definitions for SuperLU_DIST TypeMap
 */

#include "Amesos2_Superludist_TypeMap.hpp"

namespace Amesos2 {
  
  SLUD::Dtype_t TypeMap<Superludist,double>::dtype = SLUD::SLU_D;

#ifdef HAVE_TEUCHOS_COMPLEX
  SLUD::Dtype_t TypeMap<Superludist,std::complex<double> >::dtype = SLUD::SLU_Z;

  SLUD::Dtype_t TypeMap<Superludist,SLUD::Z::doublecomplex>::dtype = SLUD::SLU_Z;
#endif
  
}

#ifdef HAVE_TEUCHOS_COMPLEX
namespace std {
  ostream& operator<<(ostream& out, const SLUD::Z::doublecomplex z){
    return (out << "(" << z.r << "," << z.i << ")");
  }
}
#endif
