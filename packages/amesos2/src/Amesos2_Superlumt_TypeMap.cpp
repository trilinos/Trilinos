// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_Superlumt_TypeMap.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Fri Jul 22 11:24:39 2011
 * 
 * \brief  Definitions for SuperLU_MT TypeMap
 */


#include "Amesos2_Superlumt_TypeMap.hpp"

namespace Amesos2 {
  
  SLUMT::Dtype_t TypeMap<Superlumt,float>::dtype = SLUMT::SLU_S;

  SLUMT::Dtype_t TypeMap<Superlumt,double>::dtype = SLUMT::SLU_D;

#ifdef HAVE_TEUCHOS_COMPLEX
  SLUMT::Dtype_t TypeMap<Superlumt,std::complex<float> >::dtype = SLUMT::SLU_C;

  SLUMT::Dtype_t TypeMap<Superlumt,std::complex<double> >::dtype = SLUMT::SLU_Z;

  SLUMT::Dtype_t TypeMap<Superlumt,SLUMT::C::complex>::dtype = SLUMT::SLU_C;

  SLUMT::Dtype_t TypeMap<Superlumt,SLUMT::Z::doublecomplex>::dtype = SLUMT::SLU_Z;
#endif
  
}

#ifdef HAVE_TEUCHOS_COMPLEX
namespace std {
  ostream& operator<<(ostream& out, const SLUMT::Z::doublecomplex z){
    return (out << "(" << z.r << "," << z.i << ")");
  }

  ostream& operator<<(ostream& out, const SLUMT::C::complex c){
    return (out << "(" << c.r << "," << c.i << ")");
  }
}
#endif
