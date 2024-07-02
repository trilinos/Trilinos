// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_Superlu_TypeMap.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Fri Jul 22 11:25:30 2011
 * 
 * \brief  Definitions for SuperLU TypeMap.
 */

#include "Amesos2_Superlu_TypeMap.hpp"

namespace Amesos2 {
  
  SLU::Dtype_t TypeMap<Superlu,float>::dtype = SLU::SLU_S;

  SLU::Dtype_t TypeMap<Superlu,double>::dtype = SLU::SLU_D;

#ifdef HAVE_TEUCHOS_COMPLEX
  SLU::Dtype_t TypeMap<Superlu,std::complex<float> >::dtype = SLU::SLU_C;

  SLU::Dtype_t TypeMap<Superlu,std::complex<double> >::dtype = SLU::SLU_Z;

  SLU::Dtype_t TypeMap<Superlu,Kokkos::complex<float>>::dtype = SLU::SLU_C;

  SLU::Dtype_t TypeMap<Superlu,Kokkos::complex<double>>::dtype = SLU::SLU_Z;
#endif
  
}
