// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
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

  SLU::Dtype_t TypeMap<Superlu,SLU::C::complex>::dtype = SLU::SLU_C;

  SLU::Dtype_t TypeMap<Superlu,SLU::Z::doublecomplex>::dtype = SLU::SLU_Z;
#endif
  
}

#ifdef HAVE_TEUCHOS_COMPLEX
namespace std {
  ostream& operator<<(ostream& out, const SLU::Z::doublecomplex z){
    return (out << "(" << z.r << "," << z.i << ")");
  }

  ostream& operator<<(ostream& out, const SLU::C::complex c){
    return (out << "(" << c.r << "," << c.i << ")");
  }
}
#endif
