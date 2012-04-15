/*
// @HEADER
// ***********************************************************************
// 
//    GlobiPack: Collection of Scalar 1D globalizaton utilities
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#ifndef GLOBIPACK_TEST_LAGR_POLY_MERIT_FUNC_1D_DEF_HPP
#define GLOBIPACK_TEST_LAGR_POLY_MERIT_FUNC_1D_DEF_HPP


#include "GlobiPack_TestLagrPolyMeritFunc1D_decl.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Assert.hpp"


namespace GlobiPack {


template<typename Scalar>
TestLagrPolyMeritFunc1D<Scalar>::TestLagrPolyMeritFunc1D(
  const ArrayView<const Scalar> &alpha,
  const ArrayView<const Scalar> &phi
  )
  : alpha_(alpha), phi_(phi)
{
  TEUCHOS_ASSERT_EQUALITY(alpha.size(), phi.size());
}


// Overridden from MeritFunc1DBase

  
template<typename Scalar>
bool TestLagrPolyMeritFunc1D<Scalar>::supportsDerivEvals() const
{
  return true;
}


template<typename Scalar>
void TestLagrPolyMeritFunc1D<Scalar>::eval(
  const Scalar &alpha, const Ptr<Scalar> &phi_out,
  const Ptr<Scalar> &Dphi_out
  ) const
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  const int n = alpha_.size();

  Scalar phi = ST::zero();
  Scalar Dphi = ST::zero();

  for (int k = 0; k < n; ++k) {

    if (!is_null(phi_out)) {

      Scalar Lp_k = ST::one();
      for (int i = 0; i < n; ++i) {
        if (i!=k) {
          Lp_k *= (alpha-alpha_[i])/(alpha_[k]-alpha_[i]);
        }
      }
      
      phi += phi_[k] * Lp_k;

    }

    if (!is_null(Dphi_out)) {

      Scalar DLp_k = ST::zero();
      for (int j = 0; j < n; ++j) {
        if (j!=k) {
          Scalar DLp_k_j_prod = ST::one();
          for (int i = 0; i < n; ++i) {
            if (i!=k && i!=j) {
              DLp_k_j_prod *= (alpha-alpha_[i])/(alpha_[k]-alpha_[i]);
            }
          }
          DLp_k += DLp_k_j_prod / (alpha_[k]-alpha_[j]);
        }
      }
    
      Dphi += phi_[k] * DLp_k;

    }

  }
  
  if (!is_null(phi_out)) {
    *phi_out = phi;
  }
  
  if (!is_null(Dphi_out)) {
    *Dphi_out = Dphi;
  }

}


} // namespace GlobiPack


#endif // GLOBIPACK_TEST_LAGR_POLY_MERIT_FUNC_1D_DEF_HPP
