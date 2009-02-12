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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
