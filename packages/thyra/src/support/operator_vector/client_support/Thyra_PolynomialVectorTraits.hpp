// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_POLYNOMIAL_VECTOR_TRAITS_HPP
#define THYRA_POLYNOMIAL_VECTOR_TRAITS_HPP

#include "Thyra_VectorBase.hpp"
#include "Teuchos_PolynomialTraits.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Teuchos {

  //! Specilization of Teuchos::PolynomialTraits for %Thyra vectors.
  /*!
   * This class provides a specilization of Teuchos::PolynomialTraits for
   * Thyra::VectorBase vectors, allowing these vectors to be coefficients
   * in the Teuchos::Polynomial.
   */
  template <typename Scalar>
  class PolynomialTraits< Thyra::VectorBase<Scalar> > {
  public:

    //! Typename of coefficients
    typedef Thyra::VectorBase<Scalar> coeff_type;

    //! Typename of scalars
    typedef Scalar scalar_type;

    //! Clone a coefficient
    static inline Teuchos::RCP<coeff_type> clone(const coeff_type& c) {
      return c.clone_v();
    }

    //! Copy a coefficient
    static inline void copy(const coeff_type& x, coeff_type* y) {
      Thyra::copy(x,y);
    }

    //! Assign a scalar to a coefficient
    static inline void assign(coeff_type* y, const scalar_type& alpha) {
      Thyra::assign(y,alpha);
    }

    //! y = x + beta*y
    static inline void update(coeff_type* y, const coeff_type& x, 
            const scalar_type& beta) {
      Thyra::Vp_V(y,x,beta);
    }

  }; // class PolynomialTraits< Thyra::VectorBase<Scalar> >

} // end namespace Teuchos

#endif  // THYRA_POLYNOMIAL_VECTOR_TRAITS_HPP
