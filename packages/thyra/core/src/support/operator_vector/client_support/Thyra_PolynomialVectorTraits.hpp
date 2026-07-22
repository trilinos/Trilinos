// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      Thyra::copy(x, Teuchos::ptr(y));
    }

    //! Assign a scalar to a coefficient
    static inline void assign(coeff_type* y, const scalar_type& alpha) {
      Thyra::assign(Teuchos::ptr(y), alpha);
    }

    //! y = x + beta*y
    static inline void update(coeff_type* y, const coeff_type& x, 
            const scalar_type& beta) {
      Thyra::Vp_V(Teuchos::ptr(y), x, beta);
    }

  }; // class PolynomialTraits< Thyra::VectorBase<Scalar> >

} // end namespace Teuchos

#endif  // THYRA_POLYNOMIAL_VECTOR_TRAITS_HPP
