// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_POLYNOMIAL_TRAITS_HPP
#define TEUCHOS_POLYNOMIAL_TRAITS_HPP

#include "Teuchos_RCP.hpp"

namespace Teuchos {

  //! Traits class for polynomial coefficients in Teuchos::Polynomial.
  /*!
   * This class provides traits for implementing Teuchos::Polynomial.  The
   * default template definition here will work for any scalar type.  Any other
   * coefficient type for Teuchos::Polynomial should provide a specialization
   * of this traits class for that type that mirrors the default definition
   * below.
   */
  template <typename Scalar>
  class PolynomialTraits {
  public:

    //! Typename of coefficients
    typedef Scalar coeff_type;

    //! Typename of scalars
    typedef Scalar scalar_type;

    //! Clone a coefficient
    static inline Teuchos::RCP<coeff_type> clone(const coeff_type& c) {
      return Teuchos::rcp(new coeff_type(c));
    }

    //! Copy a coefficient
    static inline void copy(const coeff_type& x, coeff_type* y) {
      *y = x;
    }

    //! Assign a scalar to a coefficient
    static inline void assign(coeff_type* y, const scalar_type& alpha) {
     *y = alpha;
    }

    //! y = x + beta*y
    static inline void update(coeff_type* y, const coeff_type& x,
			      const scalar_type& beta) {
      *y = x + beta*(*y);
    }

  }; // class PolynomialTraits

} // end namespace Teuchos

#endif  // TEUCHOS_POLYNOMIAL_TRAITS_HPP
