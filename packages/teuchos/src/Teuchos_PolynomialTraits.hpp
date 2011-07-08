// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
