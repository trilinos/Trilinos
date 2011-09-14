//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
//@HEADER

#ifndef EPETRA_EXT_POLYNOMIAL_VECTOR_TRAITS_H
#define EPETRA_EXT_POLYNOMIAL_VECTOR_TRAITS_H

#include "Teuchos_PolynomialTraits.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

namespace Teuchos {
  
  //! Specilization of Teuchos::PolynomialTraits for %Epetra vectors.
  /*!
   * This class provides a specilization of Teuchos::PolynomialTraits for
   * Epetra_Vector vectors, allowing these vectors to be coefficients
   * in the Teuchos::Polynomial.
   */
  template <>
  class PolynomialTraits<Epetra_Vector> {
  public:

    //! Typename of coefficients
    typedef Epetra_Vector coeff_type;

    //! Typename of scalars
    typedef double scalar_type;

    //! Clone a coefficient
    static inline Teuchos::RefCountPtr<coeff_type> clone(const coeff_type& c) {
      return Teuchos::rcp(new Epetra_Vector(c));
    }

    //! Copy a coefficient
    static inline void copy(const coeff_type& x, coeff_type* y) {
      y->Scale(1.0, x);
    }

    //! Assign a scalar to a coefficient
    static inline void assign(coeff_type* y, const scalar_type& alpha) {
      y->PutScalar(alpha);
    }

    //! y = x + beta*y
    static inline void update(coeff_type* y, const coeff_type& x, 
			      const scalar_type& beta) {
      y->Update(1.0, x, beta);
    }

  }; // class PolynomialTraits<Epetra_Vector>

} // namespace Teuchos

#endif // EPETRA_EXT_POLYNOMIAL_VECTOR_TRAITS_Hx
