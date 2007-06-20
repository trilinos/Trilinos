// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
