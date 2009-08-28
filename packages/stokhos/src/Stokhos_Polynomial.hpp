// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_POLYNOMIAL_HPP
#define STOKHOS_POLYNOMIAL_HPP

#include <ostream>
#include "Teuchos_Array.hpp"

namespace Stokhos {

  template <typename T>
  class Polynomial {
  public:

    //! Constructor with all zero coefficients
    Polynomial(unsigned int deg);
    
    //! Constructor with specified coefficients
    Polynomial(const Teuchos::Array<T>& coefficients);
    
    //! Copy constructor
    Polynomial(const Polynomial& p);
    
    //! Destructor
    ~Polynomial();
    
    //! Assignment
    Polynomial& operator=(const Polynomial& p);
    
    //! Return degree
    unsigned int degree() const;
    
    //! Get coefficient
    const T& coeff(unsigned int i) const;
    
    //! Get coefficient
    T& coeff(unsigned int i);
    
    //! Get coefficient
    const T& operator[](unsigned int i) const;
    
    //! Get coefficient
    T& operator[](unsigned int i);

    //! Multiply two polynomials and put into this
    /*!
     * Sets this = alpha*a*b + beta*this
     */
    void multiply(const T& alpha, 
		  const Polynomial<T>& a,
		  const Polynomial<T>& b, 
		  const T& beta);
    
    //! Add two polynomials and put into this
    /*!
     * Sets this = alpha*a* + gamma*this;
     */
    void add(const T& alpha, 
	     const Polynomial<T>& a,
	     const T& gamma);

    //! Add two polynomials and put into this
    /*!
     * Sets this = alpha*a* + beta*b + gamma*this;
     */
    void add(const T& alpha, 
	     const Polynomial<T>& a,
	     const T& beta,
	     const Polynomial<T>& b, 
	     const T& gamma);

    void print(std::ostream& os) const;
    
  protected:
    
    //! Vector of coefficients
    Teuchos::Array<T> coeffs;
    
  }; // class Polynomial

  template <typename T> 
  std::ostream& operator << (std::ostream& os, const Polynomial<T>& p) {
    p.print(os);
    return os;
  }

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_PolynomialImp.hpp"

#endif // STOKHOS_POLYNOMIAL_HPP
