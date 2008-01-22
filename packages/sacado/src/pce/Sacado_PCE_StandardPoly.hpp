// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PCE_STANDARDPOLY_HPP
#define SACADO_PCE_STANDARDPOLY_HPP

#include <vector>
#include <ostream>

namespace Sacado {

  namespace PCE {

    template <typename T>
    class StandardPoly {
    public:

      //! Constructor with all zero coefficients
      StandardPoly(unsigned int deg);

      //! Constructor with specified coefficients
      StandardPoly(const std::vector<T>& coefficients);

      //! Copy constructor
      StandardPoly(const StandardPoly& p);

      //! Destructor
      ~StandardPoly();

      //! Assignment
      StandardPoly& operator=(const StandardPoly& p);

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
		    const StandardPoly<T>& a,
		    const StandardPoly<T>& b, 
		    const T& beta);

      //! Add two polynomials and put into this
      /*!
       * Sets this = alpha*a* + gamma*this;
       */
      void add(const T& alpha, 
	       const StandardPoly<T>& a,
	       const T& gamma);

      //! Add two polynomials and put into this
      /*!
       * Sets this = alpha*a* + beta*b + gamma*this;
       */
      void add(const T& alpha, 
	       const StandardPoly<T>& a,
	       const T& beta,
	       const StandardPoly<T>& b, 
	       const T& gamma);

      void print(std::ostream& os) const;

    protected:

      //! Vector of coefficients
      std::vector<T> coeffs;

    }; // class StandardPoly

    template <typename T> 
    std::ostream& operator << (std::ostream& os, const StandardPoly<T>& p) {
      p.print(os);
      return os;
    }

  } // Namespace PCE

} // Namespace Sacado

// Include template definitions
#include "Sacado_PCE_StandardPolyImp.hpp"

#endif // SACADO_PCE_STANDARDPOLY_HPP
