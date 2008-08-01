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

#ifndef STOKHOS_BASIS_HPP
#define STOKHOS_BASIS_HPP

#include <ostream>
#include <string>
#include "Stokhos_Polynomial.hpp"

namespace Stokhos {

  template <typename T>
  class OrthogPolyBasis {
  public:

    //! Typename of values
    typedef T value_type;

    //! Constructor
    OrthogPolyBasis() {};

    //! Destructor
    virtual ~OrthogPolyBasis() {};

    //! Return order of basis
    virtual unsigned int order() const = 0;

    //! Return dimension of basis
    virtual unsigned int dimension() const = 0;

    //! Return total size of basis
    virtual unsigned int size() const = 0;

    //! Compute norm squared of each basis element
    virtual const std::vector<T>& norm_squared() const = 0;

    //! Project a polynomial into this basis
    virtual void projectPoly(const Polynomial<T>& poly, 
			     std::vector<T>& coeffs) const = 0;

    //! Project product of two basis polynomials into this basis
    virtual void projectProduct(unsigned int i, unsigned int j,
				std::vector<T>& coeffs) const = 0;

    //! Project derivative of basis polynomial into this basis
    virtual void projectDerivative(unsigned int i, 
				   std::vector<T>& coeffs) const = 0;

    //! Write polynomial in standard basis
    virtual Polynomial<T> toStandardBasis(const T coeffs[], 
					  unsigned int n) const = 0;

    //! Evaluate basis polynomial at zero
    virtual T evaluateZero(unsigned int i) const = 0;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const std::vector<T>& point,
			       std::vector<T>& basis_pts) const {}

    //! Print basis
    virtual void print(std::ostream& os) const = 0;

    //! Get term
    virtual std::vector<unsigned int> getTerm(unsigned int i) const = 0;

    //! Get index
    virtual unsigned int 
    getIndex(const std::vector<unsigned int>& term) const = 0;

    //! Return name of basis
    virtual const std::string& getName() const = 0;

  private:

    // Prohibit copying
    OrthogPolyBasis(const OrthogPolyBasis&);

    // Prohibit Assignment
    OrthogPolyBasis& operator=(const OrthogPolyBasis& b);

  }; // class OrthogPolyBasis

  template <typename T> 
  std::ostream& operator << (std::ostream& os, const OrthogPolyBasis<T>& b) {
    b.print(os);
    return os;
  }

} // Namespace Stokhos

#endif
