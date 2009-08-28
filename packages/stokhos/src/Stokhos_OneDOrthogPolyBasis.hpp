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

#ifndef STOKHOS_ONEDORTHOGPOLYBASIS_HPP
#define STOKHOS_ONEDORTHOGPOLYBASIS_HPP

#include <ostream>
#include <string>
#include "Stokhos_Polynomial.hpp"
#include "Stokhos_Dense3Tensor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  //! Base class for 1-D orthogonal polynomials
  template <typename ordinal_type, typename value_type>
  class OneDOrthogPolyBasis {
  public:

    //! Constructor
    OneDOrthogPolyBasis() {};

    //! Destructor
    virtual ~OneDOrthogPolyBasis() {};

    //! Return order of basis
    virtual ordinal_type order() const = 0;

    //! Return total size of basis
    virtual ordinal_type size() const = 0;

    //! Compute norm squared of each basis element
    virtual const Teuchos::Array<value_type>& norm_squared() const = 0;

    //! Compute norm squared of ith element
    virtual const value_type& norm_squared(ordinal_type i) const = 0;

    //! Compute triple product tensor
    virtual Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> > getTripleProductTensor() const = 0;

    //! Compute derivative double product tensor
    virtual Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> > getDerivDoubleProductTensor() const = 0;

    //! Project a polynomial into this basis
    virtual void projectPoly(const Polynomial<value_type>& poly, 
                             Teuchos::Array<value_type>& coeffs) const = 0;

    //! Project product of two basis polynomials into this basis
    virtual void projectProduct(ordinal_type i, ordinal_type j,
                                Teuchos::Array<value_type>& coeffs) const = 0;

    //! Project derivative of basis polynomial into this basis
    virtual void projectDerivative(ordinal_type i, 
                                   Teuchos::Array<value_type>& coeffs) const = 0;

    //! Evaluate basis polynomial at zero
    virtual value_type evaluateZero(ordinal_type i) const = 0;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const value_type& point,
                               Teuchos::Array<value_type>& basis_pts) const = 0;

    //! Print basis
    virtual void print(std::ostream& os) const = 0;

    //! Return name of basis
    virtual const std::string& getName() const = 0;

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const = 0;

    //! Get sparse grid rule number
    virtual ordinal_type getRule() const = 0;

    //! Get quadrature weight factor
    virtual value_type getQuadWeightFactor() const = 0;

    //! Get quadrature point factor
    virtual value_type getQuadPointFactor() const = 0;

  private:

    // Prohibit copying
    OneDOrthogPolyBasis(const OneDOrthogPolyBasis&);

    // Prohibit Assignment
    OneDOrthogPolyBasis& operator=(const OneDOrthogPolyBasis& b);
  

  }; // class OrthogPolyBasis

  template <typename ordinal_type, typename value_type> 
  std::ostream& 
  operator << (std::ostream& os, 
	       const OneDOrthogPolyBasis<ordinal_type, value_type>& b) {
    b.print(os);
    return os;
  }

} // Namespace Stokhos

#endif
