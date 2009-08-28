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

#ifndef STOKHOS_LEGENDREBASIS_HPP
#define STOKHOS_LEGENDREBASIS_HPP

#include "Stokhos_OneDOrthogPolyBasisBase.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class LegendreBasis : 
    public OneDOrthogPolyBasisBase<ordinal_type, value_type> {
  public:

    //! Constructor
    LegendreBasis(ordinal_type p);

    //! Destructor
    ~LegendreBasis();

    //! Project a polynomial into this basis
    void projectPoly(const Polynomial<value_type>& poly, 
		     Teuchos::Array<value_type>& coeffs) const;

    //! Project derivative of basis polynomial into this basis
    void projectDerivative(ordinal_type i, 
			   Teuchos::Array<value_type>& coeffs) const;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const value_type& point,
                               Teuchos::Array<value_type>& basis_pts) const;

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

    //! Get sparse grid rule number
    virtual ordinal_type getRule() const { return 1; }

    //! Get quadrature weight factor
    virtual value_type getQuadWeightFactor() const { return 0.5; }

    //! Get quadrature point factor
    virtual value_type getQuadPointFactor() const { return 1.0; }

  private:

    // Prohibit copying
    LegendreBasis(const LegendreBasis&);

    // Prohibit Assignment
    LegendreBasis& operator=(const LegendreBasis& b);

  protected:

    //! Derivative coefficients
    Teuchos::Array< Teuchos::Array<value_type> > deriv_coeffs;

  }; // class LegendreBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_LegendreBasisImp.hpp"

#endif
