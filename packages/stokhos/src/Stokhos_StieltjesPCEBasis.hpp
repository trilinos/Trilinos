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

#ifndef STOKHOS_STIELTJESPCEBASIS_HPP
#define STOKHOS_STIELTJESPCEBASIS_HPP

#include "Stokhos_RecurrenceBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_Quadrature.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class StieltjesPCEBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    StieltjesPCEBasis(
       ordinal_type p,
       const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce,
       const Stokhos::Quadrature<ordinal_type, value_type>& quad,
       bool use_pce_quad_points,
       bool normalize = false);

    //! Destructor
    ~StieltjesPCEBasis();

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

    //! Get sparse grid rule number
    virtual ordinal_type getRule() const;

    //! Get quadrature weight factor
    virtual value_type getQuadWeightFactor() const;

    //! Get quadrature point factor
    virtual value_type getQuadPointFactor() const;

    void transformCoeffsFromStieltjes(const value_type *in, value_type *out) const;

  protected:

    //! Compute recurrence coefficients
    virtual void 
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta) const;

    //! Compute 3-term recurrence using Stieljtes procedure
    void stieltjes(ordinal_type nstart,
		   ordinal_type nfinish,
		   const Teuchos::Array<value_type>& weights,
		   const Teuchos::Array<value_type>& points,
		   Teuchos::Array<value_type>& a,
		   Teuchos::Array<value_type>& b,
		   Teuchos::Array<value_type>& nrm,
		   Teuchos::Array< Teuchos::Array<value_type> >& phi_vals) const;

    //! Compute \int \phi^2_k(t) d\lambda(t) and \int t\phi^2_k(t) d\lambda(t)
    void integrateBasisSquared(ordinal_type k, 
			       const Teuchos::Array<value_type>& a,
			       const Teuchos::Array<value_type>& b,
			       const Teuchos::Array<value_type>& weights,
			       const Teuchos::Array<value_type>& points,
			       Teuchos::Array< Teuchos::Array<value_type> >& phi_vals,
			       value_type& val1, value_type& val2) const;

    //! Evaluate polynomials via 3-term recurrence
    void evaluateRecurrence(ordinal_type k,
			    const Teuchos::Array<value_type>& a,
			    const Teuchos::Array<value_type>& b,
			    const Teuchos::Array<value_type>& points,
			    Teuchos::Array< Teuchos::Array<value_type> >& values) const;

  private:

    // Prohibit copying
    StieltjesPCEBasis(const StieltjesPCEBasis&);

    // Prohibit Assignment
    StieltjesPCEBasis& operator=(const StieltjesPCEBasis& b);

  protected:

    //! Values of PCE at quadrature points
    mutable Teuchos::Array<value_type> pce_vals;

    //! PCE quadrature weights
    Teuchos::Array<value_type> pce_weights;

    //! Values of generated polynomials at PCE quadrature points
    mutable Teuchos::Array< Teuchos::Array<value_type> > phi_vals;

    //! Use underlying pce's quadrature data
    bool use_pce_quad_points;

    //! Matrix mapping coefficients in Stieltjes basis back to original basis
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> fromStieltjesMat;

  }; // class StieltjesPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_StieltjesPCEBasisImp.hpp"

#endif
