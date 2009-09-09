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
// Questions? Contact Christopher W. Miller (cmiller@math.umd.edu).
// 
// ***********************************************************************
// @HEADER

// ************************************************************************
// Class attempts to generate an orthogonal polynomial basis for a given
// weight function and interval.  The Discritized Stilges procidure described
// in described in "On the Calculation of Rys Polynomials and Quadratures",
// Robin P. Sagar, Vedene H. Smith is used to generate the recurrence 
// coefficients.
// Please be aware that this method is not fullproof and that it appears to
// encounter trouble with strongly singular weights since Gaussan quadrature
// is used to compute the relavant integrals.  For 'nice' weight functions the
// method seems relatively robust.
// *************************************************************************
#ifndef STOKHOS_DISCRETIZEDSTIELTJESBASIS_HPP
#define STOKHOS_DISCRETIZEDSTIELTJESBASIS_HPP

#include "Stokhos_RecurrenceBasis.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class DiscretizedStieltjesBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:
    
    //! Constructor
    DiscretizedStieltjesBasis(const std::string& label, const ordinal_type& p, 
			      const value_type (*weightFn)(const value_type),
			      const value_type& leftEndPt,
			      const value_type& rightEndPt,
			      const bool normalize);
    
    //! Destructor
    ~DiscretizedStieltjesBasis();

    //! Get sparse grid rule number
    virtual ordinal_type getRule() const { return 10; }

    //! Get quadrature weight factor
    virtual value_type getQuadWeightFactor() const { return 1; }

    //! Get quadrature point factor
    virtual value_type getQuadPointFactor() const { return 1; }
    
    //!Evaluate inner product of two basis functions to test orthogonality.
    value_type eval_inner_product(const ordinal_type& order1, 
				  const ordinal_type& order2) const;

  protected:

    //! Compute recurrence coefficients
    virtual void 
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta) const;

    //! Evaluate recurrence formula
    value_type evaluateRecurrence(const value_type& point, 
				  ordinal_type order, 
				  const Teuchos::Array<value_type>& alpha,
				  const Teuchos::Array<value_type>& beta) const;

    //! Evaluates the scaled weight function.
    value_type evaluateWeight(const value_type& x) const;

    //! Quadrature functions for generating recurrance coefficients.
    value_type 
    expectedValue_tJ_nsquared(const ordinal_type& order, 
			      const Teuchos::Array<value_type>& alpha,
			      const Teuchos::Array<value_type>& beta) const;
    value_type 
    expectedValue_J_nsquared(const ordinal_type& order, 
			     const Teuchos::Array<value_type>& alpha,
			     const Teuchos::Array<value_type>& beta) const;
    
  private:

    // Prohibit copying
    DiscretizedStieltjesBasis(const DiscretizedStieltjesBasis&);

    // Prohibit Assignment
    DiscretizedStieltjesBasis& operator=(const DiscretizedStieltjesBasis& b);
    
  protected:

    //! Scale for the weight
    mutable value_type scaleFactor;

    //! Domain Params
    const value_type leftEndPt_;
    const value_type rightEndPt_;

    //!Weight function
    const value_type (*weightFn_)(const value_type);

    //! Quadrature data for discretized stieltjes procedure
    Teuchos::Array<value_type> quad_points;
    Teuchos::Array<value_type> quad_weights;
    Teuchos::Array<Teuchos::Array< value_type > > quad_values;
    
    
  }; // class DiscretizedStieltjesBasis

} // Namespace Stokhos

// Include template definitions

#include "Stokhos_DiscretizedStieltjesBasisImp.hpp"
#endif
