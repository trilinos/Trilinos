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
#ifndef STOKHOS_RECURRENCEBASIS_HPP
#define STOKHOS_RECURRENCEBASIS_HPP

#include "Stokhos_OneDOrthogPolyBasisBase.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class RecurrenceBasis : 
    public OneDOrthogPolyBasisBase<ordinal_type, value_type> {
  public:
    
    //! Constructor
    RecurrenceBasis(const ordinal_type& p , const char* label,
                    const value_type (*weightFn)(const value_type),
                    const value_type& leftEndPt,
                    const value_type& rightEndPt,const bool normalize);
    
    //! Destructor
    ~RecurrenceBasis();
    
    //! Project a polynomial into this basis (NOT IMPLIMENTED)
    virtual void projectPoly(const Polynomial<value_type>& poly, 
		     std::vector<value_type>& coeffs) const;

    //! Project derivative of basis polynomial into this basis (NOT IMPLIMENTED)
    virtual void projectDerivative(ordinal_type i, 
                             std::vector<value_type>& coeffs) const;
 
    //! Evaluates the scaled weight function.
    virtual value_type evaluateWeight(const value_type& x) const;


    ////! return vectors containing recurrance coefficients.
    virtual void getAlpha(std::vector<value_type>& alphaOut) const {alphaOut = this->alpha;}
    virtual void getBeta(std::vector<value_type>& betaOut) const {betaOut = this->beta;}
    virtual void getGamma(std::vector<value_type>& gammaOut) const {gammaOut = this->gamma;}

    //! Quadrature functions for generating recurrance coefficients.
    virtual value_type expectedValue_tJ_nsquared(const ordinal_type& order) const;
    virtual value_type expectedValue_J_nsquared(const ordinal_type& order) const;
    
    //!Evaluate inner product of two basis functions to test orthogonality.
    virtual value_type eval_inner_product(const ordinal_type& order1, const ordinal_type& order2) const;
    
    //!Evaluate p_th basis function at a given point.
    virtual value_type evaluateBasesOrder_p(const value_type& x, 
					const ordinal_type& order) const;

    //! Evaluate basis polynomials at given point
    virtual void evaluateBases(const value_type& point,
                               std::vector<value_type>& basis_pts) const;


    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  std::vector<value_type>& points,
		  std::vector<value_type>& weights,
		  std::vector< std::vector<value_type> >& values) const;

    //! Get sparse grid rule number
    virtual ordinal_type getRule() const { return 10; }

    //! Get quadrature weight factor
    virtual value_type getQuadWeightFactor() const { return 1; }

    //! Get quadrature point factor
    virtual value_type getQuadPointFactor() const { return 1; }


    virtual Teuchos::RCP< const Stokhos::Dense3Tensor<ordinal_type, value_type> > getTripleProductTensor() const;
    
  private:

    // Prohibit copying
    RecurrenceBasis(const RecurrenceBasis&);

    // Prohibit Assignment
    RecurrenceBasis& operator=(const RecurrenceBasis& b);
    
  protected:
    //! Scale for the weight
    value_type scaleFactor;

    //! Normalized?
    const bool normalize_;

    //! Domain Params
    const value_type leftEndPt_;
    const value_type rightEndPt_;

    //!Weight function
    const value_type (*weightFn_)(const value_type);
  
    //!label
    const char * label_;

    //! Recurrance coeffs
    std::vector<value_type> alpha;
    std::vector<value_type> beta;
    std::vector<value_type> gamma;
    
    
  }; // class RecurrenceBasis

} // Namespace Stokhos

// Include template definitions

#include "Stokhos_RecurrenceBasisImp.hpp"
#endif
