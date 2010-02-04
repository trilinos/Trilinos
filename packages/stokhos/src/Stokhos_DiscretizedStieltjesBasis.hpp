// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef STOKHOS_DISCRETIZEDSTIELTJESBASIS_HPP
#define STOKHOS_DISCRETIZEDSTIELTJESBASIS_HPP

#include "Stokhos_RecurrenceBasis.hpp"

namespace Stokhos {

  //! Generates three-term recurrence using the Discretized Stieltjes procedure
  /*!
   * Class generates an orthogonal polynomial basis for a given
   * weight function and interval.  The Discretized Stiltjes procedure described
   * in "On the Calculation of Rys Polynomials and Quadratures",
   * Robin P. Sagar, Vedene H. Smith is used to generate the recurrence 
   * coefficients.
   *
   * Please be aware that this method is not fullproof and that it appears to
   * encounter trouble with strongly singular weights since Gaussan quadrature
   * is used to compute the relevant integrals.  For 'nice' weight functions the
   * method seems relatively robust.
   */
  template <typename ordinal_type, typename value_type>
  class DiscretizedStieltjesBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:
    
    //! Constructor
    /*!
     * \param name name of the basis displayed when printing
     * \param p order of the basis
     * \param weightFn function pointer defining weight function
     * \param leftEndPt left end point of the domain of the weight function
     * \param rightEndPt right end point of the domain of the weight function
     * \param normalize whether polynomials should be given unit norm
     */
    DiscretizedStieltjesBasis(const std::string& name, const ordinal_type& p, 
			      value_type (*weightFn)(const value_type&),
			      const value_type& leftEndPt,
			      const value_type& rightEndPt,
			      bool normalize);
    
    //! Destructor
    ~DiscretizedStieltjesBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{ 

    //! Get sparse grid rule number as defined by Dakota's \c webbur package
    /*!
     * This method is needed for building Smolyak sparse grids out of this 
     * basis.  A rule number of 10 is not defined by the webbur package, and
     * this rule number is used internally by Stokhos::SparseGridQuadrature
     * to pass an arbitrary one-dimensional basis to that package.
     */
    virtual ordinal_type getRule() const { return 10; }

    //! Get quadrature weight factor as defined by Dakota's \c webbur package
    /*!
     * This method is needed for building Smolyak sparse grids out of this 
     * basis.
     */
    virtual value_type getQuadWeightFactor() const { return 1; }

    //! Get quadrature point factor as defined by Dakota's \c webbur package
    /*!
     * This method is needed for building Smolyak sparse grids out of this 
     * basis.
     */
    virtual value_type getQuadPointFactor() const { return 1; }

    //@}
    
    //!Evaluate inner product of two basis functions to test orthogonality.
    value_type eval_inner_product(const ordinal_type& order1, 
				  const ordinal_type& order2) const;

  protected:

    //! \name Implementation of Stokhos::RecurrenceBasis methods
    //@{ 

    //! Compute recurrence coefficients
    virtual void 
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta) const;

    //@}

    //! Evaluate 3-term recurrence formula using supplied coefficients
    value_type evaluateRecurrence(const value_type& point, 
				  ordinal_type order, 
				  const Teuchos::Array<value_type>& alpha,
				  const Teuchos::Array<value_type>& beta) const;

    //! Evaluates the scaled weight function.
    value_type evaluateWeight(const value_type& x) const;

    //! Approximates \f$\langle t\psi_k(t) \rangle\f$ where \f$k\f$ = \c order
    value_type 
    expectedValue_tJ_nsquared(const ordinal_type& order, 
			      const Teuchos::Array<value_type>& alpha,
			      const Teuchos::Array<value_type>& beta) const;

    //! Approximates \f$\langle \psi_k(t) \rangle\f$ where \f$k\f$ = \c order
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

    //! Left end point of domain
    const value_type leftEndPt_;

    //! Right end point of domain
    const value_type rightEndPt_;

    //! Weight function
    value_type (*weightFn_)(const value_type&);

    //! Quadrature points for discretized stieltjes procedure
    Teuchos::Array<value_type> quad_points;

    //! Quadrature points for discretized stieltjes procedure
    Teuchos::Array<value_type> quad_weights;

    //! Quadrature values for discretized stieltjes procedure
    Teuchos::Array<Teuchos::Array< value_type > > quad_values;
    
  }; // class DiscretizedStieltjesBasis

} // Namespace Stokhos

// Include template definitions

#include "Stokhos_DiscretizedStieltjesBasisImp.hpp"
#endif
