// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    DiscretizedStieltjesBasis(
      const std::string& name, const ordinal_type& p, 
      value_type (*weightFn)(const value_type&),
      const value_type& leftEndPt,
      const value_type& rightEndPt,
      bool normalize, 
      GrowthPolicy growth = SLOW_GROWTH);

    //! Destructor
    ~DiscretizedStieltjesBasis();
    
    //!Evaluate inner product of two basis functions to test orthogonality.
    value_type eval_inner_product(const ordinal_type& order1, 
				  const ordinal_type& order2) const;

    /*! 
     * \brief Clone this object with the option of building a higher order
     * basis.
     */
    /*!
     * This method is following the Prototype pattern (see Design Pattern's textbook).
     * The slight variation is that it allows the order of the polynomial to be modified,
     * otherwise an exact copy is formed. The use case for this is creating basis functions
     * for column indices in a spatially varying adaptive refinement context.
     */
    virtual Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > cloneWithOrder(ordinal_type p) const;

  protected:

    //! \name Implementation of Stokhos::RecurrenceBasis methods
    //@{ 

    //! Compute recurrence coefficients
    virtual bool
    computeRecurrenceCoefficients(ordinal_type n,
				  Teuchos::Array<value_type>& alpha,
				  Teuchos::Array<value_type>& beta,
				  Teuchos::Array<value_type>& delta,
				  Teuchos::Array<value_type>& gamma) const;

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

  protected:

    //! Copy constructor with specified order
    DiscretizedStieltjesBasis(const ordinal_type& p, 
			      const DiscretizedStieltjesBasis& basis);
    
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
