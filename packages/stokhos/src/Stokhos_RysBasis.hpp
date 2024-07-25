// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_RYSBASIS_HPP
#define STOKHOS_RYSBASIS_HPP

#include "Stokhos_DiscretizedStieltjesBasis.hpp"

namespace Stokhos {

  //! Rys polynomial basis
  /*!
   * Rys polynomials are polynomials orthogonal with respect to the weight
   * function
   * \f[
   *     w(x) = e^{\frac{-x^2}{2}}
   * \f]
   * defined on the interval \f$[-c,c]\f$, for a given choice of \f$c\f$.  The
   * corresponding density \f$\rho(x)\f$ is obtained by scaling \f$w(x)\f$ to
   * unit probability.
   *
   * The coefficients of the corresponding three-term recursion are generated
   * using the Discretized Stieltjes procedure implemented by
   * Stokhos::DiscretizedStieltjesBasis.
   */
  template <typename ordinal_type, typename value_type>
  class RysBasis : 
    public DiscretizedStieltjesBasis<ordinal_type, value_type> {
  public:
    
    //! Constructor
    /*!
     * \param p order of the basis
     * \param c defines domain of support of weight function
     * \param normalize whether polynomials should be given unit norm
     */
    RysBasis(ordinal_type p, value_type c, bool normalize, 
	     GrowthPolicy growth = SLOW_GROWTH) :
      DiscretizedStieltjesBasis<ordinal_type,value_type>(
	"Rys", p, rysWeight, -c, c, normalize, growth) {}
    
    //! Destructor
    ~RysBasis() {}

    //! The Rys weight function
    static value_type rysWeight(const value_type& x) { 
      return std::exp(-x*x/2.0); 
    }

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
    virtual Teuchos::RCP<OneDOrthogPolyBasis<ordinal_type,value_type> > cloneWithOrder(ordinal_type p) const
    { return Teuchos::rcp(new RysBasis<ordinal_type,value_type>(p,*this)); }

  protected:

    //! Copy constructor with specified order
    RysBasis(ordinal_type p, const RysBasis& basis) : 
      DiscretizedStieltjesBasis<ordinal_type,value_type>(p, basis) {}

  private:

    // Prohibit copying
    RysBasis(const RysBasis&);

    // Prohibit Assignment
    RysBasis& operator=(const RysBasis& b);
    
  }; // class RysBasis

} // Namespace Stokhos

#endif
