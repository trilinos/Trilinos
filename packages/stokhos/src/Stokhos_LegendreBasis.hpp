// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_LEGENDREBASIS_HPP
#define STOKHOS_LEGENDREBASIS_HPP

#include "Stokhos_RecurrenceBasis.hpp"

namespace Stokhos {

  //! Legendre polynomial basis
  /*!
   * Legendre polynomials are defined by the recurrence relationship
   * \f[
   *   \psi_{k+1}(x) = \frac{2k+1}{k+1}x\psi_{k}(x) - \frac{k}{k+1}\psi_{k-1}(x)
   * \f]
   * with \f$\psi_{-1}(x) = 0\f$ and \f$\psi_{0}(x) = 1\f$.  The corresponding
   * density function is 
   * \f[
   *   \rho(x) = \frac{1}{2}, \quad x\in[-1,1].
   * \f]
   *
   * This class implements computeRecurrenceCoefficients() using the
   * above formula.
   */
  template <typename ordinal_type, typename value_type>
  class LegendreBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param normalize whether polynomials should be given unit norm
     */
    LegendreBasis(ordinal_type p, bool normalize = false, 
		  GrowthPolicy growth = SLOW_GROWTH);

    //! Destructor
    ~LegendreBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{ 

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

    //@}

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

    //! Copy constructor with specified order
    LegendreBasis(ordinal_type p, const LegendreBasis& basis);

  private:

    // Prohibit copying
    LegendreBasis(const LegendreBasis&);

    // Prohibit Assignment
    LegendreBasis& operator=(const LegendreBasis& b);

  }; // class LegendreBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_LegendreBasisImp.hpp"

#endif
