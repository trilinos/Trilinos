// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_STIELTJESPCEBASIS_HPP
#define STOKHOS_STIELTJESPCEBASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "Stokhos_RecurrenceBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  /*! 
   * \brief Generates three-term recurrence using the Discretized Stieltjes 
   * procedure applied to a polynomial chaos expansion in another basis.
   */
  template <typename ordinal_type, typename value_type>
  class StieltjesPCEBasis : 
    public RecurrenceBasis<ordinal_type, value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param pce polynomial chaos expansion defining new density function
     * \param quad quadrature data for basis of PC expansion
     * \param use_pce_quad_points whether to use quad to define quadrature
     *        points for the new basis, or whether to use the Golub-Welsch
     *        system.
     * \param normalize whether polynomials should be given unit norm
     */
    StieltjesPCEBasis(
      ordinal_type p,
      const Teuchos::RCP<const Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
      const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
      bool use_pce_quad_points,
      bool normalize = false,
      bool project_integrals = false,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk = Teuchos::null);

    //! Destructor
    ~StieltjesPCEBasis();

    //! \name Implementation of Stokhos::OneDOrthogPolyBasis methods
    //@{

    //! Get Gauss quadrature points, weights, and values of basis at points
    virtual void 
    getQuadPoints(ordinal_type quad_order,
		  Teuchos::Array<value_type>& points,
		  Teuchos::Array<value_type>& weights,
		  Teuchos::Array< Teuchos::Array<value_type> >& values) const;

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

    //! Map expansion coefficients from this basis to original
    void transformCoeffsFromStieltjes(const value_type *in, 
				      value_type *out) const;

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

    //! Setup basis after computing recurrence coefficients
    virtual void setup();

    //@}

    //! Compute 3-term recurrence using Stieljtes procedure
    void stieltjes(ordinal_type nstart,
		   ordinal_type nfinish,
		   const Teuchos::Array<value_type>& weights,
		   const Teuchos::Array<value_type>& points,
		   Teuchos::Array<value_type>& a,
		   Teuchos::Array<value_type>& b,
		   Teuchos::Array<value_type>& nrm,
		   Teuchos::Array< Teuchos::Array<value_type> >& phi_vals) const;

    /*! 
     * \brief Compute \f$\int\psi^2_k(t) d\lambda(t)\f$ and 
     * \f$\int t\psi^2_k(t) d\lambda(t)\f$
     */
    void integrateBasisSquared(
      ordinal_type k, 
      const Teuchos::Array<value_type>& a,
      const Teuchos::Array<value_type>& b,
      const Teuchos::Array<value_type>& weights,
      const Teuchos::Array<value_type>& points,
      Teuchos::Array< Teuchos::Array<value_type> >& phi_vals,
      value_type& val1, value_type& val2) const;

    /*! 
     * \brief Compute \f$\int\psi^2_k(t) d\lambda(t)\f$ and 
     * \f$\int t\psi^2_k(t) d\lambda(t)\f$ by projecting onto original
     * PCE basis
     */
    void integrateBasisSquaredProj(
      ordinal_type k, 
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

    //! Copy constructor with specified order
    StieltjesPCEBasis(ordinal_type p, const StieltjesPCEBasis& basis);

  private:

    // Prohibit copying
    StieltjesPCEBasis(const StieltjesPCEBasis&);

    // Prohibit Assignment
    StieltjesPCEBasis& operator=(const StieltjesPCEBasis& b);

  protected:

    //! PC expansion
    Teuchos::RCP<const Stokhos::OrthogPolyApprox<ordinal_type, value_type> > pce;

    //! Quadrature object
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > quad;

    //! PCE quadrature weights
    const Teuchos::Array<value_type>& pce_weights;

    //! Values of PCE basis functions at quadrature points
    const Teuchos::Array< Teuchos::Array<value_type> >& basis_values;

    //! Values of PCE at quadrature points
    mutable Teuchos::Array<value_type> pce_vals;

    //! Values of generated polynomials at PCE quadrature points
    mutable Teuchos::Array< Teuchos::Array<value_type> > phi_vals;

    //! Use underlying pce's quadrature data
    bool use_pce_quad_points;

    //! Matrix mapping coefficients in Stieltjes basis back to original basis
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> fromStieltjesMat;

    //! Project Stieltjes integrals to original PCE basis
    bool project_integrals;

    //! PCE basis (needed for integral projection method)
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Short-hand for triple product
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

    //! Triple product tensor (needed for integral projection method)
    Teuchos::RCP<const Cijk_type> Cijk;

    //! Array store PC expansion of generated orthogonal polynomials
    mutable Teuchos::Array<value_type> phi_pce_coeffs;

  }; // class StieltjesPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_StieltjesPCEBasisImp.hpp"

#endif
