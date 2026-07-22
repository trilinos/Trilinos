// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MONOMIAL_GRAM_SCHMIDT_PCE_BASIS_HPP
#define STOKHOS_MONOMIAL_GRAM_SCHMIDT_PCE_BASIS_HPP

#include "Stokhos_GSReducedPCEBasisBase.hpp"

namespace Stokhos {

  /*! 
   * \brief Generate a basis from a given set of PCE expansions that is 
   * orthogonal with respect to the product measure induced by these expansions.
   */
  /*!
   * Given the PCE expansions, first build a non-orthogonal monomial basis.  
   * Orthogonalize this basis using Gram-Schmidt, then build a quadrature rule
   * using the simplex method.
   */
  template <typename ordinal_type, typename value_type>
  class MonomialGramSchmidtPCEBasis : 
    public GSReducedPCEBasisBase<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param pce polynomial chaos expansions defining new measure
     * \param quad quadrature data for basis defining pce
     * \param Cijk sparse triple product tensor for basis defining pce
     * \param sparse_tol tolerance for dropping terms in sparse tensors
     */
    MonomialGramSchmidtPCEBasis(
     ordinal_type p,
     const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
     const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
     const Teuchos::ParameterList& params = Teuchos::ParameterList());

    //! Destructor
    virtual ~MonomialGramSchmidtPCEBasis();

    //! \name Implementation of Stokhos::OrthogPolyBasis methods
    //@{

    //! Return string name of basis
    virtual const std::string& getName() const;

    //@}

  protected:

    //! Build the reduced basis, parameterized by total order \c max_p
    /*!
     * Returns resulting size of reduced basis
     */
    virtual ordinal_type 
    buildReducedBasis(
      ordinal_type max_p, 
      value_type threshold,
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& A, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::Array< Stokhos::MultiIndex<ordinal_type> >& terms_,
      Teuchos::Array<ordinal_type>& num_terms_,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Qp_, 
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Q_);

  private:

    // Prohibit copying
    MonomialGramSchmidtPCEBasis(const MonomialGramSchmidtPCEBasis&);

    // Prohibit Assignment
    MonomialGramSchmidtPCEBasis& operator=(const MonomialGramSchmidtPCEBasis& b);
    
  protected:

    typedef Stokhos::CompletePolynomialBasisUtils<ordinal_type,value_type> CPBUtils;
    typedef Teuchos::SerialDenseVector<ordinal_type,value_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> SDM;

    //! Name of basis
    std::string name;

  }; // class MonomialGramSchmidtPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_MonomialGramSchmidtPCEBasisImp.hpp"

#endif
