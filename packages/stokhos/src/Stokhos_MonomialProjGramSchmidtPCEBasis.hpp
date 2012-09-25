// @HEADER
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_MONOMIAL_PROJ_GRAM_SCHMIDT_PCE_BASIS_HPP
#define STOKHOS_MONOMIAL_PROJ_GRAM_SCHMIDT_PCE_BASIS_HPP

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
  class MonomialProjGramSchmidtPCEBasis : 
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
    MonomialProjGramSchmidtPCEBasis(
     ordinal_type p,
     const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
     const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
     const Teuchos::ParameterList& params = Teuchos::ParameterList());

    //! Destructor
    virtual ~MonomialProjGramSchmidtPCEBasis();

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
      Teuchos::Array< Teuchos::Array<ordinal_type> >& terms_,
      Teuchos::Array<ordinal_type>& num_terms_,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Qp_, 
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& Q_);

  private:

    // Prohibit copying
    MonomialProjGramSchmidtPCEBasis(const MonomialProjGramSchmidtPCEBasis&);

    // Prohibit Assignment
    MonomialProjGramSchmidtPCEBasis& operator=(const MonomialProjGramSchmidtPCEBasis& b);
    
  protected:

    typedef Stokhos::CompletePolynomialBasisUtils<ordinal_type,value_type> CPBUtils;
    typedef Teuchos::SerialDenseVector<ordinal_type,value_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> SDM;

    //! Name of basis
    std::string name;

  }; // class MonomialProjGramSchmidtPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_MonomialProjGramSchmidtPCEBasisImp.hpp"

#endif
