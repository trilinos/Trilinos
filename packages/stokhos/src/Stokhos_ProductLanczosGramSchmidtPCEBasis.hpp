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

#ifndef STOKHOS_PRODUCT_LANCZOS_GRAM_SCHMIDT_PCE_BASIS_HPP
#define STOKHOS_PRODUCT_LANCZOS_GRAM_SCHMIDT_PCE_BASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Stokhos_ReducedPCEBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief Generate a basis from a given set of PCE expansions that is 
   * orthogonal with respect to the product measure induced by these expansions.
   */
  /*!
   * Given the PCE expansions, first build a an orthogonal basis for each
   * compnent, then form the multivariate basis by a total-order tensor 
   * product.  The resulting basis is not necessarily orthogonal with respect
   * to the full measure.
   */
  template <typename ordinal_type, typename value_type>
  class ProductLanczosGramSchmidtPCEBasis : 
    public ReducedPCEBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param pce polynomial chaos expansions defining new measure
     * \param Cijk sparse triple product tensor for basis defining pce
     */
    ProductLanczosGramSchmidtPCEBasis(
     ordinal_type p,
     const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
     const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
     const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
     const Teuchos::ParameterList& params = Teuchos::ParameterList());

    //! Destructor
    virtual ~ProductLanczosGramSchmidtPCEBasis();

    //! \name Implementation of Stokhos::OrthogPolyBasis methods
    //@{

    //! Return order of basis
    ordinal_type order() const;

    //! Return dimension of basis
    ordinal_type dimension() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\Psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is size()-1.
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is size()-1 and \f$k=0,\dots,p\f$ where \f$p\f$
     * is the supplied \c order.
     */
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor(ordinal_type order) const;

    //! Evaluate basis polynomial \c i at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to size size()-1.
     */
    virtual void evaluateBases(
      const Teuchos::ArrayView<const value_type>& point,
      Teuchos::Array<value_type>& basis_vals) const;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const;

    //! Return string name of basis
    virtual const std::string& getName() const;

    //@}

    //! \name Implementation of Stokhos::ReducedPCEBasis methods
    //@{

    //! Transform coefficients to original basis from this basis
    virtual void 
    transformToOriginalBasis(const value_type *in, 
			     value_type *out,
			     ordinal_type ncol = 1, 
			     bool transpose = false) const;

    //! Transform coefficients from original basis to this basis
    virtual void 
    transformFromOriginalBasis(const value_type *in, 
			       value_type *out,
			       ordinal_type ncol = 1, 
			       bool transpose = false) const;

    //! Get reduced quadrature object
    virtual Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
    getReducedQuadrature() const;

    //@}

  protected:

    // Determine if a pce is linear, in that it has a total degree of at
    // most 1.  If the pce is nonlinear, return -2, or if it is constant, 
    // return -1, otherwise return the index of the variable the pce is
    // linear in, ie, if return value is i, the pce = a_{i+1}*\xi_i
    ordinal_type
    isInvariant(const Stokhos::OrthogPolyApprox<ordinal_type, value_type>& pce) const;

  private:

    // Prohibit copying
    ProductLanczosGramSchmidtPCEBasis(const ProductLanczosGramSchmidtPCEBasis&);

    // Prohibit Assignment
    ProductLanczosGramSchmidtPCEBasis& operator=(const ProductLanczosGramSchmidtPCEBasis&);
    
  protected:

    typedef Teuchos::SerialDenseVector<ordinal_type,value_type> SDV;
    typedef Teuchos::SerialDenseMatrix<ordinal_type,value_type> SDM;

    //! Name of basis
    std::string name;

    //! Algorithm parameters
    Teuchos::ParameterList params;

    //! Size of original pce basis
    ordinal_type pce_sz;
    
    //! Total order of basis
    ordinal_type p;

    //! Total dimension of basis
    ordinal_type d;

    //! Total size of basis
    ordinal_type sz;

    //! Product Lanczos basis
    Teuchos::RCP< Stokhos::CompletePolynomialBasis<ordinal_type,value_type> > tensor_lanczos_basis;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Values of transformed basis at quadrature points
    SDM Q;

    //! Coefficients of transformed basis in original basis
    SDM Qp;

    //! Reduced quadrature object
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > reduced_quad;

    //! Temporary pce used in invariant subspace calculations
    mutable Stokhos::OrthogPolyApprox<ordinal_type, value_type> tmp_pce;

  }; // class ProductLanczosGramSchmidtPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_ProductLanczosGramSchmidtPCEBasisImp.hpp"

#endif
