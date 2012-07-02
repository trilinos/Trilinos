// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_MONOMIAL_PROJ_GRAM_SCHMIDT_SIMPLEX_PCE_BASIS_HPP
#define STOKHOS_MONOMIAL_PROJ_GRAM_SCHMIDT_SIMPLEX_PCE_BASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"
#include "Stokhos_ConfigDefs.h"

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
  class MonomialProjGramSchmidtSimplexPCEBasis : 
    public ProductBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param p order of the basis
     * \param pce polynomial chaos expansions defining new measure
     * \param quad quadrature data for basis defining pce
     * \param Cijk sparse triple product tensor for basis defining pce
     * \param sparse_tol tolerance for dropping terms in sparse tensors
     */
    MonomialProjGramSchmidtSimplexPCEBasis(
     ordinal_type p,
     const Teuchos::Array< Stokhos::OrthogPolyApprox<ordinal_type, value_type> >& pce,
     const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad,
     const Teuchos::ParameterList& params = Teuchos::ParameterList());

    //! Destructor
    virtual ~MonomialProjGramSchmidtSimplexPCEBasis();

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
    virtual void evaluateBases(const Teuchos::Array<value_type>& point,
			       Teuchos::Array<value_type>& basis_vals) const;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const;

    //! Return string name of basis
    virtual const std::string& getName() const;

    //@}

    //! \name Implementation of Stokhos::ProductBasis methods
    //@{

    //! Get orders of each coordinate polynomial given an index \c i
    /*!
     * The returned array is of size \f$d\f$, where \f$d\f$ is the dimension of
     * the basis, and entry \f$l\f$ is given by \f$i_l\f$ where
     * \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual Teuchos::Array<ordinal_type> getTerm(ordinal_type i) const;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type 
    getIndex(const Teuchos::Array<ordinal_type>& term) const;

    //! Return coordinate bases
    /*!
     * Array is of size dimension().
     */
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, 
							   value_type> > > 
    getCoordinateBases() const;

    //@}

    //! Get reduced quadrature object
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
    getReducedQuadrature() const;

    //! Compute transformed PCE representation in this basis
    void computeTransformedPCE(
      ordinal_type i,
      Stokhos::OrthogPolyApprox<ordinal_type,value_type>& pce) const;

    //! Transform coefficients to original basis from this basis
    void transformToOriginalBasis(const value_type *in, value_type *out) const;

    //! Transform coefficients from original basis to this basis
    void transformFromOriginalBasis(const value_type *in, value_type *out) const;

    //! Get transformation matrix
    const Teuchos::SerialDenseMatrix<ordinal_type,value_type>&
    getTransformationMatrix() const { return Qpp; }

    //! Get reduced basis evaluated at original quadrature points
    void getBasisAtOriginalQuadraturePoints(
      Teuchos::Array< Teuchos::Array<double> >& red_basis_vals) const;

  protected:
    
    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension
     */
    void compute_terms(ordinal_type p, ordinal_type d, ordinal_type& sz,
		       Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
		       Teuchos::Array<ordinal_type>& num_terms);

    /*!
     * \brief Compute basis index given the orders for each basis
     * dimension.
     */
    ordinal_type compute_index(const Teuchos::Array<ordinal_type>& terms) const;

    void reducedQuadrature_QRCP(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values);

    /*
    void reducedQuadrature_CS(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values);
    */

    void reducedQuadrature_GLPK(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Q,
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& F,
      const Teuchos::Array<value_type>& weights,
      Teuchos::RCP< Teuchos::Array<value_type> >& reduced_weights,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_points,
      Teuchos::RCP< Teuchos::Array< Teuchos::Array<value_type> > >& reduced_values);

    ordinal_type computeRank(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& R,
      const value_type tol) const;

  private:

    // Prohibit copying
    MonomialProjGramSchmidtSimplexPCEBasis(const MonomialProjGramSchmidtSimplexPCEBasis&);

    // Prohibit Assignment
    MonomialProjGramSchmidtSimplexPCEBasis& operator=(const MonomialProjGramSchmidtSimplexPCEBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Original quadrature object
    Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> > quad;

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

    //! 2-D array of basis terms
    Teuchos::Array< Teuchos::Array<ordinal_type> > terms;

    //! Number of terms up to each order
    Teuchos::Array<ordinal_type> num_terms;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Norms of input pce's
    Teuchos::Array<value_type> pce_norms;

    //! Values of PCE basis functions at quadrature points
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> A;

    //! Values of transformed basis at quadrature points
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> Q;

    //! Coefficients of transformed basis in original basis
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> Qpp;

    //! Transition matrix between monomial and "Q" basis
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> R;

    Teuchos::Array<value_type> sigma;

    //! V^T from SVD
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> Vt;

    //! Reduced quadrature object
    Teuchos::RCP<const Stokhos::UserDefinedQuadrature<ordinal_type, value_type> > reduced_quad;

    Teuchos::LAPACK<ordinal_type,value_type> lapack;
    Teuchos::BLAS<ordinal_type,value_type> blas;

    //! Whether to print a bunch of stuff out
    bool verbose;

    //! Rank threshold
    value_type rank_threshold;

    //! Dimension reduction tolerance
    value_type reduction_tol;

    //! Orthogonalization method
    std::string orthogonalization_method;

    Teuchos::Array<ordinal_type> piv;

  }; // class MonomialProjGramSchmidtSimplexPCEBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_MonomialProjGramSchmidtSimplexPCEBasisImp.hpp"

#endif
