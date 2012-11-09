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

#ifndef STOKHOS_GRAMSCHMIDTBASIS_HPP
#define STOKHOS_GRAMSCHMIDTBASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  /*! 
   * \brief Transforms a non-orthogonal multivariate basis to an orthogonal one 
   * using the Gram-Schmit procedure.
   */
  /*!
   * Given a basis \f$\{\Psi_i\}\f$ with an inner product defined by
   * \f[
   *     (\Psi_i,\Psi_j) = \sum_{k=0}^Q w_k\Psi_i(x_k)\Psi_j(x_k)
   * \f]
   * where \f$\{x_k\}\f$ and \f$\{w_k\}\f$ are a set of \f$Q\f$ quadrature 
   * points and weights, this class generates a new basis 
   * \f$\{\tilde{\Psi}_i\}\f$ that satisfies 
   * \f$ (\Psi_i,\Psi_j) = \delta_{ij}\f$.
   *
   * NOTE:  Currently on the classical Gram-Schmidt algorithm is
   * implemented.
   */
  template <typename ordinal_type, typename value_type>
  class GramSchmidtBasis : 
    public OrthogPolyBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param basis basis defining \f$\{\Psi_i\}\f$
     * \param points quadrature points f$\{x_k\}\f$
     * \param weights quadrature weights \f$\{w_k\}\f$
     * \param sparse_tol tolerance for dropping terms in sparse tensors
     */
    GramSchmidtBasis(
     const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis,
     const Teuchos::Array< Teuchos::Array<value_type> >& points,
     const Teuchos::Array<value_type>& weights,
     const value_type& sparse_tol = 1.0e-15);

    //! Destructor
    virtual ~GramSchmidtBasis();

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

    //! Transform coefficients from original basis to this basis
    void transformCoeffs(const value_type *in, value_type *out) const;

  private:

    // Prohibit copying
    GramSchmidtBasis(const GramSchmidtBasis&);

    // Prohibit Assignment
    GramSchmidtBasis& operator=(const GramSchmidtBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Original basis (not orthogonal w.r.t. inner product)
    Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Quadrature weights defining inner product
    Teuchos::Array<value_type> weights;

    //! Quadrature basis values defining inner product
    Teuchos::Array< Teuchos::Array<value_type> > basis_values;

    //! Tolerance for computing sparse Cijk
    value_type sparse_tol;

    //! Total order of basis
    ordinal_type p;

    //! Total dimension of basis
    ordinal_type d;

    //! Total size of basis
    ordinal_type sz;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Matrix storing gram-schmidt coefficients
    Teuchos::SerialDenseMatrix<ordinal_type, value_type> gs_mat;

    //! Temporary array for basis evaluation
    mutable Teuchos::Array<value_type> basis_vals_tmp;

  }; // class GramSchmidtBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_GramSchmidtBasisImp.hpp"

#endif
