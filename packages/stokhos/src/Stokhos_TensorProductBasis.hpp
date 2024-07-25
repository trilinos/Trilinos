// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_TENSOR_PRODUCT_BASIS_HPP
#define STOKHOS_TENSOR_PRODUCT_BASIS_HPP

#include "Teuchos_RCP.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"
#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  /*!
   * \brief Multivariate orthogonal polynomial basis generated from a
   * tensor product of univariate polynomials.
   */
  /*!
   * The multivariate polynomials are given by 
   * \f[
   *     \Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)
   * \f]
   * where \f$d\f$ is the dimension of the basis and \f$i_j\leq p_j\f$,
   * where \f$p_j\f$ is the order of the $j$th basis.  The size of the basis 
   * is given by \f$(p_1+1)\dots(p_d+1)\f$.
   */
  template <typename ordinal_type, typename value_type,
	    typename coeff_compare_type = 
	    TotalOrderLess<MultiIndex<ordinal_type> > >
  class TensorProductBasis : 
    public ProductBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param bases array of 1-D coordinate bases
     * \param sparse_tol tolerance used to drop terms in sparse triple-product
     *                   tensors
     */
    TensorProductBasis(
      const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,
 value_type> > >& bases,
      const value_type& sparse_tol = 1.0e-12,
      const MultiIndex<ordinal_type>& index = MultiIndex<ordinal_type>(),
      const coeff_compare_type& coeff_compare = coeff_compare_type());

    //! Destructor
    virtual ~TensorProductBasis();

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
     * represents basis polynomial \f$l\f$ and \f$i,j,k=0,\dots,P\f$ where
     * \f$P\f$ is size()-1.
     */
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor() const;

    //! Compute linear triple product tensor where k = 0,1,..,d
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeLinearTripleProductTensor() const;

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

    //! \name Implementation of Stokhos::ProductBasis methods
    //@{

    //! Get orders of each coordinate polynomial given an index \c i
    /*!
     * The returned array is of size \f$d\f$, where \f$d\f$ is the dimension of
     * the basis, and entry \f$l\f$ is given by \f$i_l\f$ where
     * \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual const MultiIndex<ordinal_type>& term(ordinal_type i) const;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type index(const MultiIndex<ordinal_type>& term) const;

    //! Return coordinate bases
    /*!
     * Array is of size dimension().
     */
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, 
							   value_type> > > 
    getCoordinateBases() const;

    //! Return maximum order allowable for each coordinate basis
    virtual MultiIndex<ordinal_type> getMaxOrders() const;

    //@}

  private:

    // Prohibit copying
    TensorProductBasis(const TensorProductBasis&);

    // Prohibit Assignment
    TensorProductBasis& operator=(const TensorProductBasis& b);
    
  protected:

    typedef MultiIndex<ordinal_type> coeff_type;
    typedef std::map<coeff_type,ordinal_type,coeff_compare_type> coeff_set_type;
    typedef Teuchos::Array<coeff_type> coeff_map_type;

    //! Name of basis
    std::string name;

    //! Total order of basis
    ordinal_type p;

    //! Total dimension of basis
    ordinal_type d;

    //! Total size of basis
    ordinal_type sz;

    //! Array of bases
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > > bases;

    //! Tolerance for computing sparse Cijk
    value_type sparse_tol;

    //! Maximum orders for each dimension
    coeff_type max_orders;

    //! Basis set
    coeff_set_type basis_set;

    //! Basis map
    coeff_map_type basis_map;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Temporary array used in basis evaluation
    mutable Teuchos::Array< Teuchos::Array<value_type> > basis_eval_tmp;

  }; // class TensorProductBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_TensorProductBasisImp.hpp"

#endif
