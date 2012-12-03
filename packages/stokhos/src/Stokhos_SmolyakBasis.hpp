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

#ifndef STOKHOS_SMOLYAK_BASIS_HPP
#define STOKHOS_SMOLYAK_BASIS_HPP

#include <map>

#include "Teuchos_RCP.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"
#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  /*!
   * \brief Multivariate orthogonal polynomial basis generated from a
   * Smolyak sparse grid.
   */
  template <typename ordinal_type, typename value_type,
	    typename coeff_compare_type = 
	    TotalOrderLess<MultiIndex<ordinal_type> > >
  class SmolyakBasis : 
    public ProductBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param bases array of 1-D coordinate bases
     * \param sparse_tol tolerance used to drop terms in sparse triple-product
     *                   tensors
     */
    template <typename index_set_type>
    SmolyakBasis(
      const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,
 value_type> > >& bases,
      const index_set_type& index_set,
      const value_type& sparse_tol = 1.0e-12,
      const coeff_compare_type& coeff_compare = coeff_compare_type());

    //! Destructor
    virtual ~SmolyakBasis();

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

    typedef MultiIndex<ordinal_type> coeff_type;
    typedef TensorProductBasis<ordinal_type,value_type,
			       LexographicLess<coeff_type> 
			       > tensor_product_basis_type;

    //! Return number of terms in Smolyak formula
    ordinal_type getNumSmolyakTerms() const { return tp_bases.size(); }

    //! Return ith tensor product basis
    Teuchos::RCP<const tensor_product_basis_type> 
    getTensorProductBasis(ordinal_type i) const { return tp_bases[i]; }

    //! Return ith smolyak coefficient
    ordinal_type getSmolyakCoefficient(ordinal_type i) const {
      return smolyak_coeffs[i];
    }

  private:

    // Prohibit copying
    SmolyakBasis(const SmolyakBasis&);

    // Prohibit Assignment
    SmolyakBasis& operator=(const SmolyakBasis& b);
    
  protected:

    typedef std::map<coeff_type,ordinal_type,coeff_compare_type> coeff_set_type;
    typedef Teuchos::Array<coeff_type> coeff_map_type;
    typedef MultiIndex<ordinal_type> multiindex_type;

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

    //! Smolyak coefficients
    Teuchos::Array<ordinal_type> smolyak_coeffs;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Temporary array used in basis evaluation
    mutable Teuchos::Array< Teuchos::Array<value_type> > basis_eval_tmp;

    //! Tensor product bases comprising Smolyak set
    Teuchos::Array< Teuchos::RCP< tensor_product_basis_type > > tp_bases;

    //! Predicate functor for building sparse triple products
    template <typename tp_predicate_type>
    struct SmolyakPredicate {
      Teuchos::Array<tp_predicate_type> tp_preds;

      // Predicate is true if any tensor-product predicate is true
      bool operator() (const coeff_type& term) const { 
	for (ordinal_type i=0; i<tp_preds.size(); ++i)
	  if (tp_preds[i](term))
	    return true;
	return false; 
      }

    };

    //! Predicate for building sparse triple products
    SmolyakPredicate< TensorProductPredicate<ordinal_type> > sm_pred;

  }; // class SmolyakBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_SmolyakBasisImp.hpp"

#endif
