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

#ifndef STOKHOS_TOTAL_ORDER_BASIS_HPP
#define STOKHOS_TOTAL_ORDER_BASIS_HPP

#include "Teuchos_RCP.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_OneDOrthogPolyBasis.hpp"
#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  /*!
   * \brief Multivariate orthogonal polynomial basis generated from a
   * total order tensor product of univariate polynomials.
   */
  /*!
   * The multivariate polynomials are given by 
   * \f[
   *     \Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)
   * \f]
   * where \f$d\f$ is the dimension of the basis and \f$i_j\leq p_j\f$,
   * where \f$p_j\f$ is the order of the $j$th basis.  
   */
  template <typename ordinal_type, typename value_type,
	    typename coeff_compare_type = 
	    TotalOrderLess<MultiIndex<ordinal_type> > >
  class TotalOrderBasis : 
    public ProductBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param bases array of 1-D coordinate bases
     * \param sparse_tol tolerance used to drop terms in sparse triple-product
     *                   tensors
     */
    TotalOrderBasis(
      const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,
 value_type> > >& bases,
      const value_type& sparse_tol = 1.0e-12,
      const coeff_compare_type& coeff_compare = coeff_compare_type());

    //! Destructor
    virtual ~TotalOrderBasis();

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

  private:

    // Prohibit copying
    TotalOrderBasis(const TotalOrderBasis&);

    // Prohibit Assignment
    TotalOrderBasis& operator=(const TotalOrderBasis& b);
    
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

  }; // class TotalOrderBasis

  // An approach for building a sparse 3-tensor only for lexicographically
  // ordered total order basis
  // To-do:
  //   * Remove the n_choose_k() calls
  //   * Remove the loops in the Cijk_1D_Iterator::increment() functions
  //   * Store the 1-D Cijk tensors in a compressed format and eliminate
  //     the implicit searches with getValue()
  //     * Instead of looping over (i,j,k) multi-indices we could just store
  //       the 1-D Cijk tensors as an array of (i,j,k,c) tuples.
  template <typename ordinal_type, 
	    typename value_type>
  Teuchos::RCP< Sparse3Tensor<ordinal_type, value_type> >
  computeTripleProductTensorLTO(
    const TotalOrderBasis<ordinal_type, value_type,LexographicLess<MultiIndex<ordinal_type> > >& product_basis,
    bool symmetric = false) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
    TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Time");
#endif
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::Array;

    typedef MultiIndex<ordinal_type> coeff_type;
    const Array< RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases = product_basis.getCoordinateBases();
    ordinal_type d = bases.size();
    //ordinal_type p = product_basis.order();
    Array<ordinal_type> basis_orders(d);
    for (int i=0; i<d; ++i)
      basis_orders[i] = bases[i]->order();

    // Create 1-D triple products
    Array< RCP<Sparse3Tensor<ordinal_type,value_type> > > Cijk_1d(d);
    for (ordinal_type i=0; i<d; i++) {
      Cijk_1d[i] = 
	bases[i]->computeSparseTripleProductTensor(bases[i]->order()+1);
    }

    RCP< Sparse3Tensor<ordinal_type, value_type> > Cijk = 
      rcp(new Sparse3Tensor<ordinal_type, value_type>);
      
    // Create i, j, k iterators for each dimension
    typedef ProductBasisUtils::Cijk_1D_Iterator<ordinal_type> Cijk_Iterator;
    Array<Cijk_Iterator> Cijk_1d_iterators(d);
    coeff_type terms_i(d,0), terms_j(d,0), terms_k(d,0);
    Array<ordinal_type> sum_i(d,0), sum_j(d,0), sum_k(d,0);
    for (ordinal_type dim=0; dim<d; dim++) {
      Cijk_1d_iterators[dim] = Cijk_Iterator(bases[dim]->order(), symmetric);
    }

    ordinal_type I = 0;
    ordinal_type J = 0;
    ordinal_type K = 0;
    ordinal_type cnt = 0;
    bool stop = false;
    while (!stop) {

      // Fill out terms from 1-D iterators
      for (ordinal_type dim=0; dim<d; ++dim) {
	terms_i[dim] = Cijk_1d_iterators[dim].i;
	terms_j[dim] = Cijk_1d_iterators[dim].j;
	terms_k[dim] = Cijk_1d_iterators[dim].k;
      }

      // Compute global I,J,K
      /*
      ordinal_type II = lexicographicMapping(terms_i, p);
      ordinal_type JJ = lexicographicMapping(terms_j, p);
      ordinal_type KK = lexicographicMapping(terms_k, p);
      if (I != II || J != JJ || K != KK) {
	std::cout << "DIFF!!!" << std::endl;
	std::cout << terms_i << ":  I = " << I << ", II = " << II << std::endl;
	std::cout << terms_j << ":  J = " << J << ", JJ = " << JJ << std::endl;
	std::cout << terms_k << ":  K = " << K << ", KK = " << KK << std::endl;
      }
      */

      // Compute triple-product value
      value_type c = value_type(1.0);
      for (ordinal_type dim=0; dim<d; dim++) {
	c *= Cijk_1d[dim]->getValue(Cijk_1d_iterators[dim].i, 
				    Cijk_1d_iterators[dim].j, 
				    Cijk_1d_iterators[dim].k);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
	std::abs(c) <= 1.0e-12,
	std::logic_error,
	"Got 0 triple product value " << c 
	<< ", I = " << I << " = " << terms_i
	<< ", J = " << J << " = " << terms_j
	<< ", K = " << K << " = " << terms_k
	<< std::endl);
      
      // Add term to global Cijk
      Cijk->add_term(I,J,K,c);
      // Cijk->add_term(I,K,J,c);
      // Cijk->add_term(J,I,K,c);
      // Cijk->add_term(J,K,I,c);
      // Cijk->add_term(K,I,J,c);
      // Cijk->add_term(K,J,I,c);
      
      // Increment iterators to the next valid term
      ordinal_type cdim = d-1;
      bool inc = true;
      while (inc && cdim >= 0) {
	ordinal_type delta_i, delta_j, delta_k;
	bool more = 
	  Cijk_1d_iterators[cdim].increment(delta_i, delta_j, delta_k);

	// Update number of terms used for computing global index
	if (cdim == d-1) {
	  I += delta_i;
	  J += delta_j;
	  K += delta_k;
	}
	else {
	  if (delta_i > 0) {
	    for (ordinal_type ii=0; ii<delta_i; ++ii)
	      I += 
		Stokhos::n_choose_k(
		  basis_orders[cdim+1]-sum_i[cdim] -
		  (Cijk_1d_iterators[cdim].i-ii)+d-cdim,
		  d-cdim-1);
	  }
	  else {
	    for (ordinal_type ii=0; ii<-delta_i; ++ii)
	      I -= 
		Stokhos::n_choose_k(
		  basis_orders[cdim+1]-sum_i[cdim] -
		  (Cijk_1d_iterators[cdim].i+ii)+d-cdim-1,
		  d-cdim-1);
	  }

	  if (delta_j > 0) {
	    for (ordinal_type jj=0; jj<delta_j; ++jj)
	      J +=  
		Stokhos::n_choose_k(
		  basis_orders[cdim+1]-sum_j[cdim] -
		  (Cijk_1d_iterators[cdim].j-jj)+d-cdim,
		  d-cdim-1);
	  }
	  else {
	    for (ordinal_type jj=0; jj<-delta_j; ++jj)
	      J -=  
		Stokhos::n_choose_k(
		  basis_orders[cdim+1]-sum_j[cdim] -
		  (Cijk_1d_iterators[cdim].j+jj)+d-cdim-1,
		  d-cdim-1);
	  }
	  
	  if (delta_k > 0) {
	    for (ordinal_type kk=0; kk<delta_k; ++kk)
	      K +=  
		Stokhos::n_choose_k(
		  basis_orders[cdim+1]-sum_k[cdim] -
		  (Cijk_1d_iterators[cdim].k-kk)+d-cdim,
		  d-cdim-1);
	  }
	  else {
	    for (ordinal_type kk=0; kk<-delta_k; ++kk)
	      K -=  
		Stokhos::n_choose_k(
		  basis_orders[cdim+1]-sum_k[cdim] -
		  (Cijk_1d_iterators[cdim].k+kk)+d-cdim-1,
		  d-cdim-1);
	  }
	}
	
	if (!more) {
	  // If no more terms in this dimension, go to previous one
	  --cdim;
	}
	else {
	  // cdim has more terms, so reset iterators for all dimensions > cdim
	  // adjusting max order based on sum of i,j,k for previous dims
	  inc = false;

	  for (ordinal_type dim=cdim+1; dim<d; ++dim) {

	    // Update sums of orders for previous dimension
	    sum_i[dim] = sum_i[dim-1] + Cijk_1d_iterators[dim-1].i;
	    sum_j[dim] = sum_j[dim-1] + Cijk_1d_iterators[dim-1].j;
	    sum_k[dim] = sum_k[dim-1] + Cijk_1d_iterators[dim-1].k;

	    // Reset iterator for this dimension
	    Cijk_1d_iterators[dim] = 
	      Cijk_Iterator(basis_orders[dim]-sum_i[dim],
			    basis_orders[dim]-sum_j[dim],
			    basis_orders[dim]-sum_k[dim],
			    symmetric);
	  }
	}
      }
      
      if (cdim < 0)
	stop = true;
      
      cnt++;
    }
    
    Cijk->fillComplete();
    
    return Cijk;
  }

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_TotalOrderBasisImp.hpp"

#endif
