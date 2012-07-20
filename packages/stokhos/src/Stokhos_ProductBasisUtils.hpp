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

#ifndef STOKHOS_PRODUCT_BASIS_UTILS_HPP
#define STOKHOS_PRODUCT_BASIS_UTILS_HPP

#include "Teuchos_Array.hpp"

namespace Stokhos {

  /*! 
   * \brief Utilities for indexing a multi-variate complete polynomial basis 
   */
  template <typename ordinal_type, typename value_type>
  class CompletePolynomialBasisUtils {
  public:

    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension.
     */
    /*
     * Returns expansion total order.  
     * This version is for an isotropic expansion of total order \c p in
     * \c d dimensions.
     */
    static ordinal_type 
    compute_terms(ordinal_type p, ordinal_type d, 
		  ordinal_type& sz,
		  Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
		  Teuchos::Array<ordinal_type>& num_terms);

    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension
     */
    /*
     * Returns expansion total order.  
     * This version allows for anisotropy in the maximum order in each 
     * dimension.
     */
    static ordinal_type 
    compute_terms(const Teuchos::Array<ordinal_type>& basis_orders, 
		  ordinal_type& sz,
		  Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
		  Teuchos::Array<ordinal_type>& num_terms);

    /*!
     * \brief Compute basis index given the orders for each basis
     * dimension.
     */
    static ordinal_type 
    compute_index(const Teuchos::Array<ordinal_type>& term,
		  const Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
		  const Teuchos::Array<ordinal_type>& num_terms,
		  ordinal_type max_p);

  };

}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasisUtils<ordinal_type, value_type>::
compute_terms(ordinal_type p, ordinal_type d, 
	      ordinal_type& sz,
	      Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
	      Teuchos::Array<ordinal_type>& num_terms)
{
  Teuchos::Array<ordinal_type> basis_orders(d, p);
  return compute_terms(basis_orders, sz, terms, num_terms);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasisUtils<ordinal_type, value_type>::
compute_terms(const Teuchos::Array<ordinal_type>& basis_orders, 
	      ordinal_type& sz,
	      Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
	      Teuchos::Array<ordinal_type>& num_terms)
{
  // The approach here for ordering the terms is inductive on the total
  // order p.  We get the terms of total order p from the terms of total
  // order p-1 by incrementing the orders of the first dimension by 1.
  // We then increment the orders of the second dimension by 1 for all of the
  // terms whose first dimension order is 0.  We then repeat for the third
  // dimension whose first and second dimension orders are 0, and so on.
  // How this is done is most easily illustrated by an example of dimension 3:
  //
  // Order  terms   cnt  Order  terms   cnt
  //   0    0 0 0          4    4 0 0  15 5 1
  //                            3 1 0
  //   1    1 0 0  3 2 1        3 0 1
  //        0 1 0               2 2 0
  //        0 0 1               2 1 1
  //                            2 0 2
  //   2    2 0 0  6 3 1        1 3 0
  //        1 1 0               1 2 1
  //        1 0 1               1 1 2
  //        0 2 0               1 0 3
  //        0 1 1               0 4 0
  //        0 0 2               0 3 1
  //                            0 2 2
  //   3    3 0 0  10 4 1       0 1 3
  //        2 1 0               0 0 4
  //        2 0 1
  //        1 2 0
  //        1 1 1
  //        1 0 2
  //        0 3 0
  //        0 2 1
  //        0 1 2
  //        0 0 3

  // Compute total order
  ordinal_type d = basis_orders.size();
  ordinal_type p = 0;
  for (ordinal_type i=0; i<d; i++) {
    if (basis_orders[i] > p)
      p = basis_orders[i];
  }

  // Temporary array of terms grouped in terms of same order
  Teuchos::Array< Teuchos::Array< Teuchos::Array<ordinal_type> > > terms_order(p+1);

  // Store number of terms up to each order
  num_terms.resize(p+2, ordinal_type(0));

  // Set order 0
  terms_order[0].resize(1);
  terms_order[0][0].resize(d, ordinal_type(0));
  num_terms[0] = 1;

  // The array "cnt" stores the number of terms we need to increment for each
  // dimension.  
  Teuchos::Array<ordinal_type> cnt(d), cnt_next(d), term(d);
  for (ordinal_type j=0; j<d; j++) {
    if (basis_orders[j] >= 1)
      cnt[j] = 1;
    else
      cnt[j] = 0;
    cnt_next[j] = 0;
  }

  sz = 1;
  // Loop over orders
  for (ordinal_type k=1; k<=p; k++) {

    num_terms[k] = num_terms[k-1];

    // Stores the index of the term we copying
    ordinal_type prev = 0;

    // Loop over dimensions
    for (ordinal_type j=0; j<d; j++) {

      // Increment orders of cnt[j] terms for dimension j
      for (ordinal_type i=0; i<cnt[j]; i++) {
	if (terms_order[k-1][prev+i][j] < basis_orders[j]) {
	  term = terms_order[k-1][prev+i];
	  ++term[j];
	  terms_order[k].push_back(term);
	  ++sz;
	  num_terms[k]++;
	  for (ordinal_type l=0; l<=j; l++)
	    ++cnt_next[l];
	}
      }

      // Move forward to where all orders for dimension j are 0
      if (j < d-1)
	prev += cnt[j] - cnt[j+1];

    }

    // Update the number of terms we must increment for the new order
    for (ordinal_type j=0; j<d; j++) {
      cnt[j] = cnt_next[j];
      cnt_next[j] = 0;
    }

  }

  num_terms[p+1] = sz;

  // Copy into final terms array
  terms.resize(sz);
  ordinal_type i = 0;
  for (ordinal_type k=0; k<=p; k++) {
    ordinal_type num_k = terms_order[k].size();
    for (ordinal_type j=0; j<num_k; j++)
      terms[i++] = terms_order[k][j];
  }

  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::CompletePolynomialBasisUtils<ordinal_type, value_type>::
compute_index(const Teuchos::Array<ordinal_type>& term,
	      const Teuchos::Array< Teuchos::Array<ordinal_type> >& terms,
	      const Teuchos::Array<ordinal_type>& num_terms,
	      ordinal_type max_p)
{
  // The approach here for computing the index is to find the order block
  // corresponding to this term by adding up the component orders.  We then
  // do a linear search through the terms_order array for this order

  // First compute order of term
  ordinal_type d = term.size();
  ordinal_type ord = 0;
  for (ordinal_type i=0; i<d; i++)
    ord += term[i];
  TEUCHOS_TEST_FOR_EXCEPTION(ord < 0 || ord > max_p, std::logic_error,
		     "Stokhos::CompletePolynomialBasis::compute_index(): " <<
		     "Term has invalid order " << ord);

  // Now search through terms of that order to find a match
  ordinal_type k;
  if (ord == 0)
    k = 0;
  else
    k = num_terms[ord-1];
  ordinal_type k_max=num_terms[ord];
  bool found = false;
  while (k < k_max && !found) {
    bool found_term = true;
    for (ordinal_type j=0; j<d; j++) {
      found_term = found_term && (term[j] == terms[k][j]);
      if (!found_term)
	break;
    }
    found = found_term;
    ++k;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(k >= k_max && !found, std::logic_error,
		     "Stokhos::CompletePolynomialBasis::compute_index(): " <<
		     "Could not find specified term.");

  return k-1;
}

#endif
