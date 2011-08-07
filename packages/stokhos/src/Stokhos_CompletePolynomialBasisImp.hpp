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
#include "Teuchos_TimeMonitor.hpp"

template <typename ordinal_type, typename value_type>
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
CompletePolynomialBasis(
	const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases_,
	const value_type& sparse_tol_,
	bool use_old_cijk_alg_,
	const Teuchos::RCP< Teuchos::Array<value_type> >& deriv_coeffs_) :
  p(0),
  d(bases_.size()),
  sz(0),
  bases(bases_),
  basis_orders(d),
  sparse_tol(sparse_tol_),
  use_old_cijk_alg(use_old_cijk_alg_),
  deriv_coeffs(deriv_coeffs_),
  norms(),
  terms()
{

  // Compute total order
  for (ordinal_type i=0; i<d; i++) {
    basis_orders[i] = bases[i]->order();
    if (basis_orders[i] > p)
      p = basis_orders[i];
  }

  // Compute basis terms
  compute_terms();
    
  // Compute norms
  norms.resize(sz);
  value_type nrm;
  for (ordinal_type k=0; k<sz; k++) {
    nrm = value_type(1.0);
    for (ordinal_type i=0; i<d; i++)
      nrm = nrm * bases[i]->norm_squared(terms[k][i]);
    norms[k] = nrm;
  }

  // Create name
  name = "Complete polynomial basis (";
  for (ordinal_type i=0; i<d-1; i++)
    name += bases[i]->getName() + ", ";
  name += bases[d-1]->getName() + ")";

  // Allocate array for basis evaluation
  basis_eval_tmp.resize(d);
  for (ordinal_type j=0; j<d; j++)
    basis_eval_tmp[j].resize(basis_orders[j]+1);

  // Set up deriv_coeffs
  if (deriv_coeffs == Teuchos::null) {
    deriv_coeffs = Teuchos::rcp(new Teuchos::Array<value_type>(d));
    for (ordinal_type j=0; j<d; j++)
      (*deriv_coeffs)[j] = value_type(1.0);
  }
}

template <typename ordinal_type, typename value_type>
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
~CompletePolynomialBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
computeTripleProductTensor(ordinal_type order) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Fill Time");
#endif
  if (use_old_cijk_alg)
    return computeTripleProductTensorOld(order);
  else
    return computeTripleProductTensorNew(order);
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
computeTripleProductTensorOld(ordinal_type order) const
{
  // Compute Cijk = < \Psi_i \Psi_j \Psi_k >
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>);

  // Create 1-D triple products
  Teuchos::Array< Teuchos::RCP<Dense3Tensor<ordinal_type,value_type> > > Cijk_1d(d);
  for (ordinal_type i=0; i<d; i++)
    Cijk_1d[i] = bases[i]->computeTripleProductTensor();

  for (ordinal_type j=0; j<sz; j++) {
    for (ordinal_type i=0; i<sz; i++) {
      for (ordinal_type k=0; k<order; k++) {
	value_type c = value_type(1.0);
	for (ordinal_type l=0; l<d; l++)
	  c *= (*Cijk_1d[l])(terms[i][l],terms[j][l],terms[k][l]);
	if (std::abs(c/norms[i]) > sparse_tol)
	  Cijk->add_term(i,j,k,c);
      }
    }
  }

  Cijk->fillComplete();

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
computeTripleProductTensorNew(ordinal_type order) const
{
  // The algorithm for computing Cijk = < \Psi_i \Psi_j \Psi_k > here works
  // by factoring 
  // < \Psi_i \Psi_j \Psi_k > = 
  //    < \psi^1_{i_1}\psi^1_{j_1}\psi^1_{k_1} >_1 x ... x
  //    < \psi^d_{i_d}\psi^d_{j_d}\psi^d_{k_d} >_d
  // We compute the sparse triple product < \psi^l_i\psi^l_j\psi^l_k >_l
  // for each dimension l, and then compute all non-zero products of these
  // terms.  The complexity arises from iterating through all possible
  // combinations, throwing out terms that aren't in the basis and are beyond
  // the k-order limit provided
  Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
    Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>);

  // Map the specified order limit to a limit on each dimension
  // Subtract 1 to get the term for the last order we want to include,
  // add up the orders for each term to get the total order, then add 1
  Teuchos::Array<ordinal_type> term = getTerm(order-1);
  ordinal_type k_lim = 0;
  for (ordinal_type i=0; i<d; i++)
    k_lim = k_lim + term[i];
  k_lim++;

  // Create 1-D triple products
  Teuchos::Array< Teuchos::RCP<Sparse3Tensor<ordinal_type,value_type> > > Cijk_1d(d);
  for (ordinal_type i=0; i<d; i++) {
    if (k_lim <= basis_orders[i]+1)
      Cijk_1d[i] = bases[i]->computeSparseTripleProductTensor(k_lim);
    else
      Cijk_1d[i] = bases[i]->computeSparseTripleProductTensor(basis_orders[i]+1);
  }

  // Create i, j, k iterators for each dimension
  // Note:  we have to supply an initializer in the arrays of iterators to 
  // avoid checked-stl errors about singular iterators
  typedef Sparse3Tensor<ordinal_type,value_type> Cijk_type;
  typedef typename Cijk_type::k_iterator k_iterator;
  typedef typename Cijk_type::kj_iterator kj_iterator;
  typedef typename Cijk_type::kji_iterator kji_iterator;
  Teuchos::Array<k_iterator> k_iterators(d, Cijk_1d[0]->k_begin());
  Teuchos::Array<kj_iterator > j_iterators(d, Cijk_1d[0]->j_begin(k_iterators[0]));
  Teuchos::Array<kji_iterator > i_iterators(d, Cijk_1d[0]->i_begin(j_iterators[0]));
  Teuchos::Array<ordinal_type> terms_i(d), terms_j(d), terms_k(d);
  ordinal_type sum_i = 0;
  ordinal_type sum_j = 0;
  ordinal_type sum_k = 0;
  for (ordinal_type dim=0; dim<d; dim++) {
    k_iterators[dim] = Cijk_1d[dim]->k_begin();
    j_iterators[dim] = Cijk_1d[dim]->j_begin(k_iterators[dim]);
    i_iterators[dim] = Cijk_1d[dim]->i_begin(j_iterators[dim]);
    terms_i[dim] = index(i_iterators[dim]);
    terms_j[dim] = index(j_iterators[dim]);
    terms_k[dim] = index(k_iterators[dim]);
    sum_i += terms_i[dim];
    sum_j += terms_j[dim];
    sum_k += terms_k[dim];
  }

  ordinal_type I = 0;
  ordinal_type J = 0;
  ordinal_type K = 0;
  bool inc_i = true;
  bool inc_j = true;
  bool inc_k = true;
  bool stop = false;
  ordinal_type cnt = 0;
  while (!stop) {

    // Add term if it is in the basis
    if (sum_i <= p && sum_j <= p && sum_k <= p) {
      if (inc_k)
	K = compute_index(terms_k);
      if (K < order) {
	if (inc_i)
	  I = compute_index(terms_i);
	if (inc_j)
	  J = compute_index(terms_j);
	value_type c = value_type(1.0);
	for (ordinal_type dim=0; dim<d; dim++)
	  c *= value(i_iterators[dim]);
	if (std::abs(c/norms[I]) > sparse_tol)
	  Cijk->add_term(I,J,K,c);
      }
    }
    
    // Increment iterators to the next valid term
    ordinal_type cdim = 0;
    bool inc = true;
    inc_i = false;
    inc_j = false; 
    inc_k = false;
    while (inc && cdim < d) {
      ordinal_type cur_dim = cdim;
      ++i_iterators[cdim];
      inc_i = true;
      if (i_iterators[cdim] != Cijk_1d[cdim]->i_end(j_iterators[cdim])) {
	terms_i[cdim] = index(i_iterators[cdim]);
	sum_i = 0;
	for (int dim=0; dim<d; dim++)
	  sum_i += terms_i[dim];
      }
      if (i_iterators[cdim] == Cijk_1d[cdim]->i_end(j_iterators[cdim]) ||
	sum_i > p) {
	++j_iterators[cdim];
	inc_j = true;
	if (j_iterators[cdim] != Cijk_1d[cdim]->j_end(k_iterators[cdim])) {
	  terms_j[cdim] = index(j_iterators[cdim]);
	  sum_j = 0;
	  for (int dim=0; dim<d; dim++)
	    sum_j += terms_j[dim];
	}
	if (j_iterators[cdim] == Cijk_1d[cdim]->j_end(k_iterators[cdim]) ||
	  sum_j > p) {
	  ++k_iterators[cdim];
	  inc_k = true;
	  if (k_iterators[cdim] != Cijk_1d[cdim]->k_end()) {
	    terms_k[cdim] = index(k_iterators[cdim]);
	    sum_k = 0;
	    for (int dim=0; dim<d; dim++)
	      sum_k += terms_k[dim];
	  }
	  if (k_iterators[cdim] == Cijk_1d[cdim]->k_end() || sum_k > p) {
	    k_iterators[cdim] = Cijk_1d[cdim]->k_begin();
	    ++cdim;
	    terms_k[cur_dim] = index(k_iterators[cur_dim]);
	    sum_k = 0;
	    for (int dim=0; dim<d; dim++)
	      sum_k += terms_k[dim];
	  }
	  else
	    inc = false;
	  j_iterators[cur_dim] = 
	    Cijk_1d[cur_dim]->j_begin(k_iterators[cur_dim]);
	  terms_j[cur_dim] = index(j_iterators[cur_dim]);
	  sum_j = 0;
	  for (int dim=0; dim<d; dim++)
	    sum_j += terms_j[dim];
	}
	else
	  inc = false;
	i_iterators[cur_dim] = Cijk_1d[cur_dim]->i_begin(j_iterators[cur_dim]);
	terms_i[cur_dim] = index(i_iterators[cur_dim]);
	sum_i = 0;
	for (int dim=0; dim<d; dim++)
	  sum_i += terms_i[dim];
      }
      else
	inc = false;

      if (sum_i > p || sum_j > p || sum_k > p)
	inc = true;
    }

    if (cdim == d)
      stop = true;

    cnt++;
  }

  Cijk->fillComplete();

  return Cijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
computeDerivTripleProductTensor(
  const Teuchos::RCP< const Teuchos::SerialDenseMatrix<ordinal_type, value_type> >& Bij,
  const Teuchos::RCP< const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk) const
{
  // Compute Dijk = < \Psi_i \Psi_j \Psi_k' >
  Teuchos::RCP< Stokhos::Dense3Tensor<ordinal_type, value_type> > Dijk = 
    Teuchos::rcp(new Dense3Tensor<ordinal_type, value_type>(sz));
  for (ordinal_type i=0; i<sz; i++)
    for (ordinal_type j=0; j<sz; j++)
      for (ordinal_type k=0; k<sz; k++)
	(*Dijk)(i,j,k) = value_type(0.0);

  ordinal_type i,j,m;
  value_type c;
  for (ordinal_type k=0; k<sz; k++) {
    for (typename Cijk_type::k_iterator m_it=Cijk->k_begin(); 
	 m_it!=Cijk->k_end(); ++m_it) {
      m = index(m_it);
      for (typename Cijk_type::kj_iterator j_it = Cijk->j_begin(m_it); 
	   j_it != Cijk->j_end(m_it); ++j_it) {
	j = index(j_it);
	for (typename Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	     i_it != Cijk->i_end(j_it); ++i_it) {
	  i = index(i_it);
	  c = value(i_it);
	  (*Dijk)(i,j,k) += (*Bij)(m,k)*c/norms[m];
	}
      }
    }
  }

  return Dijk;
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
computeDerivDoubleProductTensor() const
{
  // Compute Bij = < \Psi_i \Psi_j' >
  Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > Bij = 
    Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type,value_type>(sz,
									 sz));
  
  // Create products
  Teuchos::Array< Teuchos::RCP<const Teuchos::SerialDenseMatrix<ordinal_type,value_type> > > Bij_1d(d);
  for (ordinal_type i=0; i<d; i++)
    Bij_1d[i] = bases[i]->computeDerivDoubleProductTensor();
  
  for (ordinal_type i=0; i<sz; i++) {
    for (ordinal_type k=0; k<sz; k++) {
      value_type t = value_type(1.0);
      value_type c = value_type(0.0);
      for (ordinal_type j=0; j<d; j++) {
	bool is_zero = false;
	for (ordinal_type l=0; l<d; l++) {
	  if (l != j && terms[i][l] != terms[k][l])
	    is_zero = true;
	  if (l != j)
	    t *= bases[l]->norm_squared(terms[k][l]);
	}
	if (!is_zero)
	  c += t*(*deriv_coeffs)[j]*(*Bij_1d[j])(terms[k][j],terms[i][j]);
      }
      (*Bij)(i,k) = c;
    }
  }

  return Bij;
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  // z = psi_{i_1}(0) * ... * psi_{i_d}(0) where i_1,...,i_d are the basis
  // terms for coefficient i
  value_type z = value_type(1.0);
  for (ordinal_type j=0; j<d; j++)
    z = z * bases[j]->evaluate(value_type(0.0), terms[i][j]);

  return z;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::Array<value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  for (ordinal_type j=0; j<d; j++)
    bases[j]->evaluateBases(point[j], basis_eval_tmp[j]);

  // Only evaluate basis upto number of terms included in basis_pts
  for (ordinal_type i=0; i<sz; i++) {
    value_type t = value_type(1.0);
    for (ordinal_type j=0; j<d; j++)
      t *= basis_eval_tmp[j][terms[i][j]];
    basis_vals[i] = t;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Complete basis of order " << p << ", dimension " << d 
     << ", and size " << sz << ".  Component bases:\n";
  for (ordinal_type i=0; i<d; i++)
    os << *bases[i];
  os << "Basis vector norms (squared):\n\t";
  for (ordinal_type i=0; i<static_cast<ordinal_type>(norms.size()); i++)
    os << norms[i] << " ";
  os << "\n";
}

template <typename ordinal_type, typename value_type>
Teuchos::Array<ordinal_type>
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  return terms[i];
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getIndex(const Teuchos::Array<ordinal_type>& term) const
{
  return compute_index(term);
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return bases;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
compute_terms()
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

  /*
  std::cout << "sz = " << sz << std::endl;
  for (ordinal_type i=0; i<sz; i++) {
    std::cout << i << ":  ";
    for (ordinal_type j=0; j<d; j++)
      std::cout << terms[i][j] << " ";
    std::cout << std::endl;
  }
  std::cout << "num_terms = ";
  for (ordinal_type i=0; i<=p; i++)
    std::cout << num_terms[i] << " ";
  std::cout << std::endl;
  */
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::CompletePolynomialBasis<ordinal_type, value_type>::
compute_index(const Teuchos::Array<ordinal_type>& term) const
{
  // The approach here for computing the index is to find the order block
  // corresponding to this term by adding up the component orders.  We then
  // do a linear search through the terms_order array for this order

  // First compute order of term
  ordinal_type ord = 0;
  for (ordinal_type i=0; i<d; i++)
    ord += term[i];
  TEST_FOR_EXCEPTION(ord < 0 || ord > p, std::logic_error,
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
  TEST_FOR_EXCEPTION(k >= k_max && !found, std::logic_error,
		     "Stokhos::CompletePolynomialBasis::compute_index(): " <<
		     "Could not find specified term.");

  return k-1;
}
