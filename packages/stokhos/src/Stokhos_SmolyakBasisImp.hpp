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
Stokhos::SmolyakBasis<ordinal_type, value_type>::
SmolyakBasis(
	const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases_,
	const value_type& sparse_tol_) :
  p(0),
  d(bases_.size()),
  sz(0),
  bases(bases_),
  basis_orders(d),
  sparse_tol(sparse_tol_),
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
  compute_terms(basis_orders, sz, terms, num_terms);
    
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
  name = "Smolyak basis (";
  for (ordinal_type i=0; i<d-1; i++)
    name += bases[i]->getName() + ", ";
  name += bases[d-1]->getName() + ")";

  // Allocate array for basis evaluation
  basis_eval_tmp.resize(d);
  for (ordinal_type j=0; j<d; j++)
    basis_eval_tmp[j].resize(basis_orders[j]+1);
}

template <typename ordinal_type, typename value_type>
Stokhos::SmolyakBasis<ordinal_type, value_type>::
~SmolyakBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::SmolyakBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::SmolyakBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::SmolyakBasis<ordinal_type, value_type>::
computeTripleProductTensor(ordinal_type order) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Fill Time");
#endif
  return computeTripleProductTensorNew(order);
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::SmolyakBasis<ordinal_type, value_type>::
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
	K = compute_index(terms_k, terms, num_terms, p);
      if (K < order) {
	if (inc_i)
	  I = compute_index(terms_i, terms, num_terms, p);
	if (inc_j)
	  J = compute_index(terms_j, terms, num_terms, p);
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
value_type
Stokhos::SmolyakBasis<ordinal_type, value_type>::
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
Stokhos::SmolyakBasis<ordinal_type, value_type>::
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
Stokhos::SmolyakBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Smolyak basis of order " << p << ", dimension " << d 
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
Stokhos::SmolyakBasis<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  return terms[i];
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type>::
getIndex(const Teuchos::Array<ordinal_type>& term) const
{
  // Currentl do a linear search.  Will have to find a better
  // data structure.
  bool found = false;
  ordinal_type k=0;
  while (!found) {
    boo matches = true;
    for (ordinal_type i=0; i<d; i++) {
      matches = matches && term[i] == terms[k][i];
      if (!matches)
	break;
    }
    if (matches)
      found = true;
    else
      k++;
  }

  return k;
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::SmolyakBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::SmolyakBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return bases;
}
