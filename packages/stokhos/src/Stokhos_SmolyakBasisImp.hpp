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

template <typename ordinal_type, typename value_type, typename ordering_type>
template <typename index_set_type,
	  typename coeff_growth_rule_type>
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
SmolyakBasis(
  const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases_,
  const index_set_type& index_set,
  const coeff_growth_rule_type& coeff_growth_rule,
  const value_type& sparse_tol_,
  const ordering_type& coeff_compare) :
  p(0),
  d(bases_.size()),
  sz(0),
  bases(bases_),
  sparse_tol(sparse_tol_),
  max_orders(d),
  basis_set(coeff_compare),
  norms()
{
  
  // Generate index set for the final Smolyak coefficients
  //
  // The Smolyak operator is given by the formula
  //
  // A = \sum_{k\in\K} \bigotimes_{i=1}^d \Delta^i_{k_i}
  //
  // where \Delta^i_0 = 0, \Delta^i_{k_i} = L^i_{k_i} - L^i_{k_i-1},
  // and K is the supplied index set.  This becomes 
  //
  // A = \sum_{k\in\tilde{K}} c(k) \bigotimes_{i=1}^d L^i_{k_i}
  //
  // for some new index set \tilde{K} and coefficient c(k).  Using the
  // formula (cf. G W Wasilkowski and H Wozniakowski, "Explicit cost bounds 
  // of algorithms for multivariate tensor product problems," 
  // Journal of Complexity (11), 1995)
  // 
  // \bigotimes_{i=1}^d \Delta^i_{k_i} = 
  //    \sum_{\alpha\in\Alpha} (-1)^{|\alpha|} 
  //        \bigotimes_{i=1}^d L^i_{k_i-\alpha_i}
  //
  // where \Alpha = {0,1}^d and |\alpha| = \alpha_1 + ... + \alpha_d, we
  // iterate over K and \Alpha, compute k-\alpha and the corresponding
  // coefficient contribution (-1)^{|\alpha|} and store these in a map.  
  // The keys of of this map with non-zero coefficients define
  // \tilde{K} and c(k).
  typedef Stokhos::TensorProductIndexSet<ordinal_type> alpha_set_type; 
  typedef Stokhos::LexographicLess<multiindex_type> index_compare;
  typedef std::map<multiindex_type,ordinal_type,index_compare> index_map_type;
  ordinal_type dim = index_set.dimension();
  alpha_set_type alpha_set(dim, 1);
  typename alpha_set_type::iterator alpha_begin = alpha_set.begin();
  typename alpha_set_type::iterator alpha_end = alpha_set.end();
  typename index_set_type::iterator index_iterator = index_set.begin();
  typename index_set_type::iterator index_end = index_set.end();
  multiindex_type diff(dim);
  index_map_type index_map;
  for (; index_iterator != index_end; ++index_iterator) {
    for (typename alpha_set_type::iterator alpha = alpha_begin; 
	 alpha != alpha_end; ++alpha) {
      bool valid_index = true;
      for (ordinal_type i=0; i<dim; ++i) {
	diff[i] = (*index_iterator)[i] - (*alpha)[i];
	if (diff[i] < 0) {
	  valid_index = false;
	  break;
	}
      }
      if (valid_index) {
	ordinal_type alpha_order = alpha->order();
	ordinal_type val;	  
	if (alpha_order % 2 == 0)
	  val = 1;
	else
	  val = -1;
	typename index_map_type::iterator index_map_iterator =
	  index_map.find(diff);
	if (index_map_iterator == index_map.end())
	  index_map[diff] = val;
	else
	  index_map_iterator->second += val;
      }
    }
  }

  // Generate tensor product bases
  typename index_map_type::iterator index_map_iterator = index_map.begin();
  typename index_map_type::iterator index_map_end = index_map.end();
  for (; index_map_iterator != index_map_end; ++index_map_iterator) {
	
    // Skip indices with zero coefficient
    if (index_map_iterator->second == 0)
      continue;

    // Apply growth rule to cofficient multi-index
    multiindex_type coeff_growth_index(dim);
    for (ordinal_type i=0; i<dim; ++i) {
      coeff_growth_index[i] = 
	coeff_growth_rule[i](index_map_iterator->first[i]);
    }

    // Build tensor product basis for given index
    Teuchos::RCP<tensor_product_basis_type> tp = 
      Teuchos::rcp(new tensor_product_basis_type(
		     bases, sparse_tol, coeff_growth_index));

    // Include coefficients in union over all sets
    for (ordinal_type i=0; i<tp->size(); ++i)
      basis_set[tp->term(i)] = ordinal_type(0);

    tp_bases.push_back(tp);
    sm_pred.tp_preds.push_back( 
      TensorProductPredicate<ordinal_type>(coeff_growth_index) );
    smolyak_coeffs.push_back(index_map_iterator->second);
  }
  sz = basis_set.size();
 
  // Generate linear odering of coefficients
  ordinal_type idx = 0;
  basis_map.resize(sz);
  for (typename coeff_set_type::iterator i = basis_set.begin(); 
       i != basis_set.end(); 
       ++i) {
    i->second = idx;
    basis_map[idx] = i->first;
    ++idx;
  }

  // Compute max coefficient orders
  for (ordinal_type i=0; i<sz; ++i) {
    for (ordinal_type j=0; j<dim; ++j)
      if (basis_map[i][j] > max_orders[j])
	max_orders[j] = basis_map[i][j];
  }

  // Resize bases to make sure they are high enough order
  for (ordinal_type i=0; i<dim; i++)
    if (bases[i]->order() < max_orders[i])
      bases[i] = bases[i]->cloneWithOrder(max_orders[i]);

  // Compute largest order
  p = 0;
  for (ordinal_type i=0; i<d; i++) {
    if (max_orders[i] > p)
      p = max_orders[i];
  }
  
  // Compute norms
  norms.resize(sz);
  value_type nrm;
  for (ordinal_type k=0; k<sz; k++) {
    nrm = value_type(1.0);
    for (ordinal_type i=0; i<d; i++)
      nrm = nrm * bases[i]->norm_squared(basis_map[k][i]);
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
    basis_eval_tmp[j].resize(max_orders[j]+1);
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
~SmolyakBasis()
{
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
const Teuchos::Array<value_type>&
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
const value_type&
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
computeTripleProductTensor(ordinal_type order) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Fill Time");
#endif

  SmolyakPredicate< TotalOrderPredicate<ordinal_type> > k_pred;
  for (ordinal_type i=0; i<sm_pred.tp_preds.size(); ++i) {
    k_pred.tp_preds.push_back(
      TotalOrderPredicate<ordinal_type>(order, sm_pred.tp_preds[i].orders) );
  }
  
  return ProductBasisUtils::computeTripleProductTensor(
    bases, basis_set, basis_map, sm_pred, k_pred, sparse_tol);
}

template <typename ordinal_type, typename value_type, typename ordering_type>
value_type
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
evaluateZero(ordinal_type i) const
{
  // z = psi_{i_1}(0) * ... * psi_{i_d}(0) where i_1,...,i_d are the basis
  // terms for coefficient i
  value_type z = value_type(1.0);
  for (ordinal_type j=0; j<d; j++)
    z = z * bases[j]->evaluate(value_type(0.0), basis_map[i][j]);

  return z;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
void
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
evaluateBases(const Teuchos::ArrayView<const value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  for (ordinal_type j=0; j<d; j++)
    bases[j]->evaluateBases(point[j], basis_eval_tmp[j]);

  // Only evaluate basis upto number of terms included in basis_pts
  for (ordinal_type i=0; i<sz; i++) {
    value_type t = value_type(1.0);
    for (ordinal_type j=0; j<d; j++)
      t *= basis_eval_tmp[j][basis_map[i][j]];
    basis_vals[i] = t;
  }
}

template <typename ordinal_type, typename value_type, typename ordering_type>
void
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
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

template <typename ordinal_type, typename value_type, typename ordering_type>
const Stokhos::MultiIndex<ordinal_type>&
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
term(ordinal_type i) const
{
  return basis_map[i];
}

template <typename ordinal_type, typename value_type, typename ordering_type>
ordinal_type
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
index(const Stokhos::MultiIndex<ordinal_type>& term) const
{
  typename coeff_set_type::const_iterator it = basis_set.find(term);
  TEUCHOS_TEST_FOR_EXCEPTION(it == basis_set.end(), std::logic_error,
			     "Invalid term " << term);
  return it->second;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
const std::string&
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
getCoordinateBases() const
{
  return bases;
}

template <typename ordinal_type, typename value_type, typename ordering_type>
Stokhos::MultiIndex<ordinal_type>
Stokhos::SmolyakBasis<ordinal_type, value_type, ordering_type>::
getMaxOrders() const
{
  return max_orders;
}
