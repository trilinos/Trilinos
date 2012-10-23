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
Stokhos::TensorProductBasis<ordinal_type, value_type>::
TensorProductBasis(
  const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases_,
  const value_type& sparse_tol_) :
  p(0),
  d(bases_.size()),
  sz(0),
  bases(bases_),
  basis_orders(d),
  sparse_tol(sparse_tol_),
  norms()
{

  // Compute largest order
  for (ordinal_type i=0; i<d; i++) {
    if (bases[i]->order() > p)
      p = bases[i]->order();
  }

  // Compute basis terms
  CPBUtils::compute_terms(basis_orders, sz, terms, num_terms);
    
  // Compute norms
  norms.resize(sz);
  value_type nrm;
  for (ordinal_type k=0; k<sz; k++) {
    const psop_coeff_type& term = psop->getCoeffTerm(k);
    nrm = value_type(1.0);
    for (ordinal_type i=0; i<d; i++)
      nrm = nrm * bases[i]->norm_squared(term[i]);
    norms[k] = nrm;
  }

  // Create name
  name = "Tensor product basis (";
  for (ordinal_type i=0; i<d-1; i++)
    name += bases[i]->getName() + ", ";
  name += bases[d-1]->getName() + ")";

  // Allocate array for basis evaluation
  basis_eval_tmp.resize(d);
  for (ordinal_type j=0; j<d; j++)
    basis_eval_tmp[j].resize(basis_orders[j]+1);
}

template <typename ordinal_type, typename value_type>
Stokhos::TensorProductBasis<ordinal_type, value_type>::
~TensorProductBasis()
{
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::TensorProductBasis<ordinal_type, value_type>::
order() const
{
  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::TensorProductBasis<ordinal_type, value_type>::
dimension() const
{
  return d;
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::TensorProductBasis<ordinal_type, value_type>::
size() const
{
  return sz;
}

template <typename ordinal_type, typename value_type>
const Teuchos::Array<value_type>&
Stokhos::TensorProductBasis<ordinal_type, value_type>::
norm_squared() const
{
  return norms;
}

template <typename ordinal_type, typename value_type>
const value_type&
Stokhos::TensorProductBasis<ordinal_type, value_type>::
norm_squared(ordinal_type i) const
{
  return norms[i];
}

template <typename ordinal_type, typename value_type>
Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
Stokhos::TensorProductBasis<ordinal_type, value_type>::
computeTripleProductTensor(ordinal_type order) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Fill Time");
#endif
  
}

template <typename ordinal_type, typename value_type>
value_type
Stokhos::TensorProductBasis<ordinal_type, value_type>::
evaluateZero(ordinal_type i) const
{
  // z = psi_{i_1}(0) * ... * psi_{i_d}(0) where i_1,...,i_d are the basis
  // terms for coefficient i
  const psop_coeff_type& term = psop->getCoeffTerm(i);
  value_type z = value_type(1.0);
  for (ordinal_type j=0; j<d; j++)
    z = z * bases[j]->evaluate(value_type(0.0), term[j]);

  return z;
}

template <typename ordinal_type, typename value_type>
void
Stokhos::TensorProductBasis<ordinal_type, value_type>::
evaluateBases(const Teuchos::Array<value_type>& point,
	      Teuchos::Array<value_type>& basis_vals) const
{
  for (ordinal_type j=0; j<d; j++)
    bases[j]->evaluateBases(point[j], basis_eval_tmp[j]);

  // Only evaluate basis upto number of terms included in basis_pts
  for (ordinal_type i=0; i<sz; i++) {
    const psop_coeff_type& term = psop->getCoeffTerm(i);
    value_type t = value_type(1.0);
    for (ordinal_type j=0; j<d; j++)
      t *= basis_eval_tmp[j][term[j]];
    basis_vals[i] = t;
  }
}

template <typename ordinal_type, typename value_type>
void
Stokhos::TensorProductBasis<ordinal_type, value_type>::
print(std::ostream& os) const
{
  os << "Tensor product basis of order " << p << ", dimension " << d 
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
Stokhos::TensorProductBasis<ordinal_type, value_type>::
getTerm(ordinal_type i) const
{
  return psop->getCoeffTerm(i);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::TensorProductBasis<ordinal_type, value_type>::
getIndex(const Teuchos::Array<ordinal_type>& term) const
{
  psop_coeff_type t(d);
  for (ordinal_type i=0; i<d; ++i)
    t[i] = term[i];
  return psop->getCoeffIndex(t);
}

template <typename ordinal_type, typename value_type>
const std::string&
Stokhos::TensorProductBasis<ordinal_type, value_type>::
getName() const
{
  return name;
}

template <typename ordinal_type, typename value_type>
Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > >
Stokhos::TensorProductBasis<ordinal_type, value_type>::
getCoordinateBases() const
{
  return bases;
}
