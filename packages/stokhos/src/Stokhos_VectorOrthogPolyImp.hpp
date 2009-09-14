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

#include "Stokhos_ProductBasis.hpp"

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
VectorOrthogPoly(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
		 const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad) :
  basis_(basis),
  quad_(quad),
  coeff_(basis->size())
{
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
VectorOrthogPoly(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis, const typename traits_type::cloner_type& cloner,
		 const Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >& quad)
  : basis_(basis),
    quad_(quad),
    coeff_(basis->size())
{
  ordinal_type sz = basis_->size();
  for (ordinal_type i=0; i<sz; i++)
    coeff_[i] = cloner.clone(i);
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
~VectorOrthogPoly()
{
}

template <typename coeff_type>
void 
Stokhos::VectorOrthogPoly<coeff_type>::
resize(ordinal_type sz) 
{
  coeff_.resize(sz);
}

template <typename coeff_type>
void 
Stokhos::VectorOrthogPoly<coeff_type>::
reserve(ordinal_type sz) 
{
  coeff_.reserve(sz);
}

template <typename coeff_type>
typename Stokhos::VectorOrthogPoly<coeff_type>::ordinal_type
Stokhos::VectorOrthogPoly<coeff_type>::
size() const 
{
  return coeff_.size();
}

template <typename coeff_type>
Teuchos::RCP<const Stokhos::OrthogPolyBasis<typename Stokhos::VectorOrthogPoly<coeff_type>::ordinal_type, typename Stokhos::VectorOrthogPoly<coeff_type>::value_type> >
Stokhos::VectorOrthogPoly<coeff_type>::
basis() const 
{
  return basis_;
}

template <typename coeff_type>
Teuchos::RCP<const Stokhos::Quadrature<typename Stokhos::VectorOrthogPoly<coeff_type>::ordinal_type, typename Stokhos::VectorOrthogPoly<coeff_type>::value_type> >
Stokhos::VectorOrthogPoly<coeff_type>::
quadrature() const 
{
  return quad_;
}

template <typename coeff_type>
Teuchos::RCP<coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
getCoeffPtr(ordinal_type i) 
{
  return coeff_[i];
}

template <typename coeff_type>
Teuchos::RCP<const coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
getCoeffPtr(ordinal_type i) const
{
  return coeff_[i];
}

template <typename coeff_type>
void
Stokhos::VectorOrthogPoly<coeff_type>::
setCoeffPtr(ordinal_type i, const Teuchos::RCP<coeff_type>& c)
{
  coeff_[i] = c;
}

template <typename coeff_type>
coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
operator[](ordinal_type i)
{ 
  return *(coeff_[i]); 
}

template <typename coeff_type>
const coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
operator[](ordinal_type i) const
{ 
  return *(coeff_[i]); 
}

template <typename coeff_type> 
coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
term(ordinal_type dimension, ordinal_type order)
{
  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> > 
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_, true);
  ordinal_type d = basis->dimension();
  Teuchos::Array<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = product_basis->getIndex(term);
  return *(coeff_[index]);
}

template <typename coeff_type> 
const coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
term(ordinal_type dimension, ordinal_type order) const
{
  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> > 
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_, true);
  ordinal_type d = basis->dimension();
  Teuchos::Array<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = product_basis->getIndex(term);
  return *(coeff_[index]);
}

template <typename coeff_type>
void
Stokhos::VectorOrthogPoly<coeff_type>::
init(const value_type& val)
{
   ordinal_type sz = coeff_.size();
  for (ordinal_type i=0; i<sz; i++)
    traits_type::init(*(coeff_[i]), val);
}

template <typename coeff_type>
void
Stokhos::VectorOrthogPoly<coeff_type>::
evaluate(const Teuchos::Array<value_type>& basis_values, coeff_type& result) const
{
  traits_type::init(result, value_type(0));
  ordinal_type sz = coeff_.size();
  for (ordinal_type i=0; i<sz; i++)
    traits_type::update(result, basis_values[i], *(coeff_[i]));
}

template <typename coeff_type>
void
Stokhos::VectorOrthogPoly<coeff_type>::
sumIntoAllTerms(const value_type& weight,
		const Teuchos::Array<value_type>& basis_values,
		const Teuchos::Array<value_type>& basis_norms,
		const coeff_type& vec)
{
  ordinal_type sz = coeff_.size();
  for (ordinal_type i=0; i<sz; i++)
    traits_type::update(*(coeff_[i]), weight*basis_values[i]/basis_norms[i], 
			vec);
}
