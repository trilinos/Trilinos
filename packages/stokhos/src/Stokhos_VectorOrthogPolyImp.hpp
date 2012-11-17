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
VectorOrthogPoly() : 
  ProductContainer<coeff_type>(),
  basis_()
{
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
VectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis,
  const Teuchos::RCP<const Epetra_BlockMap>& map) :
  ProductContainer<coeff_type>(map),
  basis_(basis)
{
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
VectorOrthogPoly(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis, 
  const Teuchos::RCP<const Epetra_BlockMap>& map,
  const typename traits_type::cloner_type& cloner)
  : ProductContainer<coeff_type>(map, cloner),
    basis_(basis)
{
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
VectorOrthogPoly(const Stokhos::VectorOrthogPoly<coeff_type>& v) : 
  ProductContainer<coeff_type>(v),
  basis_(v.basis_)
{
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
~VectorOrthogPoly()
{
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>&
Stokhos::VectorOrthogPoly<coeff_type>::
operator=(const Stokhos::VectorOrthogPoly<coeff_type>& v)
{
  if (this != &v) {
    ProductContainer<coeff_type>::operator=(v);
    basis_ = v.basis_;
  }
  return *this;
}

template <typename coeff_type>
void 
Stokhos::VectorOrthogPoly<coeff_type>::
reset(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& new_basis,
  const Teuchos::RCP<const Epetra_BlockMap>& new_map,
  const typename traits_type::cloner_type& cloner)
{
  basis_ = new_basis;
  ProductContainer<coeff_type>::reset(new_map, cloner);
}

template <typename coeff_type>
Teuchos::RCP<const Stokhos::OrthogPolyBasis<typename Stokhos::VectorOrthogPoly<coeff_type>::ordinal_type, typename Stokhos::VectorOrthogPoly<coeff_type>::value_type> >
Stokhos::VectorOrthogPoly<coeff_type>::
basis() const 
{
  return basis_;
}


template <typename coeff_type> 
coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
term(ordinal_type dimension, ordinal_type order)
{
  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> > 
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_, true);
  ordinal_type d = product_basis->dimension();
  MultiIndex<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = product_basis->index(term);
  return *(this->coeff_[this->map_->LID(index)]);
}

template <typename coeff_type> 
const coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
term(ordinal_type dimension, ordinal_type order) const
{
  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> > 
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_, true);
  ordinal_type d = product_basis->dimension();
  MultiIndex<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = product_basis->index(term);
  return *(this->coeff_[this->map_->LID(index)]);
}

template <typename coeff_type>
void
Stokhos::VectorOrthogPoly<coeff_type>::
evaluate(const Teuchos::Array<value_type>& basis_values, coeff_type& result) const
{
  traits_type::init(result, value_type(0));
  ordinal_type sz = this->coeff_.size();
  for (ordinal_type i=0; i<sz; i++)
    traits_type::update(result, basis_values[i], *(this->coeff_[i]));
}

template <typename coeff_type>
void
Stokhos::VectorOrthogPoly<coeff_type>::
sumIntoAllTerms(const value_type& weight,
		const Teuchos::Array<value_type>& basis_values,
		const Teuchos::Array<value_type>& basis_norms,
		const coeff_type& vec)
{
  ordinal_type sz = this->coeff_.size();
  int i_gid;
  for (ordinal_type i=0; i<sz; i++) {
    i_gid = this->map_->GID(i);
    traits_type::update(*(this->coeff_[i]), 
			weight*basis_values[i_gid]/basis_norms[i_gid], 
			vec);
  }
}

template <typename coeff_type>
std::ostream&
Stokhos::VectorOrthogPoly<coeff_type>::
print(std::ostream& os) const
{
  ordinal_type sz = this->coeff_.size();
  os << "Stokhos::VectorOrthogPoly of global size " 
     << this->map_->NumGlobalElements() << ", local size " << sz << " in basis "
     << "\n" << basis_->getName() << ":" << std::endl;

  Teuchos::RCP< const Stokhos::ProductBasis<ordinal_type, value_type> > 
    product_basis = Teuchos::rcp_dynamic_cast< const Stokhos::ProductBasis<ordinal_type, value_type> >(basis_);

  if (product_basis != Teuchos::null) {
    for (ordinal_type i=0; i<sz; i++) {
      const MultiIndex<ordinal_type>& trm = product_basis->term(i);
      os << "Term " << i << " (";
      for (ordinal_type j=0; j<static_cast<ordinal_type>(trm.size())-1; j++)
	os << trm[j] << ", ";
      os << trm[trm.size()-1] << "):" << std::endl;
      traits_type::print(os, *(this->coeff_[this->map_->LID(i)]));
    }
  }
  else {
    for (ordinal_type i=0; i<sz; i++) {
      os << "Term " << this->map_->GID(i) << ":" << std::endl;
      traits_type::print(os, *(this->coeff_[i]));
    }
  }

  return os;
}
