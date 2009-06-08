// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
VectorOrthogPoly(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis) :
  basis_(basis),
  coeff_(basis->size())
{
}

template <typename coeff_type>
Stokhos::VectorOrthogPoly<coeff_type>::
VectorOrthogPoly(const Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >& basis, const typename Stokhos::VectorOrthogPolyTraits<coeff_type>::cloner_type& cloner)
  : basis_(basis),
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
term(int i0, int i1, int i2, int i3, int i4,
     int i5, int i6, int i7, int i8, int i9)
{
  std::vector<ordinal_type> trm;
  ordinal_type d = basis.dimension();
  if (i0 >= 0 && d >= 1)
    trm.push_back(i0);
  if (i1 >= 0 && d >= 2)
    trm.push_back(i1);
  if (i2 >= 0 && d >= 3)
    trm.push_back(i2);
  if (i3 >= 0 && d >= 4)
    trm.push_back(i3);
  if (i4 >= 0 && d >= 5)
    trm.push_back(i4);
  if (i5 >= 0 && d >= 6)
    trm.push_back(i5);
  if (i6 >= 0 && d >= 7)
    trm.push_back(i6);
  if (i7 >= 0 && d >= 8)
    trm.push_back(i7);
  if (i8 >= 0 && d >= 9)
    trm.push_back(i8);
  if (i9 >= 0 && d >= 10)
    trm.push_back(i9);

  ordinal_type index = basis.getIndex(trm);
  return *(coeff_[index]);
}

template <typename coeff_type> 
const coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
term(int i0, int i1, int i2, int i3, int i4,
     int i5, int i6, int i7, int i8, int i9) const
{
  std::vector<ordinal_type> trm;
  ordinal_type d = basis.dimension();
  if (i0 >= 0 && d >= 1)
    trm.push_back(i0);
  if (i1 >= 0 && d >= 2)
    trm.push_back(i1);
  if (i2 >= 0 && d >= 3)
    trm.push_back(i2);
  if (i3 >= 0 && d >= 4)
    trm.push_back(i3);
  if (i4 >= 0 && d >= 5)
    trm.push_back(i4);
  if (i5 >= 0 && d >= 6)
    trm.push_back(i5);
  if (i6 >= 0 && d >= 7)
    trm.push_back(i6);
  if (i7 >= 0 && d >= 8)
    trm.push_back(i7);
  if (i8 >= 0 && d >= 9)
    trm.push_back(i8);
  if (i9 >= 0 && d >= 10)
    trm.push_back(i9);

  ordinal_type index = basis.getIndex(trm);
  return *(coeff_[index]);
}

template <typename coeff_type> 
coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
term2(ordinal_type dimension, ordinal_type order)
{
  ordinal_type d = basis.dimension();
  std::vector<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = basis.getIndex(term);
  return *(coeff_[index]);
}

template <typename coeff_type> 
const coeff_type&
Stokhos::VectorOrthogPoly<coeff_type>::
term2(ordinal_type dimension, ordinal_type order) const
{
  ordinal_type d = basis.dimension();
  std::vector<ordinal_type> term(d);
  term[dimension] = order;
  ordinal_type index = basis.getIndex(term);
  return *(coeff_[index]);
}
