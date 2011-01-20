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

template <typename coeff_type>
Stokhos::ProductContainer<coeff_type>::
ProductContainer() :
  map_(),
  coeff_()
{
}

template <typename coeff_type>
Stokhos::ProductContainer<coeff_type>::
ProductContainer(const Teuchos::RCP<const Epetra_BlockMap>& map) :
  map_(map),
  coeff_(map_->NumMyElements())
{
}

template <typename coeff_type>
Stokhos::ProductContainer<coeff_type>::
ProductContainer(const Teuchos::RCP<const Epetra_BlockMap>& map,
		 const typename traits_type::cloner_type& cloner) : 
  map_(map),
  coeff_(map_->NumMyElements())
{
  ordinal_type sz = map_->NumMyElements();
  for (ordinal_type i=0; i<sz; i++)
    coeff_[i] = cloner.clone(i);
}

template <typename coeff_type>
Stokhos::ProductContainer<coeff_type>::
ProductContainer(const Stokhos::ProductContainer<coeff_type>& v) : 
  map_(v.map_),
  coeff_(v.coeff_)
{
}

template <typename coeff_type>
Stokhos::ProductContainer<coeff_type>::
~ProductContainer()
{
}

template <typename coeff_type>
Stokhos::ProductContainer<coeff_type>&
Stokhos::ProductContainer<coeff_type>::
operator=(const Stokhos::ProductContainer<coeff_type>& v)
{
  if (this != &v) {
    map_ = v.map_;
    coeff_ = v.coeff_;
  }
  return *this;
}

template <typename coeff_type>
void 
Stokhos::ProductContainer<coeff_type>::
reset(const Teuchos::RCP<const Epetra_BlockMap>& map)
{
  map_ = map;
  ordinal_type sz = map_->NumMyElements();
  coeff_.resize(sz);
}

template <typename coeff_type>
void 
Stokhos::ProductContainer<coeff_type>::
reset(const Teuchos::RCP<const Epetra_BlockMap>& map,
      const typename traits_type::cloner_type& cloner)
{
  map_ = map;
  ordinal_type sz = map_->NumMyElements();
  coeff_.resize(sz);
  for (ordinal_type i=0; i<sz; i++)
    coeff_[i] = cloner.clone(i);
}

template <typename coeff_type>
void 
Stokhos::ProductContainer<coeff_type>::
resize(const Teuchos::RCP<const Epetra_BlockMap>& map) 
{
  map_ = map;
  coeff_.resize(map_->NumMyElements());
}

template <typename coeff_type>
void 
Stokhos::ProductContainer<coeff_type>::
reserve(ordinal_type sz) 
{
  coeff_.reserve(sz);
}

template <typename coeff_type>
typename Stokhos::ProductContainer<coeff_type>::ordinal_type
Stokhos::ProductContainer<coeff_type>::
size() const 
{
  return coeff_.size();
}

template <typename coeff_type>
Teuchos::RCP<const Epetra_BlockMap>
Stokhos::ProductContainer<coeff_type>::
map() const
{
  return map_;
}

template <typename coeff_type>
const Teuchos::Array<Teuchos::RCP<coeff_type> >&
Stokhos::ProductContainer<coeff_type>::
getCoefficients() const
{
  return coeff_;
}

template <typename coeff_type>
Teuchos::Array<Teuchos::RCP<coeff_type> >&
Stokhos::ProductContainer<coeff_type>::
getCoefficients()
{
  return coeff_;
}

template <typename coeff_type>
Teuchos::RCP<coeff_type>
Stokhos::ProductContainer<coeff_type>::
getCoeffPtr(ordinal_type i) 
{
  return coeff_[i];
}

template <typename coeff_type>
Teuchos::RCP<const coeff_type>
Stokhos::ProductContainer<coeff_type>::
getCoeffPtr(ordinal_type i) const
{
  return coeff_[i];
}

template <typename coeff_type>
void
Stokhos::ProductContainer<coeff_type>::
setCoeffPtr(ordinal_type i, const Teuchos::RCP<coeff_type>& c)
{
  coeff_[i] = c;
}

template <typename coeff_type>
coeff_type&
Stokhos::ProductContainer<coeff_type>::
operator[](ordinal_type i)
{ 
  return *(coeff_[i]); 
}

template <typename coeff_type>
const coeff_type&
Stokhos::ProductContainer<coeff_type>::
operator[](ordinal_type i) const
{ 
  return *(coeff_[i]); 
}

template <typename coeff_type>
bool 
Stokhos::ProductContainer<coeff_type>::
myGID(int i) const
{
  return map_->MyGID(i);
}

template <typename coeff_type>
void
Stokhos::ProductContainer<coeff_type>::
init(const value_type& val)
{
  ordinal_type sz = coeff_.size();
  for (ordinal_type i=0; i<sz; i++)
    traits_type::init(*(coeff_[i]), val);
}

template <typename coeff_type>
std::ostream&
Stokhos::ProductContainer<coeff_type>::
print(std::ostream& os) const
{
  Teuchos::Array<ordinal_type> trm;
  ordinal_type sz = coeff_.size();
  os << "Stokhos::ProductContainer of global size " << map->NumGlobalElements()
     << ", local size " << sz << ":" << std::endl;
  for (ordinal_type i=0; i<sz; i++) {
    os << "Term " << map->GID(i) << ":" << std::endl;
    traits_type::print(os, *(coeff_[i]));
  }
  
  return os;
}
