// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
ProductContainer(const Teuchos::RCP<const Epetra_BlockMap>& theMap) :
  map_(theMap),
  coeff_(map_->NumMyElements())
{
}

template <typename coeff_type>
Stokhos::ProductContainer<coeff_type>::
ProductContainer(const Teuchos::RCP<const Epetra_BlockMap>& theMap,
                 const typename traits_type::cloner_type& cloner) :
  map_(theMap),
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
reset(const Teuchos::RCP<const Epetra_BlockMap>& theMap)
{
  map_ = theMap;
  ordinal_type sz = map_->NumMyElements();
  coeff_.resize(sz);
}

template <typename coeff_type>
void
Stokhos::ProductContainer<coeff_type>::
reset(const Teuchos::RCP<const Epetra_BlockMap>& theMap,
      const typename traits_type::cloner_type& cloner)
{
  map_ = theMap;
  ordinal_type sz = map_->NumMyElements();
  coeff_.resize(sz);
  for (ordinal_type i=0; i<sz; i++)
    coeff_[i] = cloner.clone(i);
}

template <typename coeff_type>
void
Stokhos::ProductContainer<coeff_type>::
resize(const Teuchos::RCP<const Epetra_BlockMap>& theMap)
{
  map_ = theMap;
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
  os << "Stokhos::ProductContainer of global size " << map_->NumGlobalElements()
     << ", local size " << sz << ":" << std::endl;
  for (ordinal_type i=0; i<sz; i++) {
    os << "Term " << map_->GID(i) << ":" << std::endl;
    traits_type::print(os, *(coeff_[i]));
  }

  return os;
}
