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

#include "Stokhos_DynamicArrayTraits.hpp"

template <typename T> 
Stokhos::HermitePoly<T>::
HermitePoly() :
  coeff_(NULL), deg_(-1), len_(0) 
{
}

template <typename T> 
Stokhos::HermitePoly<T>::
HermitePoly(const T& x) :
  coeff_(), deg_(0), len_(1) 
{
  coeff_ = Stokhos::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T> 
Stokhos::HermitePoly<T>::
HermitePoly(unsigned int d, const T& x) :
  coeff_(), deg_(d), len_(d+1) 
{
  coeff_ = Stokhos::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T> 
Stokhos::HermitePoly<T>::
HermitePoly(unsigned int d) :
  coeff_(), deg_(d), len_(d+1) 
{
  coeff_ = Stokhos::ds_array<T>::get_and_fill(len_);
}

template <typename T> 
Stokhos::HermitePoly<T>::
HermitePoly(unsigned int d, unsigned int l) :
  coeff_(), deg_(d), len_(l) 
{
  coeff_ = Stokhos::ds_array<T>::get_and_fill(len_);
}

template <typename T> 
Stokhos::HermitePoly<T>::
HermitePoly(const Stokhos::HermitePoly<T>& x) :
  coeff_(), deg_(x.deg_), len_(x.deg_+1) 
{
  coeff_ = Stokhos::ds_array<T>::get_and_fill(x.coeff_, len_);
}

template <typename T> 
Stokhos::HermitePoly<T>::
~HermitePoly()
{
  if (len_ > 0)
    Stokhos::ds_array<T>::destroy_and_release(coeff_, len_);
}

template <typename T> 
Stokhos::HermitePoly<T>&
Stokhos::HermitePoly<T>::
operator=(const Stokhos::HermitePoly<T>& x) 
{
  if (len_ < x.deg_+1) {
    Stokhos::ds_array<T>::destroy_and_release(coeff_, len_);
    len_ = x.deg_+1;
    deg_ = x.deg_;
    coeff_ = Stokhos::ds_array<T>::get_and_fill(x.coeff_, len_);
  }
  else {
    deg_ = x.deg_;
    Stokhos::ds_array<T>::copy(x.coeff_, coeff_, deg_+1);
  }
  
  return *this;
}

template <typename T> 
void
Stokhos::HermitePoly<T>::
resize(unsigned int d, bool keep_coeffs)
{
  if (d+1 > len_) {
    T* c = Stokhos::ds_array<T>::get_and_fill(d+1);
    if (keep_coeffs)
      Stokhos::ds_array<T>::copy(coeff_, c, deg_+1);
    Stokhos::ds_array<T>::destroy_and_release(coeff_, len_);
    coeff_ = c;
    deg_ = d;
    len_ = d+1;
  }
  else {
    if (!keep_coeffs)
      Stokhos::ds_array<T>::zero(coeff_, deg_+1);
    deg_ = d;
  }
}

template <typename T> 
void
Stokhos::HermitePoly<T>::
reserve(unsigned int d)
{
  if (d+1 > len_) {
    T* c = Stokhos::ds_array<T>::get_and_fill(d+1);
    Stokhos::ds_array<T>::copy(coeff_, c, deg_+1);
    Stokhos::ds_array<T>::destroy_and_release(coeff_, len_);
    coeff_ = c;
    len_ = d+1;
  }
}

template <typename T>
std::ostream& 
Stokhos::operator << (std::ostream& os, const Stokhos::HermitePoly<T>& a)
{
  os << "[ ";
      
  for (unsigned int i=0; i<=a.degree(); i++) {
    os << a[i] << " ";
  }

  os << "]\n";
  return os;
}
