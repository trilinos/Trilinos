// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Sacado_DynamicArrayTraits.hpp"

template <typename ValueT, typename ScalarT>
Sacado::Fad::FadVector<ValueT,ScalarT>::
FadVector(int vec_size, int deriv_sz) :
  deriv_size_(deriv_sz), vec_(vec_size)
{
  ValueT* x = NULL;
  ValueT* dx_ = NULL;
  if (vec_size > 0) {
    x = ds_array<ValueT>::get_and_fill(vec_size);
    if (deriv_size_ > 0)
      dx_ = ds_array<ValueT>::get_and_fill(vec_size*deriv_size_);
  }
  for (int i=0; i<vec_size; i++)
    vec_[i].setMemory(deriv_size_, x+i, dx_+i*deriv_size_);
}

template <typename ValueT, typename ScalarT>
Sacado::Fad::FadVector<ValueT,ScalarT>::
FadVector(const Sacado::Fad::FadVector<ValueT,ScalarT>& fv) :
  deriv_size_(fv.deriv_size_), vec_(fv.size())
{
  int vec_size = fv.size();
  ValueT* x = NULL;
  ValueT* dx_ = NULL;
  if (vec_size > 0) {
    x = ds_array<ValueT>::get_and_fill(vec_size);
    if (deriv_size_ > 0)
      dx_ = ds_array<ValueT>::get_and_fill(vec_size*deriv_size_);
  }
  for (int i=0; i<vec_size; i++) {
    vec_[i].setMemory(deriv_size_, x+i, dx_+i*deriv_size_);
    vec_[i] = fv.vec_[i];
  }
}

template <typename ValueT, typename ScalarT>
Sacado::Fad::FadVector<ValueT,ScalarT>::
~FadVector()
{
  // Here we must destroy the value and derivative arrays
  if (vec_.size() > 0) {
    ValueT *v = vals();
    ds_array<ValueT>::destroy_and_release(v, vec_.size());
    if (deriv_size_ > 0) {
      v = dx();
      ds_array<ValueT>::destroy_and_release(v, vec_.size()*deriv_size_);
    }
  }
}

template <typename ValueT, typename ScalarT>
Sacado::Fad::FadVector<ValueT,ScalarT>&
Sacado::Fad::FadVector<ValueT,ScalarT>::
operator=(const Sacado::Fad::FadVector<ValueT,ScalarT>& fv) 
{ 
  vec_ = fv.vec_; 
  return *this; 
}

template <typename ValueT, typename ScalarT>
ValueT* 
Sacado::Fad::FadVector<ValueT,ScalarT>::
vals() 
{ 
  if (vec_.size() == 0)
    return NULL;
  return &(vec_[0].val()); 
}

template <typename ValueT, typename ScalarT>
const ValueT* 
Sacado::Fad::FadVector<ValueT,ScalarT>::
vals() const 
{ 
  if (vec_.size() == 0)
    return NULL;
  return &(vec_[0].val()); 
}

template <typename ValueT, typename ScalarT>
ValueT* 
Sacado::Fad::FadVector<ValueT,ScalarT>::
dx() 
{ 
  if (vec_.size() == 0 || deriv_size_ == 0)
    return NULL;
  return &(vec_[0].fastAccessDx(0)); 
}

template <typename ValueT, typename ScalarT>
const ValueT* 
Sacado::Fad::FadVector<ValueT,ScalarT>::
dx() const 
{ 
  if (vec_.size() == 0 || deriv_size_ == 0)
    return NULL;
  return &(vec_[0].fastAccessDx(0)); 
}
