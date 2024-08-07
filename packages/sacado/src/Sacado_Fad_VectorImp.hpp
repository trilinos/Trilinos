// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_DynamicArrayTraits.hpp"

template <typename OrdinalType, typename ValueType>
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >::
Vector(OrdinalType vec_size, OrdinalType deriv_sz, 
       VectorDerivOrientation orient) :
  deriv_size_(deriv_sz), orient_(orient), stride_(1), vec_(vec_size)
{
  ValueType* x = NULL;
  ValueType* dx_ = NULL;
  if (vec_size > 0) {
    x = ds_array<ValueType>::get_and_fill(vec_size);
    if (deriv_size_ > 0)
      dx_ = ds_array<ValueType>::get_and_fill(vec_size*deriv_size_);
  }
  if (orient_ == Row) {
    stride_ = vec_size;
    for (OrdinalType i=0; i<vec_size; i++)
      vec_[i].setMemory(deriv_size_, x+i, dx_+i, stride_);
  }
  else {
    stride_ = 1;
    for (OrdinalType i=0; i<vec_size; i++)
      vec_[i].setMemory(deriv_size_, x+i, dx_+i*deriv_size_, stride_);
  }
}

template <typename OrdinalType, typename ValueType>
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >::
Vector(const Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >& fv) :
  deriv_size_(fv.deriv_size_), orient_(fv.orient_), stride_(fv.stride_), 
  vec_(fv.size())
{
  OrdinalType vec_size = fv.size();
  ValueType* x = NULL;
  ValueType* dx_ = NULL;
  if (vec_size > 0) {
    x = ds_array<ValueType>::get_and_fill(vec_size);
    if (deriv_size_ > 0)
      dx_ = ds_array<ValueType>::get_and_fill(vec_size*deriv_size_);
  }
  if (orient_ == Row) {
    for (OrdinalType i=0; i<vec_size; i++) {
      vec_[i].setMemory(deriv_size_, x+i, dx_+i, stride_);
      vec_[i] = fv.vec_[i];
    }
  }
  else {
    for (OrdinalType i=0; i<vec_size; i++) {
      vec_[i].setMemory(deriv_size_, x+i, dx_+i*deriv_size_, stride_);
      vec_[i] = fv.vec_[i];
    }
  }
}

namespace Sacado {
namespace Fad {
template <typename OrdinalType, typename ValueType>
Vector< OrdinalType, DVFad<ValueType> >::
~Vector()
{
  // Here we must destroy the value and derivative arrays
  if (vec_.size() > 0) {
    ValueType *v = vals();
    ds_array<ValueType>::destroy_and_release(v, vec_.size());
    if (deriv_size_ > 0) {
      v = dx();
      ds_array<ValueType>::destroy_and_release(v, vec_.size()*deriv_size_);
    }
  }
}
}
}

template <typename OrdinalType, typename ValueType>
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >&
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >::
operator=(const Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >& fv) 
{ 
  vec_ = fv.vec_; 
  return *this; 
}

template <typename OrdinalType, typename ValueType>
ValueType* 
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >::
vals() 
{ 
  if (vec_.size() == 0)
    return NULL;
  return &(vec_[0].val()); 
}

template <typename OrdinalType, typename ValueType>
const ValueType* 
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >::
vals() const 
{ 
  if (vec_.size() == 0)
    return NULL;
  return &(vec_[0].val()); 
}

template <typename OrdinalType, typename ValueType>
ValueType* 
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >::
dx() 
{ 
  if (vec_.size() == 0 || deriv_size_ == 0)
    return NULL;
  return &(vec_[0].fastAccessDx(0)); 
}

template <typename OrdinalType, typename ValueType>
const ValueType* 
Sacado::Fad::Vector< OrdinalType, Sacado::Fad::DVFad<ValueType> >::
dx() const 
{ 
  if (vec_.size() == 0 || deriv_size_ == 0)
    return NULL;
  return &(vec_[0].fastAccessDx(0)); 
}
