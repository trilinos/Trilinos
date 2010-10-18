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

#include <thrust/fill.h>

template <typename ordinal_type, typename value_type>
Stokhos::CUDAStorage<ordinal_type, value_type>::
CUDAStorage(const ordinal_type& sz) : 
  coeff_(sz, value_type(0.0)) 
{
}

template <typename ordinal_type, typename value_type>
Stokhos::CUDAStorage<ordinal_type, value_type>::
CUDAStorage(const CUDAStorage<ordinal_type, value_type>& s) : 
  coeff_(s.coeff_.begin(), s.coeff_.end()) 
{
}

template <typename ordinal_type, typename value_type>
Stokhos::CUDAStorage<ordinal_type, value_type>&
Stokhos::CUDAStorage<ordinal_type, value_type>::
operator=(const Stokhos::CUDAStorage<ordinal_type, value_type>& s)
{
  if (this != &s) {
    coeff_.resize(s.coeff_.size());
    thrust::copy(s.coeff_.begin(), s.coeff_.end(), coeff_.begin());
  }
  return *this;
}

template <typename ordinal_type, typename value_type>
Stokhos::CUDAStorage<ordinal_type, value_type>::
~CUDAStorage() 
{
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CUDAStorage<ordinal_type, value_type>::
init(const value_type& v) 
{ 
  thrust::fill(coeff_.begin(), coeff_.end(), v); 
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CUDAStorage<ordinal_type, value_type>::
init(const value_type* v, const ordinal_type& sz) 
{ 
  if (sz == 0)
    thrust::copy(v, v+coeff_.size(), coeff_.begin()); 
  else
    thrust::copy(v, v+sz, coeff_.begin()); 
}

template <typename ordinal_type, typename value_type>
void
Stokhos::CUDAStorage<ordinal_type, value_type>::
load(value_type* v) 
{ 
  thrust::copy(coeff_.begin(), coeff_.end(), v); 
}

template <typename ordinal_type, typename value_type>
void 
Stokhos::CUDAStorage<ordinal_type, value_type>::
resize(const ordinal_type& sz) 
{ 
  coeff_.resize(sz); 
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::CUDAStorage<ordinal_type, value_type>::
size() const 
{ 
  return coeff_.size(); 
}

template <typename ordinal_type, typename value_type>
typename Stokhos::CUDAStorage<ordinal_type, value_type>::const_reference
Stokhos::CUDAStorage<ordinal_type, value_type>::
operator[] (const ordinal_type& i) const 
{
  return coeff_[i];
}

template <typename ordinal_type, typename value_type>
typename Stokhos::CUDAStorage<ordinal_type, value_type>::reference
Stokhos::CUDAStorage<ordinal_type, value_type>::
operator[] (const ordinal_type& i) 
{
  return coeff_[i];
}

template <typename ordinal_type, typename value_type>
typename Stokhos::CUDAStorage<ordinal_type, value_type>::const_pointer
Stokhos::CUDAStorage<ordinal_type, value_type>::
coeff() const 
{ 
  return coeff_.data();
}

template <typename ordinal_type, typename value_type>
typename Stokhos::CUDAStorage<ordinal_type, value_type>::pointer
Stokhos::CUDAStorage<ordinal_type, value_type>::
coeff() 
{ 
  return coeff_.data();
}

template <typename ordinal_type, typename value_type>
thrust::device_vector<value_type>&
Stokhos::CUDAStorage<ordinal_type, value_type>::
dev_vector() 
{ 
  return coeff_;
}

template <typename ordinal_type, typename value_type>
const thrust::device_vector<value_type>&
Stokhos::CUDAStorage<ordinal_type, value_type>::
dev_vector() const
{ 
  return coeff_;
}
