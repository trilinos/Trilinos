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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
