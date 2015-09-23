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

#include "Stokhos_Version.hpp"

namespace Sacado {
namespace ETV {

template <typename T, typename Storage> 
VectorImpl<T,Storage>::
VectorImpl() :
  th_(new storage_type(1))
{
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>::
VectorImpl(const typename VectorImpl<T,Storage>::value_type& x) :
  th_(new storage_type(1))
{
  th_->init(x);
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>::
VectorImpl(ordinal_type sz, 
	   const typename VectorImpl<T,Storage>::value_type& x) :
  th_(new storage_type(sz))
{
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>::
VectorImpl(const VectorImpl<T,Storage>& x) :
  th_(x.th_)
{
}

template <typename T, typename Storage> 
template <typename S>
VectorImpl<T,Storage>::
VectorImpl(const Expr<S>& x) :
  th_(new storage_type(x.size()))
{
  int sz = x.size();
  if (x.hasFastAccess(sz)) {
    for (int i=0; i<sz; i++)
      (*th_)[i] = x.fastAccessCoeff(i);
  }
  else {
    Stokhos_Version();
    for (int i=0; i<sz; i++)
      (*th_)[i] = x.coeff(i);
  }
}

template <typename T, typename Storage> 
void
VectorImpl<T,Storage>::
reset(ordinal_type sz_new)
{
  int sz = this->size();
  th_->resize(sz_new);
  if (sz == 1 && sz_new > sz)
    for (int i=1; i<sz_new; i++)
      (*th_)[i] = (*th_)[0];
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>& 
VectorImpl<T,Storage>::
operator=(const typename VectorImpl<T,Storage>::value_type& v) 
{
  th_.makeOwnCopy();
  th_->init(v);
  return *this;
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>& 
VectorImpl<T,Storage>::
operator=(const VectorImpl<T,Storage>& x) 
{
  th_ = x.th_;
  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
VectorImpl<T,Storage>& 
VectorImpl<T,Storage>::
operator=(const Expr<S>& x) 
{
  th_.makeOwnCopy();
  int sz = x.size();
  this->reset(sz);
  if (x.hasFastAccess(sz)) {
    for (int i=0; i<sz; i++)
      (*th_)[i] = x.fastAccessCoeff(i);
  }
  else {
    Stokhos_Version();
    for (int i=0; i<sz; i++)
      (*th_)[i] = x.coeff(i);
  }
  return *this;
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>& 
VectorImpl<T,Storage>::
operator+=(const typename VectorImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  for (int i=0; i<th_->size(); i++)
    (*th_)[i] += v;
  return *this;
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>& 
VectorImpl<T,Storage>::
operator-=(const typename VectorImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  for (int i=0; i<th_->size(); i++)
    (*th_)[i] -= v;
  return *this;
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>& 
VectorImpl<T,Storage>::
operator*=(const typename VectorImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  for (int i=0; i<th_->size(); i++)
    (*th_)[i] *= v;
  return *this;
}

template <typename T, typename Storage> 
VectorImpl<T,Storage>& 
VectorImpl<T,Storage>::
operator/=(const typename VectorImpl<T,Storage>::value_type& v)
{
  th_.makeOwnCopy();
  for (int i=0; i<th_->size(); i++)
    (*th_)[i] /= v;
  return *this;
}

template <typename T, typename Storage>
std::ostream& 
operator << (std::ostream& os, const Vector<T,Storage>& a)
{
  typedef typename Vector<T,Storage>::ordinal_type ordinal_type;

  os << "[ ";
      
  for (ordinal_type i=0; i<a.size(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}

} // namespace PCE
} // namespace Sacado
