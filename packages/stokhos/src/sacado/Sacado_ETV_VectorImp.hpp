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
VectorImpl(ordinal_type sz) :
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
    for (int i=0; i<sz; i++)
      (*th_)[i] = x.coeff(i);
  }
}

template <typename T, typename Storage> 
void
VectorImpl<T,Storage>::
reset(ordinal_type sz)
{
  th_->resize(sz);
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
  th_->resize(sz);
  if (x.hasFastAccess(sz)) {
    for (int i=0; i<sz; i++)
      (*th_)[i] = x.fastAccessCoeff(i);
  }
  else {
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
