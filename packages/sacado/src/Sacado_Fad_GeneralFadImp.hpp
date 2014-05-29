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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#include "Sacado_ConfigDefs.h"

template <typename T, typename Storage>
template <typename S>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>::GeneralFad(const Expr<S>& x) :
  Storage(x.size(), T(0.)),
  update_val_(x.updateValue())
{
  const int sz = x.size();

  if (sz) {
    if (x.hasFastAccess())
      for(int i=0; i<sz; ++i)
        this->fastAccessDx(i) = x.fastAccessDx(i);
    else
      for(int i=0; i<sz; ++i)
        this->fastAccessDx(i) = x.dx(i);
  }

  if (update_val_)
    this->val() = x.val();
}


template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
void
Sacado::Fad::GeneralFad<T,Storage>::diff(const int ith, const int n)
{
  if (this->size() != n)
    this->resize(n);

  this->zero();
  this->fastAccessDx(ith) = T(1.);

}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator=(const T& v)
{
  this->val() = v;

  if (this->size())
    this->resize(0);

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator=(
                                  const Sacado::Fad::GeneralFad<T,Storage>& x)
{
  // Copy value & dx_
  Storage::operator=(x);
  update_val_ = x.update_val_;

  return *this;
}

template <typename T, typename Storage>
template <typename S>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator=(const Expr<S>& x)
{
  const int sz = x.size();

  if (sz != this->size())
    this->resizeAndZero(sz);

  if (sz) {
    if (x.hasFastAccess())
      for(int i=0; i<sz; ++i)
        this->fastAccessDx(i) = x.fastAccessDx(i);
    else
      for(int i=0; i<sz; ++i)
        this->fastAccessDx(i) = x.dx(i);
  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() = x.val();

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator += (const T& v)
{
  if (update_val_)
    this->val() += v;

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator -= (const T& v)
{
  if (update_val_)
    this->val() -= v;

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator *= (const T& v)
{
  const int sz = this->size();

  if (update_val_)
    this->val() *= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) *= v;

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator /= (const T& v)
{
  const int sz = this->size();

  if (update_val_)
    this->val() /= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) /= v;

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::
operator += (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  if (update_val_)
    this->val() += v;

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::
operator -= (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  if (update_val_)
    this->val() -= v;

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::
operator *= (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  const int sz = this->size();

  if (update_val_)
    this->val() *= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) *= v;

  return *this;
}

template <typename T, typename Storage>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::
operator /= (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  const int sz = this->size();

  if (update_val_)
    this->val() /= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) /= v;

  return *this;
}

template <typename T, typename Storage>
template <typename S>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator += (const Sacado::Fad::Expr<S>& x)
{
  const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) += x.fastAccessDx(i);
      else
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) += x.dx(i);
    }
    else {
      this->resizeAndZero(xsz);
      if (x.hasFastAccess())
        for (int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = x.fastAccessDx(i);
      else
        for (int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = x.dx(i);
    }
  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() += x.val();

  return *this;
}

template <typename T, typename Storage>
template <typename S>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator -= (const Sacado::Fad::Expr<S>& x)
{
  const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
        for(int i=0; i<sz; ++i)
          this->fastAccessDx(i) -= x.fastAccessDx(i);
      else
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) -= x.dx(i);
    }
    else {
      this->resizeAndZero(xsz);
      if (x.hasFastAccess())
        for(int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = -x.fastAccessDx(i);
      else
        for (int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = -x.dx(i);
    }
  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() -= x.val();


  return *this;
}

template <typename T, typename Storage>
template <typename S>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator *= (const Sacado::Fad::Expr<S>& x)
{
  const int xsz = x.size(), sz = this->size();
  T xval = x.val();
  T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
        for(int i=0; i<sz; ++i)
          this->fastAccessDx(i) = v*x.fastAccessDx(i) + this->fastAccessDx(i)*xval;
      else
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) = v*x.dx(i) + this->fastAccessDx(i)*xval;
    }
    else {
      this->resizeAndZero(xsz);
      if (x.hasFastAccess())
        for(int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = v*x.fastAccessDx(i);
      else
        for (int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = v*x.dx(i);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
        this->fastAccessDx(i) *= xval;
    }
  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() *= xval;

  return *this;
}

template <typename T, typename Storage>
template <typename S>
KOKKOS_INLINE_FUNCTION
Sacado::Fad::GeneralFad<T,Storage>&
Sacado::Fad::GeneralFad<T,Storage>::operator /= (const Sacado::Fad::Expr<S>& x)
{
  const int xsz = x.size(), sz = this->size();
  T xval = x.val();
  T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
        for(int i=0; i<sz; ++i)
          this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - v*x.fastAccessDx(i) )/ (xval*xval);
      else
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - v*x.dx(i) )/ (xval*xval);
    }
    else {
      this->resizeAndZero(xsz);
      if (x.hasFastAccess())
        for(int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = - v*x.fastAccessDx(i) / (xval*xval);
      else
        for (int i=0; i<xsz; ++i)
          this->fastAccessDx(i) = -v*x.dx(i) / (xval*xval);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
        this->fastAccessDx(i) /= xval;
    }
  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() /= xval;

  return *this;
}
