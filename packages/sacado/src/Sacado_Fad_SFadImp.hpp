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

namespace Sacado {
namespace Fad {

#define FAD Expr< SFadExprTag<T,Num> >

template <typename T, int Num>
KOKKOS_INLINE_FUNCTION
Expr< SFadExprTag<T,Num> >::
Expr(const int sz, const T & x) : val_(x), update_val_(true)
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (sz != Num)
    throw "SFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
#endif

  ss_array<T>::zero(dx_, Num);
}

template <typename T, int Num>
KOKKOS_INLINE_FUNCTION
Expr< SFadExprTag<T,Num> >::
Expr(const int sz, const int i, const T & x) : val_(x), update_val_(true)
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (sz != Num)
    throw "SFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
  if (i >= Num)
    throw "SFad::SFad() Error:  Invalid derivative index.";
#endif

  ss_array<T>::zero(dx_, Num);
  dx_[i]=1.;
}

template <typename T, int Num>
KOKKOS_INLINE_FUNCTION
Expr< SFadExprTag<T,Num> >::
Expr(const Expr< SFadExprTag<T,Num> >& x) :
  val_(x.val()), update_val_(x.update_val_)
{
  for (int i=0; i<Num; i++)
    dx_[i] = x.dx_[i];
}

template <typename T, int Num>
template <typename S>
KOKKOS_INLINE_FUNCTION
Expr< SFadExprTag<T,Num> >::
Expr(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DEF)
  : update_val_(x.updateValue())
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (x.size() != Num)
    throw "SFad::SFad() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] = x.fastAccessDx(i);

  if (update_val_)
    this->val() = x.val();
  else
    this->val() = T(0.);
}


template <typename T, int Num>
KOKKOS_INLINE_FUNCTION
void
Expr< SFadExprTag<T,Num> >::
diff(const int ith, const int n)
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (n != Num)
    throw "SFad::diff() Error:  Supplied derivative dimension does not match compile time length.";
#endif

  ss_array<T>::zero(dx_, Num);
  dx_[ith] = T(1.);
}

template <typename T, int Num>
KOKKOS_INLINE_FUNCTION
void
Expr< SFadExprTag<T,Num> >::
resize(int sz)
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (sz != Num)
    throw "SFad::resize() Error:  Cannot resize fixed derivative array dimension";
#endif
}

template <typename T, int Num>
KOKKOS_INLINE_FUNCTION
Expr< SFadExprTag<T,Num> >&
Expr< SFadExprTag<T,Num> >::
operator=(const Expr< SFadExprTag<T,Num> >& x)
{
  // Copy value
  val_ = x.val_;

  // Copy dx_
  for (int i=0; i<Num; i++)
    dx_[i] = x.dx_[i];

  update_val_ = x.update_val_;

  return *this;
}

template <typename T, int Num>
template <typename S>
KOKKOS_INLINE_FUNCTION
SACADO_FAD_ENABLE_EXPR_FUNC
Expr< SFadExprTag<T,Num> >::
operator=(const Expr<S>& x)
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (x.size() != Num)
    throw "SFad::operator=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] = x.fastAccessDx(i);

  update_val_ = x.updateValue();
  if (update_val_)
    val_ = x.val();

  return *this;
}

template <typename T, int Num>
template <typename S>
KOKKOS_INLINE_FUNCTION
SACADO_FAD_ENABLE_EXPR_FUNC
Expr< SFadExprTag<T,Num> >::
operator += (const Expr<S>& x)
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (x.size() != Num)
    throw "SFad::operator+=() Error:  Attempt to assign with incompatible sizes";
#endif

  for (int i=0; i<Num; ++i)
    dx_[i] += x.fastAccessDx(i);

  update_val_ = x.updateValue();
  if (update_val_)
    val_ += x.val();

  return *this;
}

template <typename T, int Num>
template <typename S>
KOKKOS_INLINE_FUNCTION
SACADO_FAD_ENABLE_EXPR_FUNC
Expr< SFadExprTag<T,Num> >::
operator -= (const Expr<S>& x)
{
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (x.size() != Num)
    throw "SFad::operator-=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] -= x.fastAccessDx(i);

  update_val_ = x.updateValue();
  if (update_val_)
    val_ -= x.val();

  return *this;
}

template <typename T, int Num>
template <typename S>
KOKKOS_INLINE_FUNCTION
SACADO_FAD_ENABLE_EXPR_FUNC
Expr< SFadExprTag<T,Num> >::
operator *= (const Expr<S>& x)
{
  T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (x.size() != Num)
    throw "SFad::operator*=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] = val_ * x.fastAccessDx(i) + dx_[i] * xval;

  update_val_ = x.updateValue();
  if (update_val_)
    val_ *= xval;

  return *this;
}

template <typename T, int Num>
template <typename S>
KOKKOS_INLINE_FUNCTION
SACADO_FAD_ENABLE_EXPR_FUNC
Expr< SFadExprTag<T,Num> >::
operator /= (const Expr<S>& x)
{
  T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
  if (x.size() != Num)
    throw "SFad::operator/=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] = ( dx_[i]*xval - val_*x.fastAccessDx(i) )/ (xval*xval);

  update_val_ = x.updateValue();
  if (update_val_)
    val_ /= xval;

  return *this;
}

#undef FAD

} // namespace Fad
} // namespace Sacado
