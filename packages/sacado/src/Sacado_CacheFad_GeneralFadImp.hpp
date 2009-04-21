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
inline Sacado::CacheFad::GeneralFad<T,Storage>::GeneralFad(const Expr<S>& x) :
  Storage(T(0.))
{
  int sz = x.size();

  T xval = x.val(); 

  if (sz != this->size()) 
    this->resize(sz);

  if (sz) {
    if (x.hasFastAccess())
      for(int i=0; i<sz; ++i) 
	this->fastAccessDx(i) = x.fastAccessDx(i);
    else
      for(int i=0; i<sz; ++i) 
	this->fastAccessDx(i) = x.dx(i);
  }

  this->val() = xval;
}


template <typename T, typename Storage> 
inline void 
Sacado::CacheFad::GeneralFad<T,Storage>::diff(const int ith, const int n) 
{ 
  if (this->size() == 0) 
    this->resize(n);

  this->zero();
  this->fastAccessDx(ith) = T(1.);

}

template <typename T, typename Storage> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator=(const T& v) 
{
  this->val() = v;

  if (this->size()) {
    this->zero();    // We need to zero out the array for future resizes
    this->resize(0);
  }

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator=(
			     const Sacado::CacheFad::GeneralFad<T,Storage>& x) 
{
  // Copy val_ and dx_
  Storage::operator=(x);
  
  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator=(const Expr<S>& x) 
{
  int sz = x.size();

  T xval = x.val();

  if (sz != this->size()) 
    this->resize(sz);

  if (sz) {
    if (x.hasFastAccess())
      for(int i=0; i<sz; ++i)
	this->fastAccessDx(i) = x.fastAccessDx(i);
    else
      for(int i=0; i<sz; ++i)
	this->fastAccessDx(i) = x.dx(i);
  }
  
  this->val() = xval;
  
  return *this;
}

template <typename T, typename Storage> 
inline  Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator += (const T& v)
{
  this->val() += v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator -= (const T& v)
{
  this->val() -= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator *= (const T& v)
{
  int sz = this->size();

  this->val() *= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) *= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator /= (const T& v)
{
  int sz = this->size();

  this->val() /= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) /= v;

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator += (const Sacado::CacheFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  T xval = x.val();

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
      this->resize(xsz);
      if (x.hasFastAccess())
	for (int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = x.dx(i);
    }
  }

  this->val() += xval;

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator -= (const Sacado::CacheFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  T xval = x.val();

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
      this->resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = -x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = -x.dx(i);
    }
  }

  this->val() -= xval;


  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator *= (const Sacado::CacheFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();
  T xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->fastAccessDx(i) = this->val() * x.fastAccessDx(i) + this->fastAccessDx(i) * xval;
      else
	for (int i=0; i<sz; ++i)
	  this->fastAccessDx(i) = this->val() * x.dx(i) + this->fastAccessDx(i) * xval;
    }
    else {
      this->resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = this->val() * x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = this->val() * x.dx(i);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	this->fastAccessDx(i) *= xval;
    }
  }

  this->val() *= xval;

  return *this;
}

template <typename T, typename Storage>
template <typename S> 
inline Sacado::CacheFad::GeneralFad<T,Storage>& 
Sacado::CacheFad::GeneralFad<T,Storage>::operator /= (const Sacado::CacheFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();
  T xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - this->val()*x.fastAccessDx(i) )/ (xval*xval);
      else
	for (int i=0; i<sz; ++i)
	  this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - this->val()*x.dx(i) )/ (xval*xval);
    }
    else {
      this->resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = - this->val()*x.fastAccessDx(i) / (xval*xval);
      else
	for (int i=0; i<xsz; ++i)
	  this->fastAccessDx(i) = -this->val() * x.dx(i) / (xval*xval);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	this->fastAccessDx(i) /= xval;
    }
  }

  this->val() /= xval;

  return *this;
}

