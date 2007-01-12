// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::Fad::GeneralFad<T,Storage>::GeneralFad(const Expr<S>& x)
{
  int sz = x.size();

  if (sz != s_.size()) 
    s_.resize(sz);

  if (sz) {
    if (x.hasFastAccess())
      for(int i=0; i<sz; ++i) 
	s_.dx_[i] = x.fastAccessDx(i);
    else
      for(int i=0; i<sz; ++i) 
	s_.dx_[i] = x.dx(i);
  }

  val_ = x.val();
}


template <typename T, typename Storage> 
inline void 
Sacado::Fad::GeneralFad<T,Storage>::diff(const int ith, const int n) 
{ 
  if (s_.size() == 0) 
    s_.resize(n);

  s_.zero();
  s_.dx_[ith] = T(1.);

}

template <typename T, typename Storage> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator=(const T& val) 
{
  val_ = val;

  if (s_.size()) 
    s_.resize(0);

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator=(
				  const Sacado::Fad::GeneralFad<T,Storage>& x) 
{
  // Copy value
  val_ = x.val_;

  // Copy dx_
  s_.operator=(x.s_);
  
  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator=(const Expr<S>& x) 
{
  int sz = x.size();

  if (sz != s_.size()) 
    s_.resize(sz);

  if (sz) {
    if (x.hasFastAccess())
      for(int i=0; i<sz; ++i)
	s_.dx_[i] = x.fastAccessDx(i);
    else
      for(int i=0; i<sz; ++i)
	s_.dx_[i] = x.dx(i);
  }
  
  val_ = x.val();
  
  return *this;
}

template <typename T, typename Storage> 
inline  Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator += (const T& val)
{
  val_ += val;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator -= (const T& val)
{
  val_ -= val;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator *= (const T& val)
{
  int sz = s_.size();

  val_ *= val;
  for (int i=0; i<sz; ++i)
    s_.dx_[i] *= val;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator /= (const T& val)
{
  int sz = s_.size();

  val_ /= val;
  for (int i=0; i<sz; ++i)
    s_.dx_[i] /= val;

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator += (const Sacado::Fad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for (int i=0; i<sz; ++i)
	  s_.dx_[i] += x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  s_.dx_[i] += x.dx(i);
    }
    else {
      s_.resize(xsz);
      if (x.hasFastAccess())
	for (int i=0; i<xsz; ++i)
	  s_.dx_[i] = x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  s_.dx_[i] = x.dx(i);
    }
  }

  val_ += x.val();

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator -= (const Sacado::Fad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  s_.dx_[i] -= x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  s_.dx_[i] -= x.dx(i);
    }
    else {
      s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  s_.dx_[i] = -x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  s_.dx_[i] = -x.dx(i);
    }
  }

  val_ -= x.val();


  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator *= (const Sacado::Fad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();
  T xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  s_.dx_[i] = val_ * x.fastAccessDx(i) + s_.dx_[i] * xval;
      else
	for (int i=0; i<sz; ++i)
	  s_.dx_[i] = val_ * x.dx(i) + s_.dx_[i] * xval;
    }
    else {
      s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  s_.dx_[i] = val_ * x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  s_.dx_[i] = val_ * x.dx(i);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	s_.dx_[i] *= xval;
    }
  }

  val_ *= xval;

  return *this;
}

template <typename T, typename Storage>
template <typename S> 
inline Sacado::Fad::GeneralFad<T,Storage>& 
Sacado::Fad::GeneralFad<T,Storage>::operator /= (const Sacado::Fad::Expr<S>& x)
{
  int xsz = x.size(), sz = s_.size();
  T xval = x.val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (xsz) {
    if (sz) {
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  s_.dx_[i] = ( s_.dx_[i]*xval - val_*x.fastAccessDx(i) )/ (xval*xval);
      else
	for (int i=0; i<sz; ++i)
	  s_.dx_[i] = ( s_.dx_[i]*xval - val_*x.dx(i) )/ (xval*xval);
    }
    else {
      s_.resize(xsz);
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  s_.dx_[i] = - val_*x.fastAccessDx(i) / (xval*xval);
      else
	for (int i=0; i<xsz; ++i)
	  s_.dx_[i] = -val_ * x.dx(i) / (xval*xval);
    }
  }
  else {
    if (sz) {
      for (int i=0; i<sz; ++i)
	s_.dx_[i] /= xval;
    }
  }

  val_ /= xval;

  return *this;
}

