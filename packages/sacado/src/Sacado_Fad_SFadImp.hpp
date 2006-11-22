// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

template <typename T, int Num> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
Expr(const int sz, const T & x) : val_(x)
{
#ifdef SACADO_DEBUG
  if (sz != Num)
    throw "SFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
#endif

  ss_array<T>::zero(dx_, Num); 
}

template <typename T, int Num> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
Expr(const int sz, const int i, const T & x) : val_(x) 
{
#ifdef SACADO_DEBUG
  if (sz != Num)
    throw "SFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
  if (i >= Num)
    throw "SFad::SFad() Error:  Invalid derivative index.";
#endif

  ss_array<T>::zero(dx_, Num);
  dx_[i]=1.; 
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
Expr(const Expr<S>& x) : val_(x.val())
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SFad::SFad() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i) 
    dx_[i] = x.fastAccessDx(i);
}


template <typename T, int Num> 
inline void 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
diff(const int ith, const int n) 
{ 
#ifdef SACADO_DEBUG
  if (n != Num)
    throw "SFad::diff() Error:  Supplied derivative dimension does not match compile time length.";
#endif

  ss_array<T>::zero(dx_, Num);
  dx_[ith] = T(1.);
}

template <typename T, int Num> 
inline void 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
resize(int sz)
{
#ifdef SACADO_DEBUG
  if (sz != Num)
    throw "SFad::resize() Error:  Cannot resize fixed derivative array dimension";
#endif
}

template <typename T, int Num> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator=(const T& val) 
{
  val_ = val;
  ss_array<T>::zero(dx_, Num);

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator=(const Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& x) 
{
  // Copy value
  val_ = x.val_;

  // Copy dx_
  ss_array<T>::copy(x.dx_, dx_, Num);
  
  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator=(const Expr<S>& x) 
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SFad::operator=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] = x.fastAccessDx(i);
  
  val_ = x.val();
  
  return *this;
}

template <typename T, int Num> 
inline  Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator += (const T& val)
{
  val_ += val;

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator -= (const T& val)
{
  val_ -= val;

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator *= (const T& val)
{
  val_ *= val;

  for (int i=0; i<Num; ++i)
    dx_[i] *= val;

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator /= (const T& val)
{
  val_ /= val;

  for (int i=0; i<Num; ++i)
    dx_[i] /= val;

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator += (const Sacado::Fad::Expr<S>& x)
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SFad::operator+=() Error:  Attempt to assign with incompatible sizes";
#endif

  for (int i=0; i<Num; ++i)
    dx_[i] += x.fastAccessDx(i);
 
  val_ += x.val();

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator -= (const Sacado::Fad::Expr<S>& x)
{
#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SFad::operator-=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] -= x.fastAccessDx(i);

  val_ -= x.val();

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator *= (const Sacado::Fad::Expr<S>& x)
{
  T xval = x.val();

#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SFad::operator*=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] = val_ * x.fastAccessDx(i) + dx_[i] * xval;
 
  val_ *= xval;

  return *this;
}

template <typename T, int Num>
template <typename S> 
inline Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >& 
Sacado::Fad::Expr< Sacado::Fad::SFadExprTag<T,Num> >::
operator /= (const Sacado::Fad::Expr<S>& x)
{
  T xval = x.val();

#ifdef SACADO_DEBUG
  if (x.size() != Num)
    throw "SFad::operator/=() Error:  Attempt to assign with incompatible sizes";
#endif

  for(int i=0; i<Num; ++i)
    dx_[i] = ( dx_[i]*xval - val_*x.fastAccessDx(i) )/ (xval*xval);

  val_ /= xval;

  return *this;
}

