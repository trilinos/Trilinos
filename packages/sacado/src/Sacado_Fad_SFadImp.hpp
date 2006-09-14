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
template <typename S> 
inline Sacado::Fad::SFad<T,Num>::SFad(const Expr<S>& x)
{
  for(int i=0; i<Num; ++i) 
    this->dx_[i] = x.dx(i);

  this->val_ = x.val();
}


template <typename T, int Num> 
inline void 
Sacado::Fad::SFad<T,Num>::diff(const int ith, const int n) 
{ 
  // Cannot just say "this->dx_ = T(0.)" since dx_ is an array, not a std::valarray
  this->dx_[0] = T(0.);
  for(int i = 1; i < n; i++)
	this->dx_[i] = this->dx_[0];
  this->dx_[ith] = T(1.);

}

template <typename T, int Num> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator=(const T& val) 
{
  this->val_ = val;
  for (int i=0; i<Num; i++)
    this->dx_[i] = T(0.0);

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator=(const Sacado::Fad::SFad<T,Num>& x) 
{
  // Copy value
  this->val_ = x.val_;

  // Copy dx_
  for (int i=0; i<Num; i++)
    this->dx_[i] = x.dx_[i];
  
  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator=(const Expr<S>& x) 
{
  for (int i=0; i<Num; i++)
    this->dx_[i] = x.dx(i);
  
  this->val_ = x.val();
  
  return *this;
}

template <typename T, int Num> 
inline  Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator += (const T& val)
{
  this->val_ += val;

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator -= (const T& val)
{
  this->val_ -= val;

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator *= (const T& val)
{
  this->val_ *= val;
  for (int i=0; i<Num; i++)
    this->dx_[i] *= val;

  return *this;
}

template <typename T, int Num> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator /= (const T& val)
{
  this->val_ /= val;
  for (int i=0; i<Num; i++)
    this->dx_[i] /= val;

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator += (const S& x)
{
  for (int i=0; i<Num; i++)
    this->dx_[i] += x.dx(i);

  this->val_ += x.val();

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator -= (const S& x)
{
  for (int i=0; i<Num; i++)
    this->dx_[i] -= x.dx(i);

  this->val_ -= x.val();

  return *this;
}

template <typename T, int Num> 
template <typename S> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator *= (const S& x)
{
  T xval = x.val();

  for (int i=0; i<Num; i++)
    this->dx_[i] = this->val_ * x.dx(i) + this->dx_[i] * xval;

  this->val_ *= xval;

  return *this;
}

template <typename T, int Num>
template <typename S> 
inline Sacado::Fad::SFad<T,Num>& 
Sacado::Fad::SFad<T,Num>::operator /= (const S& x)
{
  T xval = x.val();

  for (int i=0; i<Num; ++i)
    this->dx_[i] = ( this->dx_[i]*xval - this->val_*x.dx(i) )/ (xval*xval);

  this->val_ /= xval;

  return *this;
}

