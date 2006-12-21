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
// @HEADER

template <typename T> 
template <typename S> 
inline Sacado::Taylor::DTaylor<T>::DTaylor(const Expr<S>& x) : 
  Expr< DTaylorImplementation<T> >(x.degree(), T(0))
{
  unsigned int d = this->degree();

  x.allocateCache(d);

  // We copy the coefficients from the highest degree to the lowest just
  // to be consistent with operator=(), even though it is not strictly 
  // necessary since "this" cannot be on the RHS in the copy constructor.
  if (x.hasFastAccess(d))
    for(int i=d; i>=0; --i) 
      this->coeff_[i] = x.fastAccessCoeff(i);
  else
    for(int i=d; i>=0; --i) 
      this->coeff_[i] = x.coeff(i);

}

template <typename T> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator=(const T& val) 
{
  this->coeff_[0] = val;

  for (unsigned int i=1; i<this->coeff_.size(); i++)
    this->coeff_[i] = T(0.);

  return *this;
}

template <typename T> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator=(const Sacado::Taylor::DTaylor<T>& x) 
{
  if (x.coeff_.size() != this->coeff_.size()) 
    this->coeff_.resize(x.coeff_.size());
  this->coeff_ = x.coeff_;
  
  return *this;
}

template <typename T> 
template <typename S> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator=(const Expr<S>& x) 
{
  unsigned int d = this->degree();
  unsigned int xd = x.degree();

  // Resize polynomial for "this" if x has greater degree
  if (xd > d) {
    this->coeff_.resize(xd+1);
    d = xd;
  }

  x.allocateCache(d);

  // Copy coefficients.  Note:  we copy from the last term down to the first
  // to take into account "this" being involved in an expression on the RHS.
  // By overwriting the degree k term, we are guarranteed not to affect any
  // of the lower degree terms.  However, this in general will force each
  // term in the expression to compute all of its coefficients at once instead
  // traversing the expression once for each degree.
  if (x.hasFastAccess(d))
    for(int i=d; i>=0; --i) 
      this->coeff_[i] = x.fastAccessCoeff(i);
    else
      for(int i=d; i>=0; --i) 
	this->coeff_[i] = x.coeff(i);
  
  return *this;
}

template <typename T> 
inline  Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator += (const T& val)
{
  this->coeff_[0] += val;

  return *this;
}

template <typename T> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator -= (const T& val)
{
  this->coeff_[0] -= val;

  return *this;
}

template <typename T> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator *= (const T& val)
{
  this->coeff_ *= val;

  return *this;
}

template <typename T> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator /= (const T& val)
{
  this->coeff_ /= val;

  return *this;
}

template <typename T> 
template <typename S> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator += (const S& x)
{
  unsigned int xd = x.degree();
  unsigned int d = this->degree();

   // Resize polynomial for "this" if x has greater degree
  if (xd > d) {
    this->resizeCoeffs(xd);
    d = xd;
  }

  x.allocateCache(d);

  if (x.hasFastAccess(d))
    for (int i=d; i>=0; --i)
      this->coeff_[i] += x.fastAccessCoeff(i);
  else
    for (int i=xd; i>=0; --i)
      this->coeff_[i] += x.coeff(i);

  return *this;
}

template <typename T> 
template <typename S> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator -= (const S& x)
{
  unsigned int xd = x.degree();
  unsigned int d = this->degree();

   // Resize polynomial for "this" if x has greater degree
  if (xd > d) {
    this->resizeCoeffs(xd);
    d = xd;
  }

  x.allocateCache(d);

  if (x.hasFastAccess(d))
    for (int i=d; i>=0; --i)
      this->coeff_[i] -= x.fastAccessCoeff(i);
  else
    for (int i=xd; i>=0; --i)
      this->coeff_[i] -= x.coeff(i);

  return *this;
}

template <typename T> 
template <typename S> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator *= (const S& x)
{
  unsigned int xd = x.degree();
  unsigned int d = this->degree();
  unsigned int dfinal = d;

   // Resize polynomial for "this" if x has greater degree
  if (xd > d) {
    this->resizeCoeffs(xd);
    dfinal = xd;
  }

  x.allocateCache(dfinal);

  if (xd) {
    if (d) {
      T tmp;
      if (x.hasFastAccess(dfinal))
	for(int i=dfinal; i>=0; --i) {
	  tmp = T(0);
	  for (int k=0; k<=i; ++k)
	    tmp += this->coeff_[k]*x.fastAccessCoeff(i-k);
	  this->coeff_[i] = tmp;
	}
      else
	for(int i=dfinal; i>=0; --i) {
	  tmp = T(0);
	  for (int k=0; k<=i; ++k)
	    tmp += this->coeff_[k]*x.coeff(i-k);
	  this->coeff_[i] = tmp;
	}
    }
    else {
      if (x.hasFastAccess(dfinal))
	for(int i=dfinal; i>=0; --i)
	  this->coeff_[i] = this->coeff_[0] * x.fastAccessCoeff(i);
      else
	for(int i=dfinal; i>=0; --i)
	  this->coeff_[i] = this->coeff_[0] * x.coeff(i);
    }
  }
  else 
    this->coeff_ *= x.coeff(0);

  return *this;
}

template <typename T>
template <typename S> 
inline Sacado::Taylor::DTaylor<T>& 
Sacado::Taylor::DTaylor<T>::operator /= (const S& x)
{
  unsigned int xd = x.degree();
  unsigned int d = this->degree();
  unsigned int dfinal = d;

   // Resize polynomial for "this" if x has greater degree
  if (xd > d) {
    this->resizeCoeffs(xd);
    dfinal = xd;
  }

  x.allocateCache(dfinal);

  if (xd) {
    std::valarray<T> tmp(this->coeff_);
    if (x.hasFastAccess(dfinal))
      for(unsigned int i=0; i<=dfinal; i++) {
	for (unsigned int k=1; k<=i; k++)
	  tmp[i] -= x.fastAccessCoeff(k)*tmp[i-k];
	tmp[i] /= x.fastAccessCoeff(0);
      }
    else
      for(unsigned int i=0; i<=dfinal; i++) {
	for (unsigned int k=1; k<=i; k++)
	  tmp[i] -= x.coeff(k)*tmp[i-k];
	tmp[i] /= x.coeff(0);
      }
    this->coeff_ = tmp;
  }
  else
    this->coeff_ /= x.coeff(0);

  return *this;
}

