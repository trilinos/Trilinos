// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

template <typename T>
template <typename S>
inline Sacado::Tay::CacheTaylor<T>::CacheTaylor(const Expr<S>& x) :
  Expr< CacheTaylorImplementation<T> >(x.degree(), T(0.))
{
  int d = this->degree();

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
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator=(const T& v)
{
  this->coeff_[0] = v;

  for (int i=1; i<this->coeff_size(); i++)
    this->coeff_[i] = T(0.);

  return *this;
}

template <typename T>
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator=(const Sacado::Tay::CacheTaylor<T>& x)
{
  if (x.coeff_size() != this->coeff_size())
    this->coeff_.resize(x.coeff_size());
  this->coeff_ = x.coeff_;

  return *this;
}

template <typename T>
template <typename S>
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator=(const Expr<S>& x)
{
  int d = this->degree();
  int xd = x.degree();

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
inline  Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator += (const T& v)
{
  this->coeff_[0] += v;

  return *this;
}

template <typename T>
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator -= (const T& v)
{
  this->coeff_[0] -= v;

  return *this;
}

template <typename T>
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator *= (const T& v)
{
  this->coeff_ *= v;

  return *this;
}

template <typename T>
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator /= (const T& v)
{
  this->coeff_ /= v;

  return *this;
}

template <typename T>
template <typename S>
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator += (const S& x)
{
  int xd = x.degree();
  int d = this->degree();

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
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator -= (const S& x)
{
  int xd = x.degree();
  int d = this->degree();

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
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator *= (const S& x)
{
  int xd = x.degree();
  int d = this->degree();
  int dfinal = d;

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
          tmp = T(0.);
          for (int k=0; k<=i; ++k)
            tmp += this->coeff_[k]*x.fastAccessCoeff(i-k);
          this->coeff_[i] = tmp;
        }
      else
        for(int i=dfinal; i>=0; --i) {
          tmp = T(0.);
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
inline Sacado::Tay::CacheTaylor<T>&
Sacado::Tay::CacheTaylor<T>::operator /= (const S& x)
{
  int xd = x.degree();
  int d = this->degree();
  int dfinal = d;

   // Resize polynomial for "this" if x has greater degree
  if (xd > d) {
    this->resizeCoeffs(xd);
    dfinal = xd;
  }

  x.allocateCache(dfinal);

  if (xd) {
    std::valarray<T> tmp(this->coeff_);
    if (x.hasFastAccess(dfinal))
      for(int i=0; i<=dfinal; i++) {
        for (int k=1; k<=i; k++)
          tmp[i] -= x.fastAccessCoeff(k)*tmp[i-k];
        tmp[i] /= x.fastAccessCoeff(0);
      }
    else
      for(int i=0; i<=dfinal; i++) {
        for (int k=1; k<=i; k++)
          tmp[i] -= x.coeff(k)*tmp[i-k];
        tmp[i] /= x.coeff(0);
      }
    this->coeff_ = tmp;
  }
  else
    this->coeff_ /= x.coeff(0);

  return *this;
}

