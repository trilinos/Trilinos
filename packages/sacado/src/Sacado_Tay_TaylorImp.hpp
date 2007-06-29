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

#include "Sacado_DynamicArrayTraits.hpp"

template <typename T> 
Sacado::Tay::Taylor<T>::TaylorData::
TaylorData() :
  coeff_(NULL), deg_(-1), len_(0) 
{
}

template <typename T> 
Sacado::Tay::Taylor<T>::TaylorData::
TaylorData(const T& x) :
  coeff_(), deg_(0), len_(1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T> 
Sacado::Tay::Taylor<T>::TaylorData::
TaylorData(unsigned int d, const T& x) :
  coeff_(), deg_(d), len_(d+1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T> 
Sacado::Tay::Taylor<T>::TaylorData::
TaylorData(unsigned int d) :
  coeff_(), deg_(d), len_(d+1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
}

template <typename T> 
Sacado::Tay::Taylor<T>::TaylorData::
TaylorData(const typename Sacado::Tay::Taylor<T>::TaylorData& x) :
  coeff_(), deg_(x.deg_), len_(x.deg_+1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(x.coeff_, len_);
}

template <typename T> 
Sacado::Tay::Taylor<T>::TaylorData::
~TaylorData()
{
  if (len_ > 0)
    Sacado::ds_array<T>::destroy_and_release(coeff_, len_);
}

template <typename T> 
typename Sacado::Tay::Taylor<T>::TaylorData&
Sacado::Tay::Taylor<T>::TaylorData::
operator=(const typename Sacado::Tay::Taylor<T>::TaylorData& x) 
{
  if (len_ < x.deg_+1) {
    Sacado::ds_array<T>::destroy_and_release(coeff_, len_);
    len_ = x.deg_+1;
    deg_ = x.deg_;
    coeff_ = Sacado::ds_array<T>::get_and_fill(x.coeff_, len_);
  }
  else {
    deg_ = x.deg_;
    Sacado::ds_array<T>::copy(x.coeff_, coeff_, deg_+1);
  }
  
  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>::
Taylor() :
  th(new TaylorData)
{
}

template <typename T> 
Sacado::Tay::Taylor<T>::
Taylor(const T& x) :
  th(new TaylorData(x))
{
}

template <typename T> 
Sacado::Tay::Taylor<T>::
Taylor(unsigned int d, const T& x) :
  th(new TaylorData(d, x))
{
}

template <typename T> 
Sacado::Tay::Taylor<T>::
Taylor(unsigned int d) :
  th(new TaylorData(d))
{
}

template <typename T> 
Sacado::Tay::Taylor<T>::
Taylor(const Sacado::Tay::Taylor<T>& x) :
  th(x.th)
{
}

template <typename T> 
Sacado::Tay::Taylor<T>::
~Taylor()
{
}

template <typename T> 
void
Sacado::Tay::Taylor<T>::
resize(unsigned int d)
{
  if (d+1 > length()) {
    Sacado::Handle<TaylorData> h(new TaylorData(d));
    th = h;
  }
  else
    Sacado::ds_array<T>::zero(coeff(), degree()+1);
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator=(const T& val) 
{
  th.makeOwnCopy();

  if (th->len_ == 0) {
    th->len_ = 1;
    th->deg_ = 0;
    th->coeff_ = Sacado::ds_array<T>::get_and_fill(th->len_);
  }

  th->coeff_[0] = val;
  Sacado::ds_array<T>::zero(th->coeff_+1, th->deg_);

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator=(const Sacado::Tay::Taylor<T>& x) 
{
  th = x.th;
  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>
Sacado::Tay::Taylor<T>::
operator+() const
{
  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T> 
Sacado::Tay::Taylor<T>::
operator-() const
{
  Taylor<T> x;
  x.th->deg_ = th->deg_;
  x.th->len_ = th->deg_+1;
  x.th->coeff_ = Sacado::ds_array<T>::get_and_fill(x.th->len_);
  for (unsigned int i=0; i<=th->deg_; i++)
    x.th->coeff_[i] = -th->coeff_[i];

  return x;
}

template <typename T> 
 Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator+=(const T& val)
{
  th.makeOwnCopy();

  th->coeff_[0] += val;

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator-=(const T& val)
{
  th.makeOwnCopy();

  th->coeff_[0] -= val;

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator*=(const T& val)
{
  th.makeOwnCopy();

  for (unsigned int i=0; i<=th->deg_; i++)
    th->coeff_[i] *= val;

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator/=(const T& val)
{
  th.makeOwnCopy();

  for (unsigned int i=0; i<=th->deg_; i++)
    th->coeff_[i] /= val;

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator+=(const Taylor<T>& x)
{
  th.makeOwnCopy();

  unsigned int d = degree();
  unsigned int xd = x.degree();
  unsigned int dmin = xd < d ? xd : d;

  unsigned int l = xd + 1;
  bool need_resize = l > length();
  T* c;
  if (need_resize) {
    c = Sacado::ds_array<T>::get_and_fill(l);
    for (unsigned int i=0; i<=d; i++)
      c[i] = fastAccessCoeff(i);
  }
  else
    c = th->coeff_;
  T* xc = x.th->coeff_;

  for (unsigned int i=0; i<=dmin; i++)
    c[i] += xc[i];
  if (th->deg_ < xd) {
    for (unsigned int i=d+1; i<=xd; i++)
      c[i] = xc[i];
    th->deg_ = xd;
  }

  if (need_resize) {
    Sacado::ds_array<T>::destroy_and_release(th->coeff_, th->len_);
    th->len_ = l;
    th->coeff_ = c;
  }

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator-=(const Taylor<T>& x)
{
  th.makeOwnCopy();

  unsigned int d = degree();
  unsigned int xd = x.degree();
  unsigned int dmin = xd < d ? xd : d;

  unsigned int l = xd + 1;
  bool need_resize = l > length();
  T* c;
  if (need_resize) {
    c = Sacado::ds_array<T>::get_and_fill(l);
    for (unsigned int i=0; i<=d; i++)
      c[i] = fastAccessCoeff(i);
  }
  else
    c = th->coeff_;
  T* xc = x.th->coeff_;

  for (unsigned int i=0; i<=dmin; i++)
    c[i] -= xc[i];
  if (d < xd) {
    for (unsigned int i=d+1; i<=xd; i++)
      c[i] = -xc[i];
    th->deg_ = xd;
  }

  if (need_resize) {
    Sacado::ds_array<T>::destroy_and_release(th->coeff_, th->len_);
    th->len_ = l;
    th->coeff_ = c;
  }

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator*=(const Taylor<T>& x)
{
  unsigned int d = degree();
  unsigned int xd = x.degree();

#ifdef SACADO_DEBUG
  if ((xd != d) && (xd != 0) && (d != 0))
    throw "Taylor Error:  Attempt to assign with incompatible degrees";
#endif
  
  T* c = th->coeff_; 
  T* xc = x.th->coeff_;
  
  if (d > 0 && xd > 0) {
    Sacado::Handle<TaylorData> h(new TaylorData(d));
    T* cc = h->coeff_;
    T tmp;
    for (unsigned int i=0; i<=d; i++) {
      tmp = T(0.0);
      for (unsigned int k=0; k<=i; ++k)
	tmp += c[k]*xc[i-k];
      cc[i] = tmp;
    }
    th = h;
  }
  else if (d > 0) {
    th.makeOwnCopy();
    for (unsigned int i=0; i<=d; i++)
      c[i] = c[i]*xc[0];
  }
  else if (xd >= 0) {
    if (length() < xd+1) {
      Sacado::Handle<TaylorData> h(new TaylorData(xd));
      T* cc = h->coeff_;
      for (unsigned int i=0; i<=xd; i++)
	cc[i] = c[0]*xc[i];
      th = h;
    }
    else {
      th.makeOwnCopy();
      for (unsigned int i=d; i>=0; i--)
	c[i] = c[0]*xc[i];
    }
  }

  return *this;
}

template <typename T> 
Sacado::Tay::Taylor<T>& 
Sacado::Tay::Taylor<T>::
operator/=(const Taylor<T>& x)
{
  unsigned int d = degree();
  unsigned int xd = x.degree();

#ifdef SACADO_DEBUG
  if ((xd != d) && (xd != 0) && (d != 0))
    throw "Taylor Error:  Attempt to assign with incompatible degrees";
#endif
  
  T* c = th->coeff_; 
  T* xc = x.th->coeff_;
  
  if (d > 0 && xd > 0) {
    Sacado::Handle<TaylorData> h(new TaylorData(d));
    T* cc = h->coeff_;
    T tmp;
    for(unsigned int i=0; i<=d; i++) {
      tmp = c[i];
      for (unsigned int k=1; k<=i; ++k)
	tmp -= xc[k]*cc[i-k];
      cc[i] = tmp / xc[0];
    }
    th = h;
  }
  else if (d > 0) {
    th.makeOwnCopy();
    for (unsigned int i=0; i<=d; i++)
      c[i] = c[i]/xc[0];
  }
  else if (xd >= 0) {
    Sacado::Handle<TaylorData> h(new TaylorData(xd));
    T* cc = h->coeff_;
    T tmp;
    cc[0] = c[0] / xc[0];
    for (unsigned int i=1; i<=xd; i++) {
      tmp = T(0.0);
      for (unsigned int k=1; k<=i; ++k)
	tmp -= xc[k]*cc[i-k];
      cc[i] = tmp / xc[0];
    }
    th = h;
  }

  return *this;
}

template <typename T>
void
Sacado::Tay::Taylor<T>::
resizeCoeffs(unsigned int len)
{
  if (th->coeff_)
    Sacado::ds_array<T>::destroy_and_release(th->coeff_, th->len_);
  th->len_ = len;
  th->coeff_ = Sacado::ds_array<T>::get_and_fill(th->len_);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator+(const Sacado::Tay::Taylor<T>& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[i] + cb[i];
  }
  else if (da > 0) {
    cc[0] = ca[0] + cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = ca[i];
  }
  else if (db >= 0) {
    cc[0] = ca[0] + cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = cb[i];
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator+(const T& a, const Sacado::Tay::Taylor<T>& b)
{
  unsigned int dc = b.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a + cb[0];
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = cb[i];

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator+(const Sacado::Tay::Taylor<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = ca[i];

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator-(const Sacado::Tay::Taylor<T>& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[i] - cb[i];
  }
  else if (da > 0) {
    cc[0] = ca[0] - cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = ca[i];
  }
  else if (db >= 0) {
    cc[0] = ca[0] - cb[0];
    for (unsigned int i=1; i<=dc; i++)
      cc[i] = -cb[i];
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator-(const T& a, const Sacado::Tay::Taylor<T>& b)
{
  unsigned int dc = b.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a - cb[0];
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = -cb[i];

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator-(const Sacado::Tay::Taylor<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = ca[i];

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator*(const Sacado::Tay::Taylor<T>& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    T tmp;
    for (unsigned int i=0; i<=dc; i++) {
      tmp = T(0.0);
      for (unsigned int k=0; k<=i; k++)
	tmp += ca[k]*cb[i-k];
      cc[i] = tmp;
    }
  }
  else if (da > 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[i]*cb[0];
  }
  else if (db >= 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[0]*cb[i];
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator*(const T& a, const Sacado::Tay::Taylor<T>& b)
{
  unsigned int dc = b.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = a*cb[i];

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator*(const Sacado::Tay::Taylor<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = ca[i]*b;

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator/(const Sacado::Tay::Taylor<T>& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    T tmp;
    for (unsigned int i=0; i<=dc; i++) {
      tmp = ca[i];
      for (unsigned int k=0; k<=i; k++)
	tmp -= cb[k]*cc[i-k];
      cc[i] = tmp / cb[0];
    }
  }
  else if (da > 0) {
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = ca[i]/cb[0];
  }
  else if (db >= 0) {
    T tmp;
    cc[0] = ca[0] / cb[0];
    for (unsigned int i=1; i<=dc; i++) {
      tmp = T(0.0);
      for (unsigned int k=0; k<=i; k++)
	tmp -= cb[k]*cc[i-k];
      cc[i] = tmp / cb[0];
    }
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator/(const T& a, const Sacado::Tay::Taylor<T>& b)
{
  unsigned int dc = b.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = a / cb[0];
  for (unsigned int i=1; i<=dc; i++) {
    tmp = T(0.0);
    for (unsigned int k=0; k<=i; k++)
      tmp -= cb[k]*cc[i-k];
    cc[i] = tmp / cb[0];
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
operator/(const Sacado::Tay::Taylor<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = ca[i]/b;

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
exp(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = std::exp(ca[0]);
  for (unsigned int i=1; i<=dc; i++) {
    tmp = T(0.0);
    for (unsigned int k=1; k<=i; k++)
      tmp += k*cc[i-k]*ca[k];
    cc[i] = tmp / i;
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
log(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = std::log(ca[0]);
  for (unsigned int i=1; i<=dc; i++) {
    tmp = i*ca[i];
    for (unsigned int k=1; k<=i-1; k++)
      tmp -= k*ca[i-k]*cc[k];
    cc[i] = tmp / (i*ca[0]);
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
log10(const Sacado::Tay::Taylor<T>& a)
{
  return log(a) / std::log(10.0);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
sqrt(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = std::sqrt(ca[0]);
  for (unsigned int i=1; i<=dc; i++) {
    tmp = ca[i];
    for (unsigned int k=1; k<=i-1; k++)
      tmp -= cc[k]*cc[i-k];
    cc[i] = tmp / (2.0*cc[0]);
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
pow(const Sacado::Tay::Taylor<T>& a,
    const Sacado::Tay::Taylor<T>& b)
{
  return exp(b*log(a));
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
pow(const T& a,
    const Sacado::Tay::Taylor<T>& b)
{
  return exp(b*log(a));
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
pow(const Sacado::Tay::Taylor<T>& a,
    const T& b)
{
  return exp(b*log(a));
}

template <typename T>
void
Sacado::Tay::
sincos(const Sacado::Tay::Taylor<T>& a,
       Sacado::Tay::Taylor<T>& s,
       Sacado::Tay::Taylor<T>& c)
{
  unsigned int dc = a.degree();
  if (s.degree() != dc)
    s.resize(dc);
  if (c.degree() != dc)
    c.resize(dc);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  T tmp1;
  T tmp2;
  cs[0] = std::sin(ca[0]);
  cc[0] = std::cos(ca[0]);
  for (unsigned int i=1; i<=dc; i++) {
    tmp1 = T(0.0);
    tmp2 = T(0.0);
    for (unsigned int k=1; k<=i; k++) {
      tmp1 += k*ca[k]*cc[i-k];
      tmp2 -= k*ca[k]*cs[i-k];
    }
    cs[i] = tmp1 / i;
    cc[i] = tmp2 / i;
  }
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
sin(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();
  Sacado::Tay::Taylor<T> s(dc);
  Sacado::Tay::Taylor<T> c(dc);
  sincos(a, s, c);

  return s;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
cos(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();
  Sacado::Tay::Taylor<T> s(dc);
  Sacado::Tay::Taylor<T> c(dc);
  sincos(a, s, c);

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
tan(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();
  Sacado::Tay::Taylor<T> s(dc);
  Sacado::Tay::Taylor<T> c(dc);
  
  sincos(a, s, c);

  return s / c;
}

template <typename T>
void
Sacado::Tay::
sinhcosh(const Sacado::Tay::Taylor<T>& a,
	 Sacado::Tay::Taylor<T>& s,
	 Sacado::Tay::Taylor<T>& c)
{
  unsigned int dc = a.degree();
  if (s.degree() != dc)
    s.resize(dc);
  if (c.degree() != dc)
    c.resize(dc);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  T tmp1;
  T tmp2;
  cs[0] = std::sinh(ca[0]);
  cc[0] = std::cosh(ca[0]);
  for (unsigned int i=1; i<=dc; i++) {
    tmp1 = T(0.0);
    tmp2 = T(0.0);
    for (unsigned int k=1; k<=i; k++) {
      tmp1 += k*ca[k]*cc[i-k];
      tmp2 += k*ca[k]*cs[i-k];
    }
    cs[i] = tmp1 / i;
    cc[i] = tmp2 / i;
  }
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
sinh(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();
  Sacado::Tay::Taylor<T> s(dc);
  Sacado::Tay::Taylor<T> c(dc);
  sinhcosh(a, s, c);

  return s;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
cosh(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();
  Sacado::Tay::Taylor<T> s(dc);
  Sacado::Tay::Taylor<T> c(dc);
  sinhcosh(a, s, c);

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
tanh(const Sacado::Tay::Taylor<T>& a)
{
  unsigned int dc = a.degree();
  Sacado::Tay::Taylor<T> s(dc);
  Sacado::Tay::Taylor<T> c(dc);
  
  sinhcosh(a, s, c);

  return s / c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
quad(const T& c0,
     const Sacado::Tay::Taylor<T>& a,
     const Sacado::Tay::Taylor<T>& b)
{
  unsigned int dc = a.degree();

  Sacado::Tay::Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = c0;
  for (unsigned int i=1; i<=dc; i++) {
    tmp = T(0.0);
    for (unsigned int k=1; k<=i; k++) 
      tmp += k*ca[k]*cb[i-k];
    cc[i] = tmp / i;
  }

  return c;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
acos(const Sacado::Tay::Taylor<T>& a)
{
  Sacado::Tay::Taylor<T> b = -1.0 / sqrt(1.0 - a*a);
  return quad(std::acos(a.coeff(0)), a, b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
asin(const Sacado::Tay::Taylor<T>& a)
{
  Sacado::Tay::Taylor<T> b = 1.0 / sqrt(1.0 - a*a);
  return quad(std::asin(a.coeff(0)), a, b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
atan(const Sacado::Tay::Taylor<T>& a)
{
  Sacado::Tay::Taylor<T> b = 1.0 / (1.0 + a*a);
  return quad(std::atan(a.coeff(0)), a, b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
atan2(const Sacado::Tay::Taylor<T>& a,
      const Sacado::Tay::Taylor<T>& b)
{
  Sacado::Tay::Taylor<T> c = atan(a/b);
  c.fastAccessCoeff(0) = atan2(a.coeff(0),b.coeff(0));
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
atan2(const T& a,
      const Sacado::Tay::Taylor<T>& b)
{
  Sacado::Tay::Taylor<T> c = atan(a/b);
  c.fastAccessCoeff(0) = atan2(a,b.coeff(0));
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
atan2(const Sacado::Tay::Taylor<T>& a,
      const T& b)
{
  Sacado::Tay::Taylor<T> c = atan(a/b);
  c.fastAccessCoeff(0) = atan2(a.coeff(0),b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
acosh(const Sacado::Tay::Taylor<T>& a)
{
  Sacado::Tay::Taylor<T> b = -1.0 / sqrt(1.0 - a*a);
  return quad(std::acosh(a.coeff(0)), a, b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
asinh(const Sacado::Tay::Taylor<T>& a)
{
  Sacado::Tay::Taylor<T> b = 1.0 / sqrt(a*a - 1.0);
  return quad(std::asinh(a.coeff(0)), a, b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
atanh(const Sacado::Tay::Taylor<T>& a)
{
  Sacado::Tay::Taylor<T> b = 1.0 / (1.0 - a*a);
  return quad(std::atanh(a.coeff(0)), a, b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
fabs(const Sacado::Tay::Taylor<T>& a)
{
  if (a.coeff(0) >= 0)
    return a;
  else
    return -a;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
abs(const Sacado::Tay::Taylor<T>& a)
{
  if (a.coeff(0) >= 0)
    return a;
  else
    return -a;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
max(const Sacado::Tay::Taylor<T>& a,
    const Sacado::Tay::Taylor<T>& b)
{
  if (a.coeff(0) >= b.coeff(0))
    return a;
  else
    return b;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
max(const T& a,
    const Sacado::Tay::Taylor<T>& b)
{
  if (a >= b.coeff(0))
    return Taylor<T>(b.degree(), a);
  else
    return b;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
max(const Sacado::Tay::Taylor<T>& a,
    const T& b)
{
  if (a.coeff(0) >= b)
    return a;
  else
    return Taylor<T>(a.degree(), b);
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
min(const Sacado::Tay::Taylor<T>& a,
    const Sacado::Tay::Taylor<T>& b)
{
  if (a.coeff(0) <= b.coeff(0))
    return a;
  else
    return b;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
min(const T& a,
    const Sacado::Tay::Taylor<T>& b)
{
  if (a <= b.coeff(0))
    return Taylor<T>(b.degree(), a);
  else
    return b;
}

template <typename T>
Sacado::Tay::Taylor<T>
Sacado::Tay::
min(const Sacado::Tay::Taylor<T>& a,
    const T& b)
{
  if (a.coeff(0) <= b)
    return a;
  else
    return Taylor<T>(a.degree(), b);
}

template <typename T>
bool
Sacado::Tay::
operator==(const Sacado::Tay::Taylor<T>& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a.coeff(0) == b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator==(const T& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a == b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator==(const Sacado::Tay::Taylor<T>& a, 
	   const T& b)
{
  return a.coeff(0) == b;
}

template <typename T>
bool
Sacado::Tay::
operator!=(const Sacado::Tay::Taylor<T>& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a.coeff(0) != b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator!=(const T& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a != b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator!=(const Sacado::Tay::Taylor<T>& a, 
	   const T& b)
{
  return a.coeff(0) != b;
}

template <typename T>
bool
Sacado::Tay::
operator<=(const Sacado::Tay::Taylor<T>& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a.coeff(0) <= b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator<=(const T& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a <= b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator<=(const Sacado::Tay::Taylor<T>& a, 
	   const T& b)
{
  return a.coeff(0) <= b;
}

template <typename T>
bool
Sacado::Tay::
operator>=(const Sacado::Tay::Taylor<T>& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a.coeff(0) >= b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator>=(const T& a, 
	   const Sacado::Tay::Taylor<T>& b)
{
  return a >= b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator>=(const Sacado::Tay::Taylor<T>& a, 
	   const T& b)
{
  return a.coeff(0) >= b;
}

template <typename T>
bool
Sacado::Tay::
operator<(const Sacado::Tay::Taylor<T>& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  return a.coeff(0) < b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator<(const T& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  return a < b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator<(const Sacado::Tay::Taylor<T>& a, 
	  const T& b)
{
  return a.coeff(0) < b;
}

template <typename T>
bool
Sacado::Tay::
operator>(const Sacado::Tay::Taylor<T>& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  return a.coeff(0) > b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator>(const T& a, 
	  const Sacado::Tay::Taylor<T>& b)
{
  return a > b.coeff(0);
}

template <typename T>
bool
Sacado::Tay::
operator>(const Sacado::Tay::Taylor<T>& a, 
	  const T& b)
{
  return a.coeff(0) > b;
}

template <typename T>
std::ostream& 
Sacado::Tay::
operator << (std::ostream& os, const Sacado::Tay::Taylor<T>& a)
{
  os << "[ ";
      
  for (unsigned int i=0; i<=a.degree(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}
