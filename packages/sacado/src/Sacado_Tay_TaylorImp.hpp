// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

#include "Sacado_ConfigDefs.h"
#include "Sacado_DynamicArrayTraits.hpp"
#include <ostream>      // for std::ostream

namespace Sacado {
namespace Tay {

template <typename T>
Taylor<T>::TaylorData::
TaylorData() :
  coeff_(NULL), deg_(-1), len_(0)
{
}

template <typename T>
Taylor<T>::TaylorData::
TaylorData(const T& x) :
  coeff_(), deg_(0), len_(1)
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T>
Taylor<T>::TaylorData::
TaylorData(int d, const T& x) :
  coeff_(), deg_(d), len_(d+1)
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T>
Taylor<T>::TaylorData::
TaylorData(int d) :
  coeff_(), deg_(d), len_(d+1)
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
}

template <typename T>
Taylor<T>::TaylorData::
TaylorData(int d, int l) :
  coeff_(), deg_(d), len_(l)
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
}

template <typename T>
Taylor<T>::TaylorData::
TaylorData(const typename Taylor<T>::TaylorData& x) :
  coeff_(), deg_(x.deg_), len_(x.deg_+1)
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(x.coeff_, len_);
}

template <typename T>
Taylor<T>::TaylorData::
~TaylorData()
{
  if (len_ > 0)
    Sacado::ds_array<T>::destroy_and_release(coeff_, len_);
}

template <typename T>
typename Taylor<T>::TaylorData&
Taylor<T>::TaylorData::
operator=(const typename Taylor<T>::TaylorData& x)
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
Taylor<T>::
Taylor() :
  th(new TaylorData)
{
}

template <typename T>
Taylor<T>::
Taylor(const T& x) :
  th(new TaylorData(x))
{
}

template <typename T>
Taylor<T>::
Taylor(const typename dummy<value_type,scalar_type>::type& x) :
  th(new TaylorData(value_type(x)))
{
}

template <typename T>
Taylor<T>::
Taylor(int d, const T& x) :
  th(new TaylorData(d, x))
{
}

template <typename T>
Taylor<T>::
Taylor(int d) :
  th(new TaylorData(d))
{
}

template <typename T>
Taylor<T>::
Taylor(const Taylor<T>& x) :
  th(x.th)
{
}

template <typename T>
Taylor<T>::
~Taylor()
{
}

template <typename T>
void
Taylor<T>::
resize(int d, bool keep_coeffs)
{
  if (d+1 > length()) {
    Sacado::Handle<TaylorData> h(new TaylorData(d));
    if (keep_coeffs)
      Sacado::ds_array<T>::copy(th->coeff_, h->coeff_, th->deg_+1);
    th = h;
  }
  else {
    th.makeOwnCopy();
    if (!keep_coeffs)
      Sacado::ds_array<T>::zero(coeff(), degree()+1);
    th->deg_ = d;
  }
}

template <typename T>
void
Taylor<T>::
reserve(int d)
{
  if (d+1 > length()) {
    Sacado::Handle<TaylorData> h(new TaylorData(th->deg_,d+1));
    Sacado::ds_array<T>::copy(th->coeff_, h->coeff_, th->deg_+1);
    th = h;
  }
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator=(const T& v)
{
  th.makeOwnCopy();

  if (th->len_ == 0) {
    th->len_ = 1;
    th->deg_ = 0;
    th->coeff_ = Sacado::ds_array<T>::get_and_fill(th->len_);
  }

  th->coeff_[0] = v;
  Sacado::ds_array<T>::zero(th->coeff_+1, th->deg_);

  return *this;
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator=(const typename dummy<value_type,scalar_type>::type& v)
{
  return operator=(value_type(v));
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator=(const Taylor<T>& x)
{
  th = x.th;
  return *this;
}

template <typename T>
Taylor<T>
Taylor<T>::
operator+() const
{
  return *this;
}

template <typename T>
Taylor<T>
Taylor<T>::
operator-() const
{
  Taylor<T> x;
  x.th->deg_ = th->deg_;
  x.th->len_ = th->deg_+1;
  x.th->coeff_ = Sacado::ds_array<T>::get_and_fill(x.th->len_);
  for (int i=0; i<=th->deg_; i++)
    x.th->coeff_[i] = -th->coeff_[i];

  return x;
}

template <typename T>
 Taylor<T>&
Taylor<T>::
operator+=(const T& v)
{
  th.makeOwnCopy();

  th->coeff_[0] += v;

  return *this;
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator-=(const T& v)
{
  th.makeOwnCopy();

  th->coeff_[0] -= v;

  return *this;
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator*=(const T& v)
{
  th.makeOwnCopy();

  for (int i=0; i<=th->deg_; i++)
    th->coeff_[i] *= v;

  return *this;
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator/=(const T& v)
{
  th.makeOwnCopy();

  for (int i=0; i<=th->deg_; i++)
    th->coeff_[i] /= v;

  return *this;
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator+=(const Taylor<T>& x)
{
  th.makeOwnCopy();

  int d = degree();
  int xd = x.degree();
  int dmin = xd < d ? xd : d;

  int l = xd + 1;
  bool need_resize = l > length();
  T* c;
  if (need_resize) {
    c = Sacado::ds_array<T>::get_and_fill(l);
    for (int i=0; i<=d; i++)
      c[i] = fastAccessCoeff(i);
  }
  else
    c = th->coeff_;
  T* xc = x.th->coeff_;

  for (int i=0; i<=dmin; i++)
    c[i] += xc[i];
  if (th->deg_ < xd) {
    for (int i=d+1; i<=xd; i++)
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
Taylor<T>&
Taylor<T>::
operator-=(const Taylor<T>& x)
{
  th.makeOwnCopy();

  int d = degree();
  int xd = x.degree();
  int dmin = xd < d ? xd : d;

  int l = xd + 1;
  bool need_resize = l > length();
  T* c;
  if (need_resize) {
    c = Sacado::ds_array<T>::get_and_fill(l);
    for (int i=0; i<=d; i++)
      c[i] = fastAccessCoeff(i);
  }
  else
    c = th->coeff_;
  T* xc = x.th->coeff_;

  for (int i=0; i<=dmin; i++)
    c[i] -= xc[i];
  if (d < xd) {
    for (int i=d+1; i<=xd; i++)
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
Taylor<T>&
Taylor<T>::
operator*=(const Taylor<T>& x)
{
  int d = degree();
  int xd = x.degree();

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
    for (int i=0; i<=d; i++) {
      tmp = T(0.0);
      for (int k=0; k<=i; ++k)
        tmp += c[k]*xc[i-k];
      cc[i] = tmp;
    }
    th = h;
  }
  else if (d > 0) {
    th.makeOwnCopy();
    for (int i=0; i<=d; i++)
      c[i] = c[i]*xc[0];
  }
  else if (xd >= 0) {
    if (length() < xd+1) {
      Sacado::Handle<TaylorData> h(new TaylorData(xd));
      T* cc = h->coeff_;
      for (int i=0; i<=xd; i++)
        cc[i] = c[0]*xc[i];
      th = h;
    }
    else {
      th.makeOwnCopy();
      for (int i=d; i>=0; i--)
        c[i] = c[0]*xc[i];
    }
  }

  return *this;
}

template <typename T>
Taylor<T>&
Taylor<T>::
operator/=(const Taylor<T>& x)
{
  int d = degree();
  int xd = x.degree();

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
    for(int i=0; i<=d; i++) {
      tmp = c[i];
      for (int k=1; k<=i; ++k)
        tmp -= xc[k]*cc[i-k];
      cc[i] = tmp / xc[0];
    }
    th = h;
  }
  else if (d > 0) {
    th.makeOwnCopy();
    for (int i=0; i<=d; i++)
      c[i] = c[i]/xc[0];
  }
  else if (xd >= 0) {
    Sacado::Handle<TaylorData> h(new TaylorData(xd));
    T* cc = h->coeff_;
    T tmp;
    cc[0] = c[0] / xc[0];
    for (int i=1; i<=xd; i++) {
      tmp = T(0.0);
      for (int k=1; k<=i; ++k)
        tmp -= xc[k]*cc[i-k];
      cc[i] = tmp / xc[0];
    }
    th = h;
  }

  return *this;
}

template <typename T>
void
Taylor<T>::
resizeCoeffs(int len)
{
  if (th->coeff_)
    Sacado::ds_array<T>::destroy_and_release(th->coeff_, th->len_);
  th->len_ = len;
  th->coeff_ = Sacado::ds_array<T>::get_and_fill(th->len_);
}

template <typename T>
Taylor<T>
operator+(const Base< Taylor<T> >& aa,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();

  int da = a.degree();
  int db = b.degree();
  int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    for (int i=0; i<=dc; i++)
      cc[i] = ca[i] + cb[i];
  }
  else if (da > 0) {
    cc[0] = ca[0] + cb[0];
    for (int i=1; i<=dc; i++)
      cc[i] = ca[i];
  }
  else if (db >= 0) {
    cc[0] = ca[0] + cb[0];
    for (int i=1; i<=dc; i++)
      cc[i] = cb[i];
  }

  return c;
}

template <typename T>
Taylor<T>
operator+(const typename Taylor<T>::value_type& a,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();

  int dc = b.degree();

  Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a + cb[0];
  for (int i=1; i<=dc; i++)
    cc[i] = cb[i];

  return c;
}

template <typename T>
Taylor<T>
operator+(const Base< Taylor<T> >& aa,
          const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();

  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (int i=1; i<=dc; i++)
    cc[i] = ca[i];

  return c;
}

template <typename T>
Taylor<T>
operator-(const Base< Taylor<T> >& aa,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();

  int da = a.degree();
  int db = b.degree();
  int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    for (int i=0; i<=dc; i++)
      cc[i] = ca[i] - cb[i];
  }
  else if (da > 0) {
    cc[0] = ca[0] - cb[0];
    for (int i=1; i<=dc; i++)
      cc[i] = ca[i];
  }
  else if (db >= 0) {
    cc[0] = ca[0] - cb[0];
    for (int i=1; i<=dc; i++)
      cc[i] = -cb[i];
  }

  return c;
}

template <typename T>
Taylor<T>
operator-(const typename Taylor<T>::value_type& a,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();

  int dc = b.degree();

  Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a - cb[0];
  for (int i=1; i<=dc; i++)
    cc[i] = -cb[i];

  return c;
}

template <typename T>
Taylor<T>
operator-(const Base< Taylor<T> >& aa,
          const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();

  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (int i=1; i<=dc; i++)
    cc[i] = ca[i];

  return c;
}

template <typename T>
Taylor<T>
operator*(const Base< Taylor<T> >& aa,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();

  int da = a.degree();
  int db = b.degree();
  int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    T tmp;
    for (int i=0; i<=dc; i++) {
      tmp = T(0.0);
      for (int k=0; k<=i; k++)
        tmp += ca[k]*cb[i-k];
      cc[i] = tmp;
    }
  }
  else if (da > 0) {
    for (int i=0; i<=dc; i++)
      cc[i] = ca[i]*cb[0];
  }
  else if (db >= 0) {
    for (int i=0; i<=dc; i++)
      cc[i] = ca[0]*cb[i];
  }

  return c;
}

template <typename T>
Taylor<T>
operator*(const typename Taylor<T>::value_type& a,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();

  int dc = b.degree();

  Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  for (int i=0; i<=dc; i++)
    cc[i] = a*cb[i];

  return c;
}

template <typename T>
Taylor<T>
operator*(const Base< Taylor<T> >& aa,
          const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();

  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (int i=0; i<=dc; i++)
    cc[i] = ca[i]*b;

  return c;
}

template <typename T>
Taylor<T>
operator/(const Base< Taylor<T> >& aa,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();

  int da = a.degree();
  int db = b.degree();
  int dc = da > db ? da : db;

#ifdef SACADO_DEBUG
  if ((da != db) && (da != 0) && (db != 0))
    throw "operator+():  Arguments have incompatible degrees!";
#endif

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
    T tmp;
    for (int i=0; i<=dc; i++) {
      tmp = ca[i];
      for (int k=0; k<=i; k++)
        tmp -= cb[k]*cc[i-k];
      cc[i] = tmp / cb[0];
    }
  }
  else if (da > 0) {
    for (int i=0; i<=dc; i++)
      cc[i] = ca[i]/cb[0];
  }
  else if (db >= 0) {
    T tmp;
    cc[0] = ca[0] / cb[0];
    for (int i=1; i<=dc; i++) {
      tmp = T(0.0);
      for (int k=0; k<=i; k++)
        tmp -= cb[k]*cc[i-k];
      cc[i] = tmp / cb[0];
    }
  }

  return c;
}

template <typename T>
Taylor<T>
operator/(const typename Taylor<T>::value_type& a,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();

  int dc = b.degree();

  Taylor<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = a / cb[0];
  for (int i=1; i<=dc; i++) {
    tmp = T(0.0);
    for (int k=0; k<=i; k++)
      tmp -= cb[k]*cc[i-k];
    cc[i] = tmp / cb[0];
  }

  return c;
}

template <typename T>
Taylor<T>
operator/(const Base< Taylor<T> >& aa,
          const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();

  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (int i=0; i<=dc; i++)
    cc[i] = ca[i]/b;

  return c;
}

template <typename T>
Taylor<T>
exp(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = std::exp(ca[0]);
  for (int i=1; i<=dc; i++) {
    tmp = T(0.0);
    for (int k=1; k<=i; k++)
      tmp += k*cc[i-k]*ca[k];
    cc[i] = tmp / i;
  }

  return c;
}

template <typename T>
Taylor<T>
log(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = std::log(ca[0]);
  for (int i=1; i<=dc; i++) {
    tmp = i*ca[i];
    for (int k=1; k<=i-1; k++)
      tmp -= k*ca[i-k]*cc[k];
    cc[i] = tmp / (i*ca[0]);
  }

  return c;
}

template <typename T>
Taylor<T>
log10(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  return log(a) / std::log(10.0);
}

template <typename T>
Taylor<T>
sqrt(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = std::sqrt(ca[0]);
  for (int i=1; i<=dc; i++) {
    tmp = ca[i];
    for (int k=1; k<=i-1; k++)
      tmp -= cc[k]*cc[i-k];
    cc[i] = tmp / (2.0*cc[0]);
  }

  return c;
}

template <typename T>
Taylor<T>
cbrt(const Base< Taylor<T> >& aa)
{
  return pow(aa, typename Taylor<T>::value_type(1.0/3.0));
}

template <typename T>
Taylor<T>
pow(const Base< Taylor<T> >& aa,
    const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  return exp(b*log(a));
}

template <typename T>
Taylor<T>
pow(const typename Taylor<T>::value_type& a,
    const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();
  return exp(b*std::log(a));
}

template <typename T>
Taylor<T>
pow(const Base< Taylor<T> >& aa,
    const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();
  return exp(b*log(a));
}

template <typename T>
void
sincos(const Base< Taylor<T> >& aa,
       Taylor<T>& s,
       Taylor<T>& c)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  if (s.degree() != dc)
    s.resize(dc, false);
  if (c.degree() != dc)
    c.resize(dc, false);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  T tmp1;
  T tmp2;
  cs[0] = std::sin(ca[0]);
  cc[0] = std::cos(ca[0]);
  for (int i=1; i<=dc; i++) {
    tmp1 = T(0.0);
    tmp2 = T(0.0);
    for (int k=1; k<=i; k++) {
      tmp1 += k*ca[k]*cc[i-k];
      tmp2 -= k*ca[k]*cs[i-k];
    }
    cs[i] = tmp1 / i;
    cc[i] = tmp2 / i;
  }
}

template <typename T>
Taylor<T>
sin(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  Taylor<T> s(dc);
  Taylor<T> c(dc);
  sincos(a, s, c);

  return s;
}

template <typename T>
Taylor<T>
cos(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  Taylor<T> s(dc);
  Taylor<T> c(dc);
  sincos(a, s, c);

  return c;
}

template <typename T>
Taylor<T>
tan(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  Taylor<T> s(dc);
  Taylor<T> c(dc);

  sincos(a, s, c);

  return s / c;
}

template <typename T>
void
sinhcosh(const Base< Taylor<T> >& aa,
         Taylor<T>& s,
         Taylor<T>& c)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  if (s.degree() != dc)
    s.resize(dc, false);
  if (c.degree() != dc)
    c.resize(dc, false);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  T tmp1;
  T tmp2;
  cs[0] = std::sinh(ca[0]);
  cc[0] = std::cosh(ca[0]);
  for (int i=1; i<=dc; i++) {
    tmp1 = T(0.0);
    tmp2 = T(0.0);
    for (int k=1; k<=i; k++) {
      tmp1 += k*ca[k]*cc[i-k];
      tmp2 += k*ca[k]*cs[i-k];
    }
    cs[i] = tmp1 / i;
    cc[i] = tmp2 / i;
  }
}

template <typename T>
Taylor<T>
sinh(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  Taylor<T> s(dc);
  Taylor<T> c(dc);
  sinhcosh(a, s, c);

  return s;
}

template <typename T>
Taylor<T>
cosh(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  Taylor<T> s(dc);
  Taylor<T> c(dc);
  sinhcosh(a, s, c);

  return c;
}

template <typename T>
Taylor<T>
tanh(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  int dc = a.degree();
  Taylor<T> s(dc);
  Taylor<T> c(dc);

  sinhcosh(a, s, c);

  return s / c;
}

template <typename T>
Taylor<T>
quad(const typename Taylor<T>::value_type& c0,
     const Base< Taylor<T> >& aa,
     const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  int dc = a.degree();

  Taylor<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  T tmp;
  cc[0] = c0;
  for (int i=1; i<=dc; i++) {
    tmp = T(0.0);
    for (int k=1; k<=i; k++)
      tmp += k*ca[k]*cb[i-k];
    cc[i] = tmp / i;
  }

  return c;
}

template <typename T>
Taylor<T>
acos(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  Taylor<T> b = -1.0 / sqrt(1.0 - a*a);
  return quad(std::acos(a.coeff(0)), a, b);
}

template <typename T>
Taylor<T>
asin(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  Taylor<T> b = 1.0 / sqrt(1.0 - a*a);
  return quad(std::asin(a.coeff(0)), a, b);
}

template <typename T>
Taylor<T>
atan(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  Taylor<T> b = 1.0 / (1.0 + a*a);
  return quad(std::atan(a.coeff(0)), a, b);
}

template <typename T>
Taylor<T>
atan2(const Base< Taylor<T> >& aa,
      const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();

  Taylor<T> c = atan(a/b);
  c.fastAccessCoeff(0) = atan2(a.coeff(0),b.coeff(0));
  return c;
}

template <typename T>
Taylor<T>
atan2(const typename Taylor<T>::value_type& a,
      const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();

  Taylor<T> c = atan(a/b);
  c.fastAccessCoeff(0) = atan2(a,b.coeff(0));
  return c;
}

template <typename T>
Taylor<T>
atan2(const Base< Taylor<T> >& aa,
      const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();

  Taylor<T> c = atan(a/b);
  c.fastAccessCoeff(0) = atan2(a.coeff(0),b);
  return c;
}

template <typename T>
Taylor<T>
acosh(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  Taylor<T> b = -1.0 / sqrt(1.0 - a*a);
  return quad(acosh(a.coeff(0)), a, b);
}

template <typename T>
Taylor<T>
asinh(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  Taylor<T> b = 1.0 / sqrt(a*a - 1.0);
  return quad(asinh(a.coeff(0)), a, b);
}

template <typename T>
Taylor<T>
atanh(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  Taylor<T> b = 1.0 / (1.0 - a*a);
  return quad(atanh(a.coeff(0)), a, b);
}

template <typename T>
Taylor<T>
fabs(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  if (a.coeff(0) >= 0)
    return a;
  else
    return -a;
}

template <typename T>
Taylor<T>
abs(const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  if (a.coeff(0) >= 0)
    return a;
  else
    return -a;
}

template <typename T>
Taylor<T>
max(const Base< Taylor<T> >& aa,
    const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();

  if (a.coeff(0) >= b.coeff(0))
    return a;
  else
    return b;
}

template <typename T>
Taylor<T>
max(const typename Taylor<T>::value_type& a,
    const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();

  if (a >= b.coeff(0))
    return Taylor<T>(b.degree(), a);
  else
    return b;
}

template <typename T>
Taylor<T>
max(const Base< Taylor<T> >& aa,
    const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();

  if (a.coeff(0) >= b)
    return a;
  else
    return Taylor<T>(a.degree(), b);
}

template <typename T>
Taylor<T>
min(const Base< Taylor<T> >& aa,
    const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();

  if (a.coeff(0) <= b.coeff(0))
    return a;
  else
    return b;
}

template <typename T>
Taylor<T>
min(const typename Taylor<T>::value_type& a,
    const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();

  if (a <= b.coeff(0))
    return Taylor<T>(b.degree(), a);
  else
    return b;
}

template <typename T>
Taylor<T>
min(const Base< Taylor<T> >& aa,
    const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();

  if (a.coeff(0) <= b)
    return a;
  else
    return Taylor<T>(a.degree(), b);
}

template <typename T>
bool
operator==(const Base< Taylor<T> >& aa,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  return a.coeff(0) == b.coeff(0);
}

template <typename T>
bool
operator==(const typename Taylor<T>::value_type& a,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();
  return a == b.coeff(0);
}

template <typename T>
bool
operator==(const Base< Taylor<T> >& aa,
           const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();
  return a.coeff(0) == b;
}

template <typename T>
bool
operator!=(const Base< Taylor<T> >& aa,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  return a.coeff(0) != b.coeff(0);
}

template <typename T>
bool
operator!=(const typename Taylor<T>::value_type& a,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();
  return a != b.coeff(0);
}

template <typename T>
bool
operator!=(const Base< Taylor<T> >& aa,
           const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();
  return a.coeff(0) != b;
}

template <typename T>
bool
operator<=(const Base< Taylor<T> >& aa,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  return a.coeff(0) <= b.coeff(0);
}

template <typename T>
bool
operator<=(const typename Taylor<T>::value_type& a,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();
  return a <= b.coeff(0);
}

template <typename T>
bool
operator<=(const Base< Taylor<T> >& aa,
           const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();
  return a.coeff(0) <= b;
}

template <typename T>
bool
operator>=(const Base< Taylor<T> >& aa,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  return a.coeff(0) >= b.coeff(0);
}

template <typename T>
bool
operator>=(const typename Taylor<T>::value_type& a,
           const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();
  return a >= b.coeff(0);
}

template <typename T>
bool
operator>=(const Base< Taylor<T> >& aa,
           const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();
  return a.coeff(0) >= b;
}

template <typename T>
bool
operator<(const Base< Taylor<T> >& aa,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  return a.coeff(0) < b.coeff(0);
}

template <typename T>
bool
operator<(const typename Taylor<T>::value_type& a,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();
  return a < b.coeff(0);
}

template <typename T>
bool
operator<(const Base< Taylor<T> >& aa,
          const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();
  return a.coeff(0) < b;
}

template <typename T>
bool
operator>(const Base< Taylor<T> >& aa,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& a = aa.derived();
  const Taylor<T>& b = bb.derived();
  return a.coeff(0) > b.coeff(0);
}

template <typename T>
bool
operator>(const typename Taylor<T>::value_type& a,
          const Base< Taylor<T> >& bb)
{
  const Taylor<T>& b = bb.derived();
  return a > b.coeff(0);
}

template <typename T>
bool
operator>(const Base< Taylor<T> >& aa,
          const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& a = aa.derived();
  return a.coeff(0) > b;
}

template <typename T>
bool toBool(const Taylor<T>& x) {
  bool is_zero = true;
  for (int i=0; i<=x.degree(); i++)
    is_zero = is_zero && (x.coeff(i) == 0.0);
  return !is_zero;
}

template <typename T>
inline bool
operator && (const Base< Taylor<T> >& xx1, const Base< Taylor<T> >& xx2)
{
  const Taylor<T>& x1 = xx1.derived();
  const Taylor<T>& x2 = xx2.derived();
  return toBool(x1) && toBool(x2);
}

template <typename T>
inline bool
operator && (const typename Taylor<T>::value_type& a,
             const Base< Taylor<T> >& xx2)
{
  const Taylor<T>& x2 = xx2.derived();
  return a && toBool(x2);
}

template <typename T>
inline bool
operator && (const Base< Taylor<T> >& xx1,
             const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& x1 = xx1.derived();
  return toBool(x1) && b;
}

template <typename T>
inline bool
operator || (const Base< Taylor<T> >& xx1, const Base< Taylor<T> >& xx2)
{
  const Taylor<T>& x1 = xx1.derived();
  const Taylor<T>& x2 = xx2.derived();
  return toBool(x1) || toBool(x2);
}

template <typename T>
inline bool
operator || (const typename Taylor<T>::value_type& a,
             const Base< Taylor<T> >& xx2)
{
  const Taylor<T>& x2 = xx2.derived();
  return a || toBool(x2);
}

template <typename T>
inline bool
operator || (const Base< Taylor<T> >& xx1,
             const typename Taylor<T>::value_type& b)
{
  const Taylor<T>& x1 = xx1.derived();
  return toBool(x1) || b;
}

template <typename T>
std::ostream&
operator << (std::ostream& os, const Base< Taylor<T> >& aa)
{
  const Taylor<T>& a = aa.derived();
  os << "[ ";

  for (int i=0; i<=a.degree(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]";
  return os;
}

} // namespace Tay
} // namespace Sacado
