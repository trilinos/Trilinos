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

#include "Teuchos_TestForException.hpp"
#include "Sacado_ConfigDefs.h"
#include "Sacado_DynamicArrayTraits.hpp"
#include <ostream>	// for std::ostream

namespace Sacado {
namespace PCE {

// Initialize static data
template <typename T> typename Hermite<T>::ws_type Hermite<T>::workspace(1);

template <typename T> 
Hermite<T>::HermiteData::
HermiteData() :
  coeff_(NULL), deg_(-1), len_(0) 
{
}

template <typename T> 
Hermite<T>::HermiteData::
HermiteData(const T& x) :
  coeff_(), deg_(0), len_(1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T> 
Hermite<T>::HermiteData::
HermiteData(unsigned int d, const T& x) :
  coeff_(), deg_(d), len_(d+1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
  coeff_[0] = x;
}

template <typename T> 
Hermite<T>::HermiteData::
HermiteData(unsigned int d) :
  coeff_(), deg_(d), len_(d+1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
}

template <typename T> 
Hermite<T>::HermiteData::
HermiteData(unsigned int d, unsigned int l) :
  coeff_(), deg_(d), len_(l) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(len_);
}

template <typename T> 
Hermite<T>::HermiteData::
HermiteData(const typename Hermite<T>::HermiteData& x) :
  coeff_(), deg_(x.deg_), len_(x.deg_+1) 
{
  coeff_ = Sacado::ds_array<T>::get_and_fill(x.coeff_, len_);
}

template <typename T> 
Hermite<T>::HermiteData::
~HermiteData()
{
  if (len_ > 0)
    Sacado::ds_array<T>::destroy_and_release(coeff_, len_);
}

template <typename T> 
typename Hermite<T>::HermiteData&
Hermite<T>::HermiteData::
operator=(const typename Hermite<T>::HermiteData& x) 
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
Hermite<T>::
Hermite() :
  th(new HermiteData)
{
}

template <typename T> 
Hermite<T>::
Hermite(const T& x) :
  th(new HermiteData(x))
{
}

template <typename T> 
Hermite<T>::
Hermite(unsigned int d, const T& x) :
  th(new HermiteData(d, x))
{
}

template <typename T> 
Hermite<T>::
Hermite(unsigned int d) :
  th(new HermiteData(d))
{
}

template <typename T> 
Hermite<T>::
Hermite(const Hermite<T>& x) :
  th(x.th)
{
}

template <typename T> 
Hermite<T>::
~Hermite()
{
}

template <typename T> 
void
Hermite<T>::
resize(unsigned int d, bool keep_coeffs)
{
  if (d+1 > length()) {
    Sacado::Handle<HermiteData> h(new HermiteData(d));
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
Hermite<T>::
reserve(unsigned int d)
{
  if (d+1 > length()) {
    Sacado::Handle<HermiteData> h(new HermiteData(th->deg_,d+1));
    Sacado::ds_array<T>::copy(th->coeff_, h->coeff_, th->deg_+1);
    th = h;
  }
}

template <typename T> 
void
Hermite<T>::
initWorkspace(unsigned int d)
{
  workspace.resize(d+1);
}

template <typename T> 
StandardPoly<T>
Hermite<T>::
toStandardBasis() const
{
  const typename Hermite<T>::ws_type::tp_type& Cijk = 
    Hermite<T>::workspace.getTripleProduct();
  const typename Hermite<T>::ws_type::tp_type::basis_type& basis = 
    Cijk.getBasis();
  return basis.toStandardBasis(th->coeff_, th->deg_+1);
}

template <typename T> 
Hermite<T>& 
Hermite<T>::
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
Hermite<T>& 
Hermite<T>::
operator=(const Hermite<T>& x) 
{
  th = x.th;
  return *this;
}

template <typename T> 
Hermite<T>
Hermite<T>::
operator+() const
{
  return *this;
}

template <typename T> 
Hermite<T> 
Hermite<T>::
operator-() const
{
  Hermite<T> x;
  x.th->deg_ = th->deg_;
  x.th->len_ = th->deg_+1;
  x.th->coeff_ = Sacado::ds_array<T>::get_and_fill(x.th->len_);
  for (unsigned int i=0; i<=th->deg_; i++)
    x.th->coeff_[i] = -th->coeff_[i];

  return x;
}

template <typename T> 
 Hermite<T>& 
Hermite<T>::
operator+=(const T& v)
{
  th.makeOwnCopy();

  th->coeff_[0] += v;

  return *this;
}

template <typename T> 
Hermite<T>& 
Hermite<T>::
operator-=(const T& v)
{
  th.makeOwnCopy();

  th->coeff_[0] -= v;

  return *this;
}

template <typename T> 
Hermite<T>& 
Hermite<T>::
operator*=(const T& v)
{
  th.makeOwnCopy();

  for (unsigned int i=0; i<=th->deg_; i++)
    th->coeff_[i] *= v;

  return *this;
}

template <typename T> 
Hermite<T>& 
Hermite<T>::
operator/=(const T& v)
{
  th.makeOwnCopy();

  for (unsigned int i=0; i<=th->deg_; i++)
    th->coeff_[i] /= v;

  return *this;
}

template <typename T> 
Hermite<T>& 
Hermite<T>::
operator+=(const Hermite<T>& x)
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
Hermite<T>& 
Hermite<T>::
operator-=(const Hermite<T>& x)
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
Hermite<T>& 
Hermite<T>::
operator*=(const Hermite<T>& x)
{
#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::operator*=()";
  TEST_FOR_EXCEPTION((x.degree() != degree()) && (x.degree() != 0) && 
		     (degree() != 0), std::logic_error,
		     func << ":  Attempt to assign with incompatible degrees");
#endif
  unsigned int d = degree();
  unsigned int xd = x.degree();
  
  T* c = th->coeff_; 
  T* xc = x.th->coeff_;
  
  if (d > 0 && xd > 0) {
#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < d+1, std::logic_error,
		       func << ":  Workspace size (" 
		       << Hermite<T>::workspace.size() 
		       << ") is too small for computation (" << d+1 
		       << " needed).");
#endif
    Sacado::Handle<HermiteData> h(new HermiteData(d));
    T* cc = h->coeff_;
    T tmp;
    const typename Hermite<T>::ws_type::tp_type& Cijk = 
      workspace.getTripleProduct();

    for (unsigned int k=0; k<=d; k++) {
      tmp = T(0.0);
      for (unsigned int i=0; i<=d; i++) {
	for (unsigned int j=0; j<=xd; j++)
	  tmp += Cijk.value(i,j,k)*c[i]*xc[j];
      }
      cc[k] = tmp / Cijk.norm_squared(k);
    }
    th = h;
  }
  else if (d > 0) {
    th.makeOwnCopy();
    for (unsigned int i=0; i<=d; i++)
      c[i] = c[i]*xc[0];
  }
  else  {
    if (length() < xd+1) {
      Sacado::Handle<HermiteData> h(new HermiteData(xd));
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
Hermite<T>& 
Hermite<T>::
operator/=(const Hermite<T>& x)
{
#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::operator/=()";
  TEST_FOR_EXCEPTION((x.degree() != degree()) && (x.degree() != 0) && 
		     (degree() != 0), std::logic_error,
		     func << ":  Attempt to assign with incompatible degrees");
#endif

  unsigned int d = degree();
  unsigned int xd = x.degree();
  
  T* c = th->coeff_; 
  T* xc = x.th->coeff_;
  
  if (xd > 0) {
    unsigned int dd = (xd > d) ? xd : d;

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dd+1, std::logic_error,
		       func << ":  Workspace size (" << 
		       Hermite<T>::workspace.size() 
		       << ") is too small for computation (" << dd+1 
		       << " needed).");
#endif

    Sacado::Handle<HermiteData> h(new HermiteData(dd));
    T* cc = h->coeff_;

    const typename Hermite<T>::ws_type::matrix_type& A = workspace.getMatrix();
    const typename Hermite<T>::ws_type::matrix_type& b = workspace.getRHS();
    const typename Hermite<T>::ws_type::tp_type& Cijk = 
      workspace.getTripleProduct();
    
    // Fill A
    for (unsigned int i=0; i<=dd; i++) {
      for (unsigned int j=0; j<=dd; j++) {
	A(i,j) = 0.0;
	for (unsigned int k=0; k<=xd; k++) {
	  A(i,j) += Cijk.value(i,j,k)*xc[k];
	}
	A(i,j) /= Cijk.norm_squared(i);
      }
    }
    
    // Fill b
    for (unsigned int i=0; i<=d; i++)
      b(i,0) = c[i];
    for (unsigned int i=d+1; i<=dd; i++)
      b(i,0) = T(0.0);

    // Solve system
    int info = Hermite<T>::workspace.solve(dd+1, 1);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

    // Get coefficients
    for (unsigned int i=0; i<=dd; i++)
      cc[i] = b(i,0);

    th = h;
  }
  else {
    th.makeOwnCopy();
    for (unsigned int i=0; i<=d; i++)
      c[i] = c[i]/xc[0];
  }

  return *this;
}

template <typename T>
void
Hermite<T>::
resizeCoeffs(unsigned int len)
{
  if (th->coeff_)
    Sacado::ds_array<T>::destroy_and_release(th->coeff_, th->len_);
  th->len_ = len;
  th->coeff_ = Sacado::ds_array<T>::get_and_fill(th->len_);
}

template <typename T>
Hermite<T>
operator+(const Hermite<T>& a, 
	  const Hermite<T>& b)
{
#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::operator+()";
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

  Hermite<T> c(dc);
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
Hermite<T>
operator+(const T& a, const Hermite<T>& b)
{
  unsigned int dc = b.degree();

  Hermite<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a + cb[0];
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = cb[i];

  return c;
}

template <typename T>
Hermite<T>
operator+(const Hermite<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] + b;
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = ca[i];

  return c;
}

template <typename T>
Hermite<T>
operator-(const Hermite<T>& a, 
	  const Hermite<T>& b)
{
#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::operator-()";
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

  Hermite<T> c(dc);
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
Hermite<T>
operator-(const T& a, const Hermite<T>& b)
{
  unsigned int dc = b.degree();

  Hermite<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  cc[0] = a - cb[0];
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = -cb[i];

  return c;
}

template <typename T>
Hermite<T>
operator-(const Hermite<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  cc[0] = ca[0] - b;
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = ca[i];

  return c;
}

template <typename T>
Hermite<T>
operator*(const Hermite<T>& a, 
	  const Hermite<T>& b)
{
#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::operator*()";
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (da > 0 && db > 0) {
#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
		       func << ":  Workspace size (" 
		       << Hermite<T>::workspace.size() 
		       << ") is too small for computation (" << dc+1 
		       << " needed).");
#endif

    T tmp;
    const typename Hermite<T>::ws_type::tp_type& Cijk = 
      Hermite<T>::workspace.getTripleProduct();

    for (unsigned int k=0; k<=dc; k++) {
      tmp = T(0.0);
      for (unsigned int i=0; i<=da; i++) {
	for (unsigned int j=0; j<=db; j++)
	  tmp += Cijk.value(i,j,k)*ca[i]*cb[j];
      }
      cc[k] = tmp / Cijk.norm_squared(k);
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
Hermite<T>
operator*(const T& a, const Hermite<T>& b)
{
  unsigned int dc = b.degree();

  Hermite<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = a*cb[i];

  return c;
}

template <typename T>
Hermite<T>
operator*(const Hermite<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = ca[i]*b;

  return c;
}

template <typename T>
Hermite<T>
operator/(const Hermite<T>& a, 
	  const Hermite<T>& b)
{
#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::operator/()";
  TEST_FOR_EXCEPTION((a.degree() != b.degree()) && (a.degree() != 0) && 
		     (b.degree() != 0), 
		     std::logic_error,
		     func << ":  Arguments have incompatible degrees!");
#endif

  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (db > 0) {

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
		       func << ":  Workspace size (" 
		       << Hermite<T>::workspace.size() 
		       << ") is too small for computation (" << dc+1 
		       << " needed).");
#endif

    typename Hermite<T>::ws_type::matrix_type& A = 
      Hermite<T>::workspace.getMatrix();
    typename Hermite<T>::ws_type::matrix_type& B = 
      Hermite<T>::workspace.getRHS();
    const typename Hermite<T>::ws_type::tp_type& Cijk = 
      Hermite<T>::workspace.getTripleProduct();
    
    // Fill A
    for (unsigned int i=0; i<=dc; i++) {
      for (unsigned int j=0; j<=dc; j++) {
	A(i,j) = 0.0;
	for (unsigned int k=0; k<=db; k++) {
	  A(i,j) += Cijk.value(i,j,k)*cb[k];
	}
	A(i,j) /= Cijk.norm_squared(i);
      }
    }
    
    // Fill b
    for (unsigned int i=0; i<=da; i++)
      B(i,0) = ca[i];
    for (unsigned int i=da+1; i<=dc; i++)
      B(i,0) = T(0.0);

    // Solve system
    int info = Hermite<T>::workspace.solve(dc+1, 1);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

    // Get coefficients
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = B(i,0);
  }
  else {
    for (unsigned int i=0; i<=da; i++)
      cc[i] = ca[i]/cb[0];
  }

  return c;
}

template <typename T>
Hermite<T>
operator/(const T& a, const Hermite<T>& b)
{
  unsigned int dc = b.degree();

  Hermite<T> c(dc);
  const T* cb = b.coeff();
  T* cc = c.coeff();

  if (dc > 0) {
#ifdef SACADO_DEBUG
    const char* func = "Sacado::PCE::Hermite::operator/()";
    TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
		       func << ":  Workspace size (" 
		       << Hermite<T>::workspace.size() 
		       << ") is too small for computation (" << dc+1 
		       << " needed).");
#endif

    typename Hermite<T>::ws_type::matrix_type& A = 
      Hermite<T>::workspace.getMatrix();
    typename Hermite<T>::ws_type::matrix_type& B = 
      Hermite<T>::workspace.getRHS();
    const typename Hermite<T>::ws_type::tp_type& Cijk = 
      Hermite<T>::workspace.getTripleProduct();
    
    // Fill A
    for (unsigned int i=0; i<=dc; i++) {
      for (unsigned int j=0; j<=dc; j++) {
	A(i,j) = 0.0;
	for (unsigned int k=0; k<=dc; k++) {
	  A(i,j) += Cijk.value(i,j,k)*cb[k];
	}
	A(i,j) /= Cijk.norm_squared(i);
      }
    }

    // Fill b
    B(0,0) = a;
    for (unsigned int i=1; i<=dc; i++)
      B(i,0) = T(0.0);

    // Solve system
    int info = Hermite<T>::workspace.solve(dc+1, 1);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

    // Get coefficients
    for (unsigned int i=0; i<=dc; i++)
      cc[i] = B(i,0);
  }
  else 
    cc[0] = a / cb[0];

  return c;
}

template <typename T>
Hermite<T>
operator/(const Hermite<T>& a, const T& b)
{
  unsigned int dc = a.degree();

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  for (unsigned int i=0; i<=dc; i++)
    cc[i] = ca[i]/b;

  return c;
}

template <typename T>
Hermite<T>
exp(const Hermite<T>& a)
{
  unsigned int dc = a.degree();

#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::exp()";
  TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
		     func << ":  Workspace size (" 
		     << Hermite<T>::workspace.size() 
		     << ") is too small for computation (" << dc+1 
		     << " needed).");
#endif

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  typename Hermite<T>::ws_type::matrix_type& A = 
    Hermite<T>::workspace.getMatrix();
  typename Hermite<T>::ws_type::matrix_type& b = 
    Hermite<T>::workspace.getRHS();
  const typename Hermite<T>::ws_type::tp_type& Cijk = 
    Hermite<T>::workspace.getTripleProduct();
  const typename Hermite<T>::ws_type::tp_type::basis_type& basis = 
    Cijk.getBasis();

  // Fill A and b
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      A(i-1,j-1) = 0.0;
      for (unsigned int k=1; k<=dc; k++) {
	if (i+j >= k) {
	  A(i-1,j-1) += Cijk.value(i-1,j,k-1)*basis.derivCoeff(k)*ca[k];
	}
      }
      A(i-1,j-1) /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
    }
    A(i-1,i-1) -= T(1.0);
    b(i-1,0) = -ca[i];
  }

  // Solve system
  int info = Hermite<T>::workspace.solve(dc, 1);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

  // Compute degree-0 coefficient
  T s = basis.getBasisPoly(0).coeff(0) * ca[0];
  T t = T(1.0);
  for (unsigned int i=1; i<=dc; i++) {
    s += basis.getBasisPoly(i).coeff(0) * ca[i];
    t += basis.getBasisPoly(i).coeff(0) * b(i-1,0);
  }
  s = std::exp(s);
  cc[0] = (s/t);

  // Compute remaining coefficients
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = b(i-1,0) * cc[0];

  cc[0] /= basis.getBasisPoly(0).coeff(0);

  return c;
}

template <typename T>
Hermite<T>
log(const Hermite<T>& a)
{
  unsigned int dc = a.degree();

#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::exp()";
  TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
		     func << ":  Workspace size (" 
		     << Hermite<T>::workspace.size() 
		     << ") is too small for computation (" << dc+1 
		     << " needed).");
#endif

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  T* cc = c.coeff();

  typename Hermite<T>::ws_type::matrix_type& A = 
    Hermite<T>::workspace.getMatrix();
  typename Hermite<T>::ws_type::matrix_type& b = 
    Hermite<T>::workspace.getRHS();
  const typename Hermite<T>::ws_type::tp_type& Cijk = 
    Hermite<T>::workspace.getTripleProduct();
  const typename Hermite<T>::ws_type::tp_type::basis_type& basis = 
    Cijk.getBasis();

  // Fill A and b
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      A(i-1,j-1) = 0.0;
      for (unsigned int k=0; k<=dc; k++) {
	if (i+j-2 >= k) {
	  A(i-1,j-1) += Cijk.value(i-1,j-1,k)*basis.derivCoeff(j)*ca[k];
	}
      }
      A(i-1,j-1) /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
    }
    b(i-1,0) = ca[i];
  }

  // Solve system
  int info = Hermite<T>::workspace.solve(dc, 1);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

  // Compute degree-0 coefficient
  T s = basis.getBasisPoly(0).coeff(0) * ca[0];
  T t = T(0.0);
  for (unsigned int i=1; i<=dc; i++) {
    s += basis.getBasisPoly(i).coeff(0) * ca[i];
    t += basis.getBasisPoly(i).coeff(0) * b(i-1,0);
  }
  cc[0] = (std::log(s) - t) / basis.getBasisPoly(0).coeff(0);

  // Compute remaining coefficients
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = b(i-1,0);

  return c;
}

template <typename T>
Hermite<T>
log10(const Hermite<T>& a)
{
  return log(a) / std::log(10.0);
}

template <typename T>
Hermite<T>
sqrt(const Hermite<T>& a)
{
  return exp(T(0.5)*log(a));
}

template <typename T>
Hermite<T>
pow(const Hermite<T>& a,
    const Hermite<T>& b)
{
  return exp(b*log(a));
}

template <typename T>
Hermite<T>
pow(const T& a,
    const Hermite<T>& b)
{
  return exp(b*std::log(a));
}

template <typename T>
Hermite<T>
pow(const Hermite<T>& a,
    const T& b)
{
  return exp(b*log(a));
}

template <typename T>
void
sincos(const Hermite<T>& a,
       Hermite<T>& s,
       Hermite<T>& c)
{
  unsigned int dc = a.degree();

#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::sincos()";
  TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
		     func << ":  Workspace size (" 
		     << Hermite<T>::workspace.size() 
		     << ") is too small for computation (" << dc+1 
		     << " needed).");
#endif

  if (s.degree() != dc)
    s.resize(dc, false);
  if (c.degree() != dc)
    c.resize(dc, false);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  typename Hermite<T>::ws_type::matrix_type& A = 
    Hermite<T>::workspace.getMatrix();
  typename Hermite<T>::ws_type::matrix_type& b = 
    Hermite<T>::workspace.getRHS();
  const typename Hermite<T>::ws_type::tp_type& Cijk = 
    Hermite<T>::workspace.getTripleProduct();
  const typename Hermite<T>::ws_type::tp_type::basis_type& basis = 
    Cijk.getBasis();

  // Fill A and b
  A.putScalar(T(0.0));
  b.putScalar(T(0.0));
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      T tmp = T(0.0);
      for (unsigned int k=1; k<=dc; k++) {
	if (i+j >= k) {
	  tmp += Cijk.value(i-1,j,k-1)*basis.derivCoeff(k)*ca[k];
	}
      }
      tmp /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
      A(i-1,j-1+dc) = -tmp;
      A(i-1+dc,j-1) = tmp;
    }
    A(i-1,i-1) = T(1.0);
    A(i-1+dc,i-1+dc) = T(1.0);
    b(i-1,0) = ca[i];
    b(i-1+dc,1) = ca[i];
  }

  // Solve system
  int info = Hermite<T>::workspace.solve(2*dc, 2);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

  // Compute degree-0 coefficients
  T t = basis.getBasisPoly(0).coeff(0) * ca[0];
  T a00 = T(0.0);
  T a01 = T(0.0);
  T a10 = T(0.0);
  T a11 = T(0.0);
  T b0 = b(0,0);
  T b1 = b(1,0);
  for (unsigned int i=1; i<=dc; i++) {
    t += basis.getBasisPoly(i).coeff(0) * ca[i];
    a00 += basis.getBasisPoly(i).coeff(0) * b(i-1,1);
    a01 += basis.getBasisPoly(i).coeff(0) * b(i-1,0);
    a10 += basis.getBasisPoly(i).coeff(0) * b(i-1+dc,1);
    a11 += basis.getBasisPoly(i).coeff(0) * b(i-1+dc,0);
  }
  A(0,0) = T(1.0) - a00;
  A(0,1) = a01;
  A(1,0) = -a10;
  A(1,1) = T(1.0) + a11;
  b(0,0) = std::sin(t);
  b(1,0) = std::cos(t);
  info = Hermite<T>::workspace.solve(2, 1);
#ifdef SACADO_DEBUG
  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for (2x2) solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in (2x2) LU factorization is exactly zero");
#endif
  cs[0] = b(0,0);
  cc[0] = b(1,0);

  // Compute remaining coefficients
  b(0,0) = b0;
  b(1,0) = b1;
  for (unsigned int i=1; i<=dc; i++) {
    cs[i] = cc[0]*b(i-1,0) - cs[0]*b(i-1,1);
    cc[i] = cc[0]*b(i-1+dc,0) - cs[0]*b(i-1+dc,1);
  }

  cs[0] /= basis.getBasisPoly(0).coeff(0);
  cc[0] /= basis.getBasisPoly(0).coeff(0);
}

template <typename T>
Hermite<T>
sin(const Hermite<T>& a)
{
  unsigned int dc = a.degree();
  Hermite<T> s(dc);
  Hermite<T> c(dc);
  sincos(a, s, c);

  return s;
}

template <typename T>
Hermite<T>
cos(const Hermite<T>& a)
{
  unsigned int dc = a.degree();
  Hermite<T> s(dc);
  Hermite<T> c(dc);
  sincos(a, s, c);

  return c;
}

template <typename T>
Hermite<T>
tan(const Hermite<T>& a)
{
  unsigned int dc = a.degree();
  Hermite<T> s(dc);
  Hermite<T> c(dc);
  
  sincos(a, s, c);

  return s / c;
}

template <typename T>
void
sinhcosh(const Hermite<T>& a,
	 Hermite<T>& s,
	 Hermite<T>& c)
{
    unsigned int dc = a.degree();

#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::sinhcosh()";
  TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
		     func << ":  Workspace size (" 
		     << Hermite<T>::workspace.size() 
		     << ") is too small for computation (" << dc+1 
		     << " needed).");
#endif

  if (s.degree() != dc)
    s.resize(dc, false);
  if (c.degree() != dc)
    c.resize(dc, false);

  const T* ca = a.coeff();
  T* cs = s.coeff();
  T* cc = c.coeff();

  typename Hermite<T>::ws_type::matrix_type& A = 
    Hermite<T>::workspace.getMatrix();
  typename Hermite<T>::ws_type::matrix_type& b = 
    Hermite<T>::workspace.getRHS();
  const typename Hermite<T>::ws_type::tp_type& Cijk = 
    Hermite<T>::workspace.getTripleProduct();
  const typename Hermite<T>::ws_type::tp_type::basis_type& basis = 
    Cijk.getBasis();

  // Fill A and b
  A.putScalar(T(0.0));
  b.putScalar(T(0.0));
  for (unsigned int i=1; i<=dc; i++) {
    for (unsigned int j=1; j<=dc; j++) {
      T tmp = T(0.0);
      for (unsigned int k=1; k<=dc; k++) {
	if (i+j >= k) {
	  tmp += Cijk.value(i-1,j,k-1)*basis.derivCoeff(k)*ca[k];
	}
      }
      tmp /= Cijk.norm_squared(i-1)*basis.derivCoeff(i);
      A(i-1,j-1+dc) = -tmp;
      A(i-1+dc,j-1) = -tmp;
    }
    A(i-1,i-1) = T(1.0);
    A(i-1+dc,i-1+dc) = T(1.0);
    b(i-1,0) = ca[i];
    b(i-1+dc,1) = ca[i];
  }

  // Solve system
  int info = Hermite<T>::workspace.solve(2*dc, 2);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

  // Compute degree-0 coefficients
  T t = basis.getBasisPoly(0).coeff(0) * ca[0];
  T a00 = T(1.0);
  T a01 = T(0.0);
  T a10 = T(0.0);
  T a11 = T(1.0);
  T b0 = b(0,0);
  T b1 = b(1,0);
  for (unsigned int i=1; i<=dc; i++) {
    t += basis.getBasisPoly(i).coeff(0) * ca[i];
    a00 += basis.getBasisPoly(i).coeff(0) * b(i-1,1);
    a01 += basis.getBasisPoly(i).coeff(0) * b(i-1,0);
    a10 += basis.getBasisPoly(i).coeff(0) * b(i-1+dc,1);
    a11 += basis.getBasisPoly(i).coeff(0) * b(i-1+dc,0);
  }
  A(0,0) = a00;
  A(0,1) = a01;
  A(1,0) = a10;
  A(1,1) = a11;
  b(0,0) = std::sinh(t);
  b(1,0) = std::cosh(t);
  info = Hermite<T>::workspace.solve(2, 1);
#ifdef SACADO_DEBUG
  TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		     func << ":  Argument " << info 
		     << " for (2x2) solve had illegal value");
  TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		     func << ":  Diagonal entry " << info 
		     << " in (2x2) LU factorization is exactly zero");
#endif
  cs[0] = b(0,0);
  cc[0] = b(1,0);

  // Compute remaining coefficients
  b(0,0) = b0;
  b(1,0) = b1;
  for (unsigned int i=1; i<=dc; i++) {
    cs[i] = cc[0]*b(i-1,0) + cs[0]*b(i-1,1);
    cc[i] = cc[0]*b(i-1+dc,0) + cs[0]*b(i-1+dc,1);
  }

  cs[0] /= basis.getBasisPoly(0).coeff(0);
  cc[0] /= basis.getBasisPoly(0).coeff(0);
}

template <typename T>
Hermite<T>
sinh(const Hermite<T>& a)
{
  unsigned int dc = a.degree();
  Hermite<T> s(dc);
  Hermite<T> c(dc);
  sinhcosh(a, s, c);

  return s;
}

template <typename T>
Hermite<T>
cosh(const Hermite<T>& a)
{
  unsigned int dc = a.degree();
  Hermite<T> s(dc);
  Hermite<T> c(dc);
  sinhcosh(a, s, c);

  return c;
}

template <typename T>
Hermite<T>
tanh(const Hermite<T>& a)
{
  unsigned int dc = a.degree();
  Hermite<T> s(dc);
  Hermite<T> c(dc);
  
  sinhcosh(a, s, c);

  return s / c;
}

template <typename T, typename OpT>
Hermite<T>
quad(const OpT& quad_func,
     const Hermite<T>& a,
     const Hermite<T>& b)
{
  unsigned int dc = a.degree();

  Hermite<T> c(dc);
  const T* ca = a.coeff();
  const T* cb = b.coeff();
  T* cc = c.coeff();

#ifdef SACADO_DEBUG
  const char* func = "Sacado::PCE::Hermite::quad()";
  TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc, std::logic_error,
		     func << ":  Workspace size (" 
		     << Hermite<T>::workspace.size() 
		     << ") is too small for computation (" << dc+1 
		     << " needed).");
#endif

  typename Hermite<T>::ws_type::matrix_type& A = 
    Hermite<T>::workspace.getMatrix();
  typename Hermite<T>::ws_type::matrix_type& B = 
    Hermite<T>::workspace.getRHS();
  const typename Hermite<T>::ws_type::tp_type& Cijk = 
    Hermite<T>::workspace.getTripleProduct();
  const typename Hermite<T>::ws_type::tp_type::basis_type& basis = 
    Cijk.getBasis();
    
    // Fill A
    for (unsigned int i=1; i<=dc; i++) {
      for (unsigned int j=1; j<=dc; j++) {
	A(i-1,j-1) = 0.0;
	for (unsigned int k=0; k<=dc; k++) {
	  A(i-1,j-1) += Cijk.value(i-1,j-1,k)*cb[k];
	}
	A(i-1,j-1) = A(i-1,j-1)*basis.derivCoeff(j) / 
	  (Cijk.norm_squared(i-1)*basis.derivCoeff(i));
      }
    }
    
    // Fill b
    for (unsigned int i=1; i<=dc; i++)
      B(i-1,0) = ca[i];

    // Solve system
    int info = Hermite<T>::workspace.solve(dc, 1);

#ifdef SACADO_DEBUG
    TEST_FOR_EXCEPTION(info < 0, std::logic_error,
		       func << ":  Argument " << info 
		            << " for solve had illegal value");
    TEST_FOR_EXCEPTION(info > 0, std::logic_error,
		       func << ":  Diagonal entry " << info 
		            << " in LU factorization is exactly zero");
#endif

  // Compute degree-0 coefficient
  T s = basis.getBasisPoly(0).coeff(0) * ca[0];
  T t = T(0.0);
  for (unsigned int i=1; i<=dc; i++) {
    s += basis.getBasisPoly(i).coeff(0) * ca[i];
    t += basis.getBasisPoly(i).coeff(0) * B(i-1,0);
  }
  cc[0] = (quad_func(s) - t) / basis.getBasisPoly(0).coeff(0);

  // Get coefficients
  for (unsigned int i=1; i<=dc; i++)
    cc[i] = B(i-1,0);

  return c;
}

template <typename T>
struct acos_quad_func { 
  T operator() (const T& a) const { return std::acos(a); } 
};

template <typename T>
struct asin_quad_func { 
  T operator() (const T& a) const { return std::asin(a); } 
};

template <typename T>
struct atan_quad_func { 
  T operator() (const T& a) const { return std::atan(a); } 
};

template <typename T>
struct acosh_quad_func { 
  T operator() (const T& a) const { return std::log(a+std::sqrt(a*a-T(1.0))); }
};

template <typename T>
struct asinh_quad_func { 
  T operator() (const T& a) const { return std::log(a+std::sqrt(a*a+T(1.0))); }
};

template <typename T>
struct atanh_quad_func { 
  T operator() (const T& a) const { return 0.5*std::log((T(1.0)+a)/(T(1.0)-a)); } 
};

template <typename T>
Hermite<T>
acos(const Hermite<T>& a)
{
  Hermite<T> b = -sqrt(T(1.0) - a*a);
  return quad(acos_quad_func<T>(), a, b);
}

template <typename T>
Hermite<T>
asin(const Hermite<T>& a)
{
  Hermite<T> b = sqrt(T(1.0) - a*a);
  return quad(asin_quad_func<T>(), a, b);
}

template <typename T>
Hermite<T>
atan(const Hermite<T>& a)
{
  Hermite<T> b = T(1.0) + a*a;
  return quad(atan_quad_func<T>(), a, b);
}

// template <typename T>
// Hermite<T>
// atan2(const Hermite<T>& a,
//       const Hermite<T>& b)
// {
//   Hermite<T> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b.coeff(0));
// }

// template <typename T>
// Hermite<T>
// atan2(const T& a,
//       const Hermite<T>& b)
// {
//   Hermite<T> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a,b.coeff(0));
// }

// template <typename T>
// Hermite<T>
// atan2(const Hermite<T>& a,
//       const T& b)
// {
//   Hermite<T> c = atan(a/b);
//   c.fastAccessCoeff(0) = atan2(a.coeff(0),b);
// }

template <typename T>
Hermite<T>
acosh(const Hermite<T>& a)
{
  Hermite<T> b = sqrt(a*a - T(1.0));
  return quad(acosh_quad_func<T>(), a, b);
//   return std::log(a + std::sqrt(a*a - T(1.0)));
}

template <typename T>
Hermite<T>
asinh(const Hermite<T>& a)
{
  Hermite<T> b = sqrt(a*a + T(1.0));
  return quad(asinh_quad_func<T>(), a, b);
//   return std::log(a + std::sqrt(a*a + T(1.0)));
}

template <typename T>
Hermite<T>
atanh(const Hermite<T>& a)
{
  Hermite<T> b = 1.0 - a*a;
  return quad(atanh_quad_func<T>(), a, b);
//   return T(0.5) * std::log((T(1.0) + a)/(T(1.0) - a));
}

template <typename T>
Hermite<T>
fabs(const Hermite<T>& a)
{
  if (a.coeff(0) >= 0)
    return a;
  else
    return -a;
}

template <typename T>
Hermite<T>
abs(const Hermite<T>& a)
{
  if (a.coeff(0) >= 0)
    return a;
  else
    return -a;
}

// template <typename T>
// Hermite<T>
// deriv(const Hermite<T>& a)
// {
//   unsigned int dc = a.degree();

// #ifdef SACADO_DEBUG
//   const char* func = "Sacado::PCE::Hermite::deriv()";
//   TEST_FOR_EXCEPTION(Hermite<T>::workspace.size() < dc+1, std::logic_error,
// 		     func << ":  Workspace size (" 
// 		     << Hermite<T>::workspace.size() 
// 		     << ") is too small for computation (" << dc+1 
// 		     << " needed).");
// #endif

//   Hermite<T> c(dc);
//   const T* ca = a.coeff();
//   T* cc = c.coeff();

//   const typename Hermite<T>::ws_type::tp_type& Cijk = 
//     Hermite<T>::workspace.getTripleProduct();

//   for (unsigned int i=0; i<dc; i++)
//     cc[i] = T(2.0)*(i+1)*ca[i+1] / Cijk.norm_squared(i);
//   cc[dc] = T(0.0);

//   return c;
// }

template <typename T>
Hermite<T>
max(const Hermite<T>& a,
    const Hermite<T>& b)
{
  if (a.coeff(0) >= b.coeff(0))
    return a;
  else
    return b;
}

template <typename T>
Hermite<T>
max(const T& a,
    const Hermite<T>& b)
{
  if (a >= b.coeff(0))
    return Hermite<T>(b.degree(), a);
  else
    return b;
}

template <typename T>
Hermite<T>
max(const Hermite<T>& a,
    const T& b)
{
  if (a.coeff(0) >= b)
    return a;
  else
    return Hermite<T>(a.degree(), b);
}

template <typename T>
Hermite<T>
min(const Hermite<T>& a,
    const Hermite<T>& b)
{
  if (a.coeff(0) <= b.coeff(0))
    return a;
  else
    return b;
}

template <typename T>
Hermite<T>
min(const T& a,
    const Hermite<T>& b)
{
  if (a <= b.coeff(0))
    return Hermite<T>(b.degree(), a);
  else
    return b;
}

template <typename T>
Hermite<T>
min(const Hermite<T>& a,
    const T& b)
{
  if (a.coeff(0) <= b)
    return a;
  else
    return Hermite<T>(a.degree(), b);
}

template <typename T>
bool
operator==(const Hermite<T>& a, 
	   const Hermite<T>& b)
{
  return a.coeff(0) == b.coeff(0);
}

template <typename T>
bool
operator==(const T& a, 
	   const Hermite<T>& b)
{
  return a == b.coeff(0);
}

template <typename T>
bool
operator==(const Hermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) == b;
}

template <typename T>
bool
operator!=(const Hermite<T>& a, 
	   const Hermite<T>& b)
{
  return a.coeff(0) != b.coeff(0);
}

template <typename T>
bool
operator!=(const T& a, 
	   const Hermite<T>& b)
{
  return a != b.coeff(0);
}

template <typename T>
bool
operator!=(const Hermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) != b;
}

template <typename T>
bool
operator<=(const Hermite<T>& a, 
	   const Hermite<T>& b)
{
  return a.coeff(0) <= b.coeff(0);
}

template <typename T>
bool
operator<=(const T& a, 
	   const Hermite<T>& b)
{
  return a <= b.coeff(0);
}

template <typename T>
bool
operator<=(const Hermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) <= b;
}

template <typename T>
bool
operator>=(const Hermite<T>& a, 
	   const Hermite<T>& b)
{
  return a.coeff(0) >= b.coeff(0);
}

template <typename T>
bool
operator>=(const T& a, 
	   const Hermite<T>& b)
{
  return a >= b.coeff(0);
}

template <typename T>
bool
operator>=(const Hermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) >= b;
}

template <typename T>
bool
operator<(const Hermite<T>& a, 
	  const Hermite<T>& b)
{
  return a.coeff(0) < b.coeff(0);
}

template <typename T>
bool
operator<(const T& a, 
	  const Hermite<T>& b)
{
  return a < b.coeff(0);
}

template <typename T>
bool
operator<(const Hermite<T>& a, 
	  const T& b)
{
  return a.coeff(0) < b;
}

template <typename T>
bool
operator>(const Hermite<T>& a, 
	  const Hermite<T>& b)
{
  return a.coeff(0) > b.coeff(0);
}

template <typename T>
bool
operator>(const T& a, 
	  const Hermite<T>& b)
{
  return a > b.coeff(0);
}

template <typename T>
bool
operator>(const Hermite<T>& a, 
	  const T& b)
{
  return a.coeff(0) > b;
}

template <typename T>
std::ostream& 
operator << (std::ostream& os, const Hermite<T>& a)
{
  os << "[ ";
      
  for (unsigned int i=0; i<=a.degree(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}

} // namespace PCE
} // namespace Sacado
