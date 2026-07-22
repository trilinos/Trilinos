// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_SerialDenseMatrix.hpp"


namespace Sacado {
namespace UQ {

/*
template <typename Storage>
KOKKOS_INLINE_FUNCTION
typename PCE<Storage>::value_type
PCE<Storage>::
evaluate(const Teuchos::Array<typename PCE<Storage>::value_type>& point) const
{
  return s_.evaluate(point);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
typename PCE<Storage>::value_type
PCE<Storage>::
evaluate(
  const Teuchos::Array<typename PCE<Storage>::value_type>& point,
  const Teuchos::Array<typename PCE<Storage>::value_type>& bvals) const
{
  return s_.evaluate(point, bvals);
}
*/

template <typename Storage>
KOKKOS_INLINE_FUNCTION
typename PCE<Storage>::value_type
PCE<Storage>::
standard_deviation() const {
  value_type s = 0.0;
  const ordinal_type sz = this->size();
  const_pointer c = this->coeff();
  for (ordinal_type i=1; i<sz; ++i)
    s += c[i]*c[i];
  return std::sqrt(s);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
typename PCE<Storage>::value_type
PCE<Storage>::
two_norm_squared() const {
  value_type s = 0.0;
  const ordinal_type sz = this->size();
  const_pointer c = this->coeff();
  for (ordinal_type i=0; i<sz; ++i)
    s += c[i]*c[i];
  return s;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
typename PCE<Storage>::value_type
PCE<Storage>::
inner_product(const PCE& x) const {
  value_type s = 0.0;
  const ordinal_type sz = this->size();
  const ordinal_type xsz = x.size();
  const ordinal_type n = sz < xsz ? sz : xsz;
  const_pointer c = this->coeff();
  const_pointer xc = x.coeff();
  for (ordinal_type i=0; i<n; ++i)
    s += c[i]*xc[i];
  return s;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
PCE<Storage>::
isEqualTo(const PCE& x) const {
  typedef IsEqual<value_type> IE;
  const ordinal_type sz = this->size();
  if (x.size() != sz) return false;
  bool eq = true;
  for (ordinal_type i=0; i<sz; i++)
    eq = eq && IE::eval(x.coeff(i), this->coeff(i));
  return eq;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
PCE<Storage>::
isEqualTo(const PCE& x) const volatile {
  typedef IsEqual<value_type> IE;
  const ordinal_type sz = this->size();
  if (x.size() != sz) return false;
  bool eq = true;
  for (ordinal_type i=0; i<sz; i++)
    eq = eq && IE::eval(x.coeff(i), this->coeff(i));
  return eq;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator=(const typename PCE<Storage>::value_type v)
{
  s_.init(value_type(0));
  s_[0] = v;
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
/*volatile*/ PCE<Storage>&
PCE<Storage>::
operator=(const typename PCE<Storage>::value_type v) volatile
{
  s_.init(value_type(0));
  s_[0] = v;
  return const_cast<PCE&>(*this);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator=(const PCE<Storage>& x)
{
  if (this != &x) {
    if (!s_.is_view())
      cijk_ = x.cijk_;
    if (!s_.is_view() && is_constant(x)) {
      s_.resize(1);
      s_[0] = x.s_[0];
    }
    else
      s_ = x.s_;

    // For DyamicStorage as a view (is_owned=false), we need to set
    // the trailing entries when assigning a constant vector (because
    // the copy constructor in this case doesn't reset the size of this)
    if (s_.size() > x.s_.size())
      for (ordinal_type i=x.s_.size(); i<s_.size(); i++)
        s_[i] = value_type(0);
  }
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator=(const volatile PCE<Storage>& x)
{
  if (this != &x) {
    if (!s_.is_view())
      cijk_ = const_cast<const my_cijk_type&>(x.cijk_);
    if (!s_.is_view() && is_constant(x)) {
      s_.resize(1);
      s_[0] = x.s_[0];
    }
    else
      s_ = x.s_;

    // For DyamicStorage as a view (is_owned=false), we need to set
    // the trailing entries when assigning a constant vector (because
    // the copy constructor in this case doesn't reset the size of this)
    if (s_.size() > x.s_.size())
      for (ordinal_type i=x.s_.size(); i<s_.size(); i++)
        s_[i] = value_type(0);
  }
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
/*volatile*/ PCE<Storage>&
PCE<Storage>::
operator=(const PCE<Storage>& x) volatile
{
  if (this != &x) {
    if (!s_.is_view())
      const_cast<my_cijk_type&>(cijk_) = x.cijk_;
    if (!s_.is_view() && is_constant(x)) {
      s_.resize(1);
      s_[0] = x.s_[0];
    }
    else
      s_ = x.s_;

    // For DyamicStorage as a view (is_owned=false), we need to set
    // the trailing entries when assigning a constant vector (because
    // the copy constructor in this case doesn't reset the size of this)
    if (s_.size() > x.s_.size())
      for (ordinal_type i=x.s_.size(); i<s_.size(); i++)
        s_[i] = value_type(0);
  }
  return const_cast<PCE&>(*this);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
/*volatile*/ PCE<Storage>&
PCE<Storage>::
operator=(const volatile PCE<Storage>& x) volatile
{
  if (this != &x) {
    if (!s_.is_view())
      const_cast<my_cijk_type&>(cijk_) =
        const_cast<const my_cijk_type&>(x.cijk_);
    if (!s_.is_view() && is_constant(x)) {
      s_.resize(1);
      s_[0] = x.s_[0];
    }
    else
      s_ = x.s_;

    // For DyamicStorage as a view (is_owned=false), we need to set
    // the trailing entries when assigning a constant vector (because
    // the copy constructor in this case doesn't reset the size of this)
    if (s_.size() > x.s_.size())
      for (ordinal_type i=x.s_.size(); i<s_.size(); i++)
        s_[i] = value_type(0);
  }
  return const_cast<PCE&>(*this);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
PCE<Storage>::
operator-() const
{
  const ordinal_type sz = this->size();
  PCE<Storage> x(cijk_, sz);
  pointer xc = x.coeff();
  const_pointer cc = this->coeff();
  for (ordinal_type i=0; i<sz; ++i)
    xc[i] = -cc[i];
  return x;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
PCE<Storage>::
operator-() const volatile
{
  const ordinal_type sz = this->size();
  PCE<Storage> x(const_cast<const my_cijk_type&>(cijk_), sz);
  pointer xc = x.coeff();
  const_volatile_pointer cc = this->coeff();
  for (ordinal_type i=0; i<sz; ++i)
    xc[i] = -cc[i];
  return x;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator*=(const typename PCE<Storage>::value_type v)
{
  pointer cc = this->coeff();
  const ordinal_type sz = this->size();
  for (ordinal_type i=0; i<sz; ++i)
    cc[i] *= v;
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
/*volatile*/ PCE<Storage>&
PCE<Storage>::
operator*=(const typename PCE<Storage>::value_type v) volatile
{
  volatile_pointer cc = this->coeff();
  const ordinal_type sz = this->size();
  for (ordinal_type i=0; i<sz; ++i)
    cc[i] *= v;
  return const_cast<PCE&>(*this);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator/=(const typename PCE<Storage>::value_type v)
{
  pointer cc = this->coeff();
  const ordinal_type sz = this->size();
  for (ordinal_type i=0; i<sz; ++i)
    cc[i] /= v;
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
/*volatile*/ PCE<Storage>&
PCE<Storage>::
operator/=(const typename PCE<Storage>::value_type v) volatile
{
  volatile pointer cc = this->coeff();
  const ordinal_type sz = this->size();
  for (ordinal_type i=0; i<sz; ++i)
    cc[i] /= v;
  return const_cast<PCE&>(*this);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator+=(const PCE<Storage>& x)
{
  const ordinal_type xsz = x.size();
  const ordinal_type sz = this->size();
  if (xsz > sz) {
    this->reset(x.cijk_, xsz);
  }
  const_pointer xc = x.coeff();
  pointer cc = this->coeff();
  for (ordinal_type i=0; i<xsz; ++i)
    cc[i] += xc[i];
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator-=(const PCE<Storage>& x)
{
  const ordinal_type xsz = x.size();
  const ordinal_type sz = this->size();
  if (xsz > sz) {
    this->reset(x.cijk_, xsz);
  }
  const_pointer xc = x.coeff();
  pointer cc = this->coeff();
  for (ordinal_type i=0; i<xsz; ++i)
    cc[i] -= xc[i];
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator*=(const PCE<Storage>& x)
{
  const ordinal_type xsz = x.size();
  const ordinal_type sz = this->size();
  const ordinal_type csz = sz > xsz ? sz : xsz;


  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    sz != xsz && sz != 1 && xsz != 1, std::logic_error,
    "Sacado::UQ::PCE::operator*=(): input sizes do not match");
  )

  if (cijk_.is_empty() && !x.cijk_.is_empty())
    cijk_ = x.cijk_;

  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    cijk_.is_empty() && csz != 1, std::logic_error,
    "Sacado::UQ::PCE::operator*(): empty cijk but expansion size > 1");
  )

  if (csz > sz)
    s_.resize(csz);

  const_pointer xc = x.coeff();
  pointer cc = this->coeff();
  if (xsz == 1) {
    const value_type xcz = xc[0];
    for (ordinal_type i=0; i<sz; ++i)
      cc[i] *= xcz;
  }
  else if (sz == 1) {
    const value_type ccz = cc[0];
    for (ordinal_type i=0; i<xsz; ++i)
      cc[i] = ccz*xc[i];
  }
  else {
    PCE<Storage> y(cijk_, csz);
    pointer yc = y.coeff();
    for (ordinal_type i=0; i<csz; ++i) {
      const cijk_size_type num_entry = cijk_.num_entry(i);
      const cijk_size_type entry_beg = cijk_.entry_begin(i);
      const cijk_size_type entry_end = entry_beg + num_entry;
      value_type ytmp = 0;
      for (cijk_size_type entry = entry_beg; entry < entry_end; ++entry) {
        const cijk_size_type j = cijk_.coord(entry,0);
        const cijk_size_type k = cijk_.coord(entry,1);
        ytmp += cijk_.value(entry) * ( cc[j] * xc[k] + cc[k] * xc[j] );
      }
      yc[i] = ytmp;
    }
    s_ = y.s_;
  }
  return *this;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>&
PCE<Storage>::
operator/=(const PCE<Storage>& x)
{
  const ordinal_type xsz = x.size();
  const ordinal_type sz = this->size();
  const ordinal_type csz = sz > xsz ? sz : xsz;

  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    sz != xsz && sz != 1 && xsz != 1, std::logic_error,
    "Sacado::UQ::PCE::operator/=(): input sizes do not match");
  )

  if (cijk_.is_empty() && !x.cijk_.is_empty())
    cijk_ = x.cijk_;

  if (csz > sz)
    s_.resize(csz);

  const_pointer xc = x.coeff();
  pointer cc = this->coeff();

  KOKKOS_IF_ON_DEVICE(
  const value_type xcz = xc[0];
    for (ordinal_type i=0; i<sz; ++i)
      cc[i] /= xcz;
  )

  KOKKOS_IF_ON_HOST(
  if (xsz == 1) {//constant denom
    const value_type xcz = xc[0];
    for (ordinal_type i=0; i<sz; ++i)
      cc[i] /= xcz;
  }
  else {

    PCE<Storage> y(cijk_, csz);
    CG_divide(*this, x, y);
    s_ = y.s_;
  }
  )

  return *this;

}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator+(const PCE<Storage>& a, const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::my_cijk_type my_cijk_type;
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type asz = a.size();
  const ordinal_type bsz = b.size();
  const ordinal_type csz = asz > bsz ? asz : bsz;
  my_cijk_type c_cijk = asz > bsz ? a.cijk() : b.cijk();

  PCE<Storage> c(c_cijk, csz);
  const_pointer ac = a.coeff();
  const_pointer bc = b.coeff();
  pointer cc = c.coeff();
  if (asz > bsz) {
    for (ordinal_type i=0; i<bsz; ++i)
      cc[i] = ac[i] + bc[i];
    for (ordinal_type i=bsz; i<asz; ++i)
      cc[i] = ac[i];
  }
  else {
    for (ordinal_type i=0; i<asz; ++i)
      cc[i] = ac[i] + bc[i];
    for (ordinal_type i=asz; i<bsz; ++i)
      cc[i] = bc[i];
  }

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator+(const typename PCE<Storage>::value_type& a,
          const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type bsz = b.size();
  PCE<Storage> c(b.cijk(), bsz);
  const_pointer bc = b.coeff();
  pointer cc = c.coeff();
  cc[0] = a + bc[0];
  for (ordinal_type i=1; i<bsz; ++i)
    cc[i] = bc[i];
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator+(const PCE<Storage>& a,
          const typename PCE<Storage>::value_type& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type asz = a.size();
  PCE<Storage> c(a.cijk(), asz);
  const_pointer ac = a.coeff();
  pointer cc = c.coeff();
  cc[0] = ac[0] + b;
  for (ordinal_type i=1; i<asz; ++i)
    cc[i] = ac[i];
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator-(const PCE<Storage>& a, const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::my_cijk_type my_cijk_type;
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type asz = a.size();
  const ordinal_type bsz = b.size();
  const ordinal_type csz = asz > bsz ? asz : bsz;
  my_cijk_type c_cijk = asz > bsz ? a.cijk() : b.cijk();

  PCE<Storage> c(c_cijk, csz);
  const_pointer ac = a.coeff();
  const_pointer bc = b.coeff();
  pointer cc = c.coeff();
  if (asz > bsz) {
    for (ordinal_type i=0; i<bsz; ++i)
      cc[i] = ac[i] - bc[i];
    for (ordinal_type i=bsz; i<asz; ++i)
      cc[i] = ac[i];
  }
  else {
    for (ordinal_type i=0; i<asz; ++i)
      cc[i] = ac[i] - bc[i];
    for (ordinal_type i=asz; i<bsz; ++i)
      cc[i] = -bc[i];
  }

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator-(const typename PCE<Storage>::value_type& a,
          const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type bsz = b.size();
  PCE<Storage> c(b.cijk(), bsz);
  const_pointer bc = b.coeff();
  pointer cc = c.coeff();
  cc[0] = a - bc[0];
  for (ordinal_type i=1; i<bsz; ++i)
    cc[i] = -bc[i];
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator-(const PCE<Storage>& a,
          const typename PCE<Storage>::value_type& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type asz = a.size();
  PCE<Storage> c(a.cijk(), asz);
  const_pointer ac = a.coeff();
  pointer cc = c.coeff();
  cc[0] = ac[0] - b;
  for (ordinal_type i=1; i<asz; ++i)
    cc[i] = ac[i];
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator*(const PCE<Storage>& a, const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::my_cijk_type my_cijk_type;
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;
  typedef typename PCE<Storage>::value_type value_type;
  typedef typename my_cijk_type::size_type cijk_size_type;

  const ordinal_type asz = a.size();
  const ordinal_type bsz = b.size();
  const ordinal_type csz = asz > bsz ? asz : bsz;

  KOKKOS_IF_ON_HOST((
  TEUCHOS_TEST_FOR_EXCEPTION(
    asz != bsz && asz != 1 && bsz != 1, std::logic_error,
    "Sacado::UQ::PCE::operator*(): input sizes do not match");
  ))

  my_cijk_type c_cijk = a.cijk().is_empty() ? b.cijk() : a.cijk();

  KOKKOS_IF_ON_HOST((
  TEUCHOS_TEST_FOR_EXCEPTION(
    c_cijk.is_empty() && csz != 1, std::logic_error,
    "Sacado::UQ::PCE::operator*(): empty cijk but expansion size > 1");
  ))

  PCE<Storage> c(c_cijk, csz);
  const_pointer ac = a.coeff();
  const_pointer bc = b.coeff();
  pointer cc = c.coeff();

  if (asz == 1) {
    const value_type acz = ac[0];
    for (ordinal_type i=0; i<csz; ++i)
      cc[i] = acz * bc[i];
  }
  else if (bsz == 1) {
    const value_type bcz = bc[0];
    for (ordinal_type i=0; i<csz; ++i)
      cc[i] = ac[i] * bcz;
  }
  else {
    for (ordinal_type i=0; i<csz; ++i) {
      const cijk_size_type num_entry = c_cijk.num_entry(i);
      const cijk_size_type entry_beg = c_cijk.entry_begin(i);
      const cijk_size_type entry_end = entry_beg + num_entry;
      value_type ytmp = 0;
      for (cijk_size_type entry = entry_beg; entry < entry_end; ++entry) {
        const cijk_size_type j = c_cijk.coord(entry,0);
        const cijk_size_type k = c_cijk.coord(entry,1);
        ytmp += c_cijk.value(entry) * ( ac[j] * bc[k] + ac[k] * bc[j] );
      }
      cc[i] = ytmp;
    }
  }

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator*(const typename PCE<Storage>::value_type& a,
          const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type bsz = b.size();
  PCE<Storage> c(b.cijk(), bsz);
  const_pointer bc = b.coeff();
  pointer cc = c.coeff();
  for (ordinal_type i=0; i<bsz; ++i)
    cc[i] = a*bc[i];
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator*(const PCE<Storage>& a,
          const typename PCE<Storage>::value_type& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type asz = a.size();
  PCE<Storage> c(a.cijk(), asz);
  const_pointer ac = a.coeff();
  pointer cc = c.coeff();
  for (ordinal_type i=0; i<asz; ++i)
    cc[i] = ac[i]*b;
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator/(const PCE<Storage>& a, const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;
  typedef typename PCE<Storage>::value_type value_type;
  typedef typename PCE<Storage>::my_cijk_type my_cijk_type;

  const ordinal_type asz = a.size();
  const ordinal_type bsz = b.size();
  const ordinal_type csz = asz > bsz ? asz : bsz;
  
  KOKKOS_IF_ON_HOST((
  TEUCHOS_TEST_FOR_EXCEPTION(
  asz != bsz && asz != 1 && bsz != 1, std::logic_error,
  "Sacado::UQ::PCE::operator/(): input sizes do not match");
  ))
  my_cijk_type c_cijk = asz == bsz || asz >1 ? a.cijk() : b.cijk();

  PCE<Storage> c(c_cijk, csz);

  KOKKOS_IF_ON_HOST((
  const_pointer ac = a.coeff();
  pointer cc = c.coeff();
  value_type bcz = b.fastAccessCoeff(0);
  for (ordinal_type i=0; i<asz; ++i)
    cc[i] = ac[i]/bcz;
  ))

  KOKKOS_IF_ON_HOST((
  if (bsz == 1) {//constant denom
    const_pointer ac = a.coeff();
    const_pointer bc = b.coeff();
    pointer cc = c.coeff();
    const value_type bcz = bc[0];
    for (ordinal_type i=0; i<csz; ++i)
      cc[i] = ac[i] / bcz;
  }
  else {
   CG_divide(a,b,c);
  }
  ))

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator/(const typename PCE<Storage>::value_type& a,
          const PCE<Storage>& b)
{
  //Creat a 0-th order PCE for a
  PCE<Storage> a_pce(a);
  return operator/(a_pce,b);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
operator/(const PCE<Storage>& a,
          const typename PCE<Storage>::value_type& b)
{
  typedef typename PCE<Storage>::pointer pointer;
  typedef typename PCE<Storage>::const_pointer const_pointer;
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  const ordinal_type asz = a.size();
  PCE<Storage> c(a.cijk(), asz);
  const_pointer ac = a.coeff();
  pointer cc = c.coeff();
  for (ordinal_type i=0; i<asz; ++i)
    cc[i] = ac[i]/b;
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
exp(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::exp():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::exp( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
log(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::log():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::log( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
log10(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::log10():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::log10( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
sqrt(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::sqrt():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::sqrt( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
cbrt(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::cbrt():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::cbrt( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
pow(const PCE<Storage>& a, const PCE<Storage>& b)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1 || b.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::pow():  arguments have size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::pow(a.fastAccessCoeff(0), b.fastAccessCoeff(0));

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
pow(const typename PCE<Storage>::value_type& a,
    const PCE<Storage>& b)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    b.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::pow():  arguments have size != 1");
  )

  PCE<Storage> c(b.cijk(), 1);
  c.fastAccessCoeff(0) = std::pow(a, b.fastAccessCoeff(0));

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
pow(const PCE<Storage>& a,
    const typename PCE<Storage>::value_type& b)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::pow():  arguments have size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::pow(a.fastAccessCoeff(0), b);

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
sin(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::sin():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::sin( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
cos(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::cos():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::cos( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
tan(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::tan():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::tan( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
sinh(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::sinh():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::sinh( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
cosh(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::cosh():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::cosh( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
tanh(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::tanh():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::tanh( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
acos(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::acos():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::acos( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
asin(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::asin():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::asin( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
atan(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::atan():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::atan( a.fastAccessCoeff(0) );

  return c;
}

/*
template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
acosh(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::acosh():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::acosh( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
asinh(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::asinh():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::asinh( a.fastAccessCoeff(0) );

  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
atanh(const PCE<Storage>& a)
{
  KOKKOS_IF_ON_HOST(
  TEUCHOS_TEST_FOR_EXCEPTION(
    a.size() != 1, std::logic_error,
    "Sacado::UQ::PCE::atanh():  argument has size != 1");
  )

  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::atanh( a.fastAccessCoeff(0) );

  return c;
}
*/

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
fabs(const PCE<Storage>& a)
{
  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = a.two_norm();
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
abs(const PCE<Storage>& a)
{
  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = a.two_norm();
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
ceil(const PCE<Storage>& a)
{
  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = std::ceil( a.fastAccessCoeff(0) );
  return c;
}

// template <typename Storage>
// KOKKOS_INLINE_FUNCTION
// PCE<Storage>
// max(const PCE<Storage>& a, const PCE<Storage>& b)
// {
//   if (a.two_norm() >= b.two_norm())
//     return a;
//   return b;
// }

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
max(const typename PCE<Storage>::value_type& a,
    const PCE<Storage>& b)
{
  if (a >= b.two_norm()) {
    PCE<Storage> c(b.cijk(), 1);
    c.fastAccessCoeff(0) = a;
    return c;
  }
  return b;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
max(const PCE<Storage>& a,
    const typename PCE<Storage>::value_type& b)
{
  if (a.two_norm() >= b)
    return a;
  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = b;
  return c;
}

// template <typename Storage>
// KOKKOS_INLINE_FUNCTION
// PCE<Storage>
// min(const PCE<Storage>& a, const PCE<Storage>& b)
// {
//   if (a.two_norm() <= b.two_norm())
//     return a;
//   return b;
// }

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
min(const typename PCE<Storage>::value_type& a,
    const PCE<Storage>& b)
{
  if (a <= b.two_norm()) {
    PCE<Storage> c(b.cijk(), 1);
    c.fastAccessCoeff(0) = a;
    return c;
  }
  return b;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
PCE<Storage>
min(const PCE<Storage>& a,
    const typename PCE<Storage>::value_type& b)
{
  if (a.two_norm() <= b)
    return a;
  PCE<Storage> c(a.cijk(), 1);
  c.fastAccessCoeff(0) = b;
  return c;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator==(const PCE<Storage>& a, const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::ordinal_type ordinal_type;
  const ordinal_type asz = a.size();
  const ordinal_type bsz = b.size();
  const ordinal_type n = asz > bsz ? asz : bsz;
  for (ordinal_type i=0; i<n; i++)
    if (a.coeff(i) != b.coeff(i))
      return false;
  return true;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator==(const typename PCE<Storage>::value_type& a,
           const PCE<Storage>& b)
{
  typedef typename PCE<Storage>::ordinal_type ordinal_type;
  const ordinal_type n = b.size();
  if (a != b.coeff(0))
    return false;
  for (ordinal_type i=1; i<n; i++)
    if (b.fastAccessCoeff(i) != typename PCE<Storage>::value_type(0.0))
      return false;
  return true;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator==(const PCE<Storage>& a,
           const typename PCE<Storage>::value_type& b)
{
  typedef typename PCE<Storage>::ordinal_type ordinal_type;
  const ordinal_type n = a.size();
  if (a.coeff(0) != b)
    return false;
  for (ordinal_type i=1; i<n; i++)
    if (a.fastAccessCoeff(i) != typename PCE<Storage>::value_type(0.0))
      return false;
  return true;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator!=(const PCE<Storage>& a, const PCE<Storage>& b)
{
  return !(a == b);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator!=(const typename PCE<Storage>::value_type& a,
           const PCE<Storage>& b)
{
  return !(a == b);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator!=(const PCE<Storage>& a,
           const typename PCE<Storage>::value_type& b)
{
  return !(a == b);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator<=(const PCE<Storage>& a, const PCE<Storage>& b)
{
  return a.two_norm() <= b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator<=(const typename PCE<Storage>::value_type& a,
           const PCE<Storage>& b)
{
  return a <= b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator<=(const PCE<Storage>& a,
           const typename PCE<Storage>::value_type& b)
{
  return a.two_norm() <= b;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator>=(const PCE<Storage>& a, const PCE<Storage>& b)
{
  return a.two_norm() >= b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator>=(const typename PCE<Storage>::value_type& a,
           const PCE<Storage>& b)
{
  return a >= b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator>=(const PCE<Storage>& a,
           const typename PCE<Storage>::value_type& b)
{
  return a.two_norm() >= b;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator<(const PCE<Storage>& a, const PCE<Storage>& b)
{
  return a.two_norm() < b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator<(const typename PCE<Storage>::value_type& a,
          const PCE<Storage>& b)
{
  return a < b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator<(const PCE<Storage>& a,
          const typename PCE<Storage>::value_type& b)
{
  return a.two_norm() < b;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator>(const PCE<Storage>& a, const PCE<Storage>& b)
{
  return a.two_norm() > b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator>(const typename PCE<Storage>::value_type& a,
          const PCE<Storage>& b)
{
  return a > b.two_norm();
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator>(const PCE<Storage>& a,
          const typename PCE<Storage>::value_type& b)
{
  return a.two_norm() > b;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool toBool(const PCE<Storage>& x) {
  typedef typename PCE<Storage>::ordinal_type ordinal_type;
  bool is_zero = true;
  const ordinal_type sz = x.size();
  for (ordinal_type i=0; i<sz; i++)
    is_zero = is_zero && (x.fastAccessCoeff(i) == 0.0);
  return !is_zero;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator && (const PCE<Storage>& x1, const PCE<Storage>& x2)
{
  return toBool(x1) && toBool(x2);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator && (const typename PCE<Storage>::value_type& a,
             const PCE<Storage>& x2)
{
  return a && toBool(x2);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator && (const PCE<Storage>& x1,
             const typename PCE<Storage>::value_type& b)
{
  return toBool(x1) && b;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator || (const PCE<Storage>& x1, const PCE<Storage>& x2)
{
  return toBool(x1) || toBool(x2);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator || (const typename PCE<Storage>::value_type& a,
             const PCE<Storage>& x2)
{
  return a || toBool(x2);
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
bool
operator || (const PCE<Storage>& x1,
             const typename PCE<Storage>::value_type& b)
{
  return toBool(x1) || b;
}

template <typename Storage>
std::ostream&
operator << (std::ostream& os, const PCE<Storage>& a)
{
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  os << "[ ";

  for (ordinal_type i=0; i<a.size(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]";
  return os;
}

template <typename Storage>
std::istream&
operator >> (std::istream& is, PCE<Storage>& a)
{
  typedef typename PCE<Storage>::ordinal_type ordinal_type;

  // Read in the opening "["
  char bracket;
  is >> bracket;

  for (ordinal_type i=0; i<a.size(); i++) {
    is >> a.fastAccessCoeff(i);
  }

  // Read in the closing "]"

  is >> bracket;
  return is;
}

template <typename Storage>
KOKKOS_INLINE_FUNCTION
void
CG_divide(const PCE<Storage>& a, const PCE<Storage>& b, PCE<Storage>& c) {
    typedef typename PCE<Storage>::ordinal_type ordinal_type;
    typedef typename PCE<Storage>::value_type value_type;

    const ordinal_type size = c.size();

    //Needed scalars
    value_type alpha, beta, rTz, rTz_old, resid;

    //Needed temporary PCEs 
    PCE<Storage> r(a.cijk(),size);
    PCE<Storage> p(a.cijk(),size);
    PCE<Storage> bp(a.cijk(),size);
    PCE<Storage> z(a.cijk(),size);

    //compute residual = a - b*c 
    r =  a - b*c;
    z = r/b.coeff(0);
    p = z;
    resid = r.two_norm();
    //Compute <r,z>=rTz (L2 inner product)
    rTz = r.inner_product(z);
    ordinal_type k = 0;
    value_type tol = 1e-6;
    while ( resid > tol && k < 100){
      bp = b*p;
      //Compute alpha = <r,z>/<p,b*p>
      alpha = rTz/p.inner_product(bp);
      //Update the solution c = c + alpha*p
      c = c + alpha*p;
      rTz_old = rTz;
      //Compute the new residual r = r - alpha*b*p
      r = r - alpha*bp;
      resid = r.two_norm();
      //Compute beta = rTz_new/ rTz_old and new p
      z = r/b.coeff(0);
      rTz = r.inner_product(z);
      beta = rTz/rTz_old;
      p = z + beta*p;
      k++;
   }
 }

} // namespace UQ
} // namespace Sacado
