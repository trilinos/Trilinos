// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Sacado_DynamicArrayTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_ConstantOrthogPolyExpansion.hpp"

namespace Sacado {
namespace PCE {

template <typename T, typename Storage> 
OrthogPoly<T,Storage>::
OrthogPoly() :
  expansion_(),
  th(new Stokhos::OrthogPolyApprox<int,value_type,Storage>)
{
  const_expansion_ = Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<int,T>);
  expansion_ = const_expansion_;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>::
OrthogPoly(const typename OrthogPoly<T,Storage>::value_type& x) :
  expansion_(),
  th(new Stokhos::OrthogPolyApprox<int,value_type,Storage>(Teuchos::null, 1, &x))
{
  const_expansion_ = Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<int,T>);
  expansion_ = const_expansion_;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>::
OrthogPoly(const Teuchos::RCP<expansion_type>& expansion) :
  expansion_(expansion),
  th(new Stokhos::OrthogPolyApprox<int,value_type,Storage>(expansion_->getBasis()))
{
  const_expansion_ = Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<int,T>);
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>::
OrthogPoly(const Teuchos::RCP<expansion_type>& expansion,
	   ordinal_type sz) :
  expansion_(expansion),
  th(new Stokhos::OrthogPolyApprox<int,value_type,Storage>(expansion_->getBasis(), sz))
{
  const_expansion_ = Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<int,T>);
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>::
OrthogPoly(const OrthogPoly<T,Storage>& x) :
  expansion_(x.expansion_),
  th(x.th)
{
  const_expansion_ = Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<int,T>);
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>::
~OrthogPoly()
{
}

template <typename T, typename Storage> 
void
OrthogPoly<T,Storage>::
reset(const Teuchos::RCP<expansion_type>& expansion)
{
  expansion_ = expansion;
  th->reset(expansion_->getBasis());
}

template <typename T, typename Storage> 
void
OrthogPoly<T,Storage>::
reset(const Teuchos::RCP<expansion_type>& expansion, ordinal_type sz)
{
  expansion_ = expansion;
  th->reset(expansion_->getBasis(), sz);
}

template <typename T, typename Storage> 
typename OrthogPoly<T,Storage>::value_type
OrthogPoly<T,Storage>::
evaluate(const Teuchos::Array<typename OrthogPoly<T,Storage>::value_type>& point) const
{
  return th->evaluate(point);
}

template <typename T, typename Storage> 
typename OrthogPoly<T,Storage>::value_type
OrthogPoly<T,Storage>::
evaluate(
  const Teuchos::Array<typename OrthogPoly<T,Storage>::value_type>& point,
  const Teuchos::Array<typename OrthogPoly<T,Storage>::value_type>& bvals) const
{
  return th->evaluate(point, bvals);
}

template <typename T, typename Storage> 
bool 
OrthogPoly<T,Storage>::
isEqualTo(const OrthogPoly& x) const {
  typedef IsEqual<value_type> IE;
  if (x.size() != this->size()) return false;
  // Allow expansions to be different if their size is 1 and one
  // of them is a constant expansion
  if (expansion_ != x.expansion_) {
    if (x.size() != 1)
      return false;
    if ((expansion_ != const_expansion_) && 
	(x.expansion_ != x.const_expansion_))
      return false;
  }
  bool eq = true;
  for (int i=0; i<this->size(); i++)
    eq = eq && IE::eval(x.coeff(i), this->coeff(i));
  return eq;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator=(const typename OrthogPoly<T,Storage>::value_type& v) 
{
  th.makeOwnCopy();
  *th = v;
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator=(const OrthogPoly<T,Storage>& x) 
{
  expansion_ = x.expansion_;
  th = x.th;
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>
OrthogPoly<T,Storage>::
operator+() const
{
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage> 
OrthogPoly<T,Storage>::
operator-() const
{
  OrthogPoly<T,Storage> x(expansion_);
  expansion_->unaryMinus(*(x.th), *th);
  return x;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator+=(const typename OrthogPoly<T,Storage>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->plusEqual(*th, v);
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator-=(const typename OrthogPoly<T,Storage>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->minusEqual(*th, v);
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator*=(const typename OrthogPoly<T,Storage>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->timesEqual(*th, v);
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator/=(const typename OrthogPoly<T,Storage>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->divideEqual(*th, v);
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator+=(const OrthogPoly<T,Storage>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->plusEqual(*th, *x.th);
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator-=(const OrthogPoly<T,Storage>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->minusEqual(*th, *x.th);
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator*=(const OrthogPoly<T,Storage>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->timesEqual(*th, *x.th);
  return *this;
}

template <typename T, typename Storage> 
OrthogPoly<T,Storage>& 
OrthogPoly<T,Storage>::
operator/=(const OrthogPoly<T,Storage>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->divideEqual(*th, *x.th);
  return *this;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator+(const OrthogPoly<T,Storage>& a, 
	  const OrthogPoly<T,Storage>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T,Storage> c(e, 0);
  e->plus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	  b.getOrthogPolyApprox());

  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator+(const typename OrthogPoly<T,Storage>::value_type& a, 
	  const OrthogPoly<T,Storage>& b)
{
  OrthogPoly<T,Storage> c(b.expansion(), 0);
  b.expansion()->plus(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator+(const OrthogPoly<T,Storage>& a, 
	  const typename OrthogPoly<T,Storage>::value_type& b)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->plus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator-(const OrthogPoly<T,Storage>& a, 
	  const OrthogPoly<T,Storage>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T,Storage> c(e, 0);
  e->minus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	   b.getOrthogPolyApprox());

  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator-(const typename OrthogPoly<T,Storage>::value_type& a, 
	  const OrthogPoly<T,Storage>& b)
{
  OrthogPoly<T,Storage> c(b.expansion(), 0);
  b.expansion()->minus(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator-(const OrthogPoly<T,Storage>& a, 
	  const typename OrthogPoly<T,Storage>::value_type& b)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->minus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator*(const OrthogPoly<T,Storage>& a, 
	  const OrthogPoly<T,Storage>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T,Storage> c(e, 0);
  e->times(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	   b.getOrthogPolyApprox());

  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator*(const typename OrthogPoly<T,Storage>::value_type& a, 
	  const OrthogPoly<T,Storage>& b)
{
  OrthogPoly<T,Storage> c(b.expansion(), 0);
  b.expansion()->times(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator*(const OrthogPoly<T,Storage>& a, 
	  const typename OrthogPoly<T,Storage>::value_type& b)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->times(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator/(const OrthogPoly<T,Storage>& a, 
	  const OrthogPoly<T,Storage>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T,Storage> c(e, 0);
  e->divide(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	    b.getOrthogPolyApprox());

  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator/(const typename OrthogPoly<T,Storage>::value_type& a, 
	  const OrthogPoly<T,Storage>& b)
{
  OrthogPoly<T,Storage> c(b.expansion(), 0);
  b.expansion()->divide(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
operator/(const OrthogPoly<T,Storage>& a, 
	  const typename OrthogPoly<T,Storage>::value_type& b)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->divide(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
exp(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->exp(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
log(const OrthogPoly<T,Storage>& a)
{
  TEUCHOS_FUNC_TIME_MONITOR("LOG");
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  {
    TEUCHOS_FUNC_TIME_MONITOR("OPA LOG");
    a.expansion()->log(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  }

  OrthogPoly<T,Storage> d(c);
  return d;
}

template <typename T, typename Storage>
void
log(OrthogPoly<T,Storage>& c, const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> d(a.expansion(), 0);
  a.expansion()->log(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
log10(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->log10(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
sqrt(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->sqrt(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
cbrt(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->cbrt(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
pow(const OrthogPoly<T,Storage>& a, 
    const OrthogPoly<T,Storage>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T,Storage> c(e, 0);
  e->pow(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	 b.getOrthogPolyApprox());

  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
pow(const T& a,
    const OrthogPoly<T,Storage>& b)
{
  OrthogPoly<T,Storage> c(b.expansion(), 0);
  b.expansion()->pow(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
pow(const OrthogPoly<T,Storage>& a,
    const T& b)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->pow(c.getOrthogPolyApprox(),a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
sin(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->sin(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
cos(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->cos(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
tan(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->tan(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
sinh(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->sinh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
cosh(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->cosh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
tanh(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->tanh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
acos(const OrthogPoly<T,Storage>& a)
{
   OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->acos(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
asin(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->asin(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
atan(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->atan(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
acosh(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->acosh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
asinh(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->asinh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
atanh(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->atanh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
fabs(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->fabs(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
abs(const OrthogPoly<T,Storage>& a)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->abs(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
max(const OrthogPoly<T,Storage>& a,
    const OrthogPoly<T,Storage>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T,Storage> c(e, 0);
  e->max(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	 b.getOrthogPolyApprox());

  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
max(const typename OrthogPoly<T,Storage>::value_type& a,
    const OrthogPoly<T,Storage>& b)
{
  OrthogPoly<T,Storage> c(b.expansion(), 0);
  b.expansion()->max(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
max(const OrthogPoly<T,Storage>& a,
    const typename OrthogPoly<T,Storage>::value_type& b)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->max(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
min(const OrthogPoly<T,Storage>& a,
    const OrthogPoly<T,Storage>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T,Storage>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T,Storage> c(e, 0);
  e->min(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	 b.getOrthogPolyApprox());

  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
min(const typename OrthogPoly<T,Storage>::value_type& a,
    const OrthogPoly<T,Storage>& b)
{
  OrthogPoly<T,Storage> c(b.expansion(), 0);
  b.expansion()->min(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T, typename Storage>
OrthogPoly<T,Storage>
min(const OrthogPoly<T,Storage>& a,
    const typename OrthogPoly<T,Storage>::value_type& b)
{
  OrthogPoly<T,Storage> c(a.expansion(), 0);
  a.expansion()->min(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T, typename Storage>
bool
operator==(const OrthogPoly<T,Storage>& a, 
	   const OrthogPoly<T,Storage>& b)
{
  int n = std::max(a.size(), b.size());
  for (int i=0; i<n; i++)
    if (a.coeff(i) != b.coeff(i))
      return false;
  return true;
}

template <typename T, typename Storage>
bool
operator==(const typename OrthogPoly<T,Storage>::value_type& a, 
	   const OrthogPoly<T,Storage>& b)
{
  if (a != b.coeff(0))
    return false;
  for (int i=1; i<b.size(); i++)
    if (b.coeff(i) != T(0.0))
      return false;
  return true;
}

template <typename T, typename Storage>
bool
operator==(const OrthogPoly<T,Storage>& a, 
	   const typename OrthogPoly<T,Storage>::value_type& b)
{
  if (a.coeff(0) != b)
    return false;
  for (int i=1; i<a.size(); i++)
    if (a.coeff(i) != T(0.0))
      return false;
  return true;
}

template <typename T, typename Storage>
bool
operator!=(const OrthogPoly<T,Storage>& a, 
	   const OrthogPoly<T,Storage>& b)
{
  return !(a == b);
}

template <typename T, typename Storage>
bool
operator!=(const typename OrthogPoly<T,Storage>::value_type& a, 
	   const OrthogPoly<T,Storage>& b)
{
  return !(a == b);
}

template <typename T, typename Storage>
bool
operator!=(const OrthogPoly<T,Storage>& a, 
	   const typename OrthogPoly<T,Storage>::value_type& b)
{
  return !(a == b);
}

template <typename T, typename Storage>
bool
operator<=(const OrthogPoly<T,Storage>& a, 
	   const OrthogPoly<T,Storage>& b)
{
  return a.two_norm() <= b.two_norm();
}

template <typename T, typename Storage>
bool
operator<=(const typename OrthogPoly<T,Storage>::value_type& a, 
	   const OrthogPoly<T,Storage>& b)
{
  return a <= b.two_norm();
}

template <typename T, typename Storage>
bool
operator<=(const OrthogPoly<T,Storage>& a, 
	   const typename OrthogPoly<T,Storage>::value_type& b)
{
  return a.two_norm() <= b;
}

template <typename T, typename Storage>
bool
operator>=(const OrthogPoly<T,Storage>& a, 
	   const OrthogPoly<T,Storage>& b)
{
  return a.two_norm() >= b.two_norm();
}

template <typename T, typename Storage>
bool
operator>=(const typename OrthogPoly<T,Storage>::value_type& a, 
	   const OrthogPoly<T,Storage>& b)
{
  return a >= b.two_norm();
}

template <typename T, typename Storage>
bool
operator>=(const OrthogPoly<T,Storage>& a, 
	   const typename OrthogPoly<T,Storage>::value_type& b)
{
  return a.two_norm() >= b;
}

template <typename T, typename Storage>
bool
operator<(const OrthogPoly<T,Storage>& a, 
	  const OrthogPoly<T,Storage>& b)
{
  return a.two_norm() < b.two_norm();
}

template <typename T, typename Storage>
bool
operator<(const typename OrthogPoly<T,Storage>::value_type& a, 
	  const OrthogPoly<T,Storage>& b)
{
  return a < b.two_norm();
}

template <typename T, typename Storage>
bool
operator<(const OrthogPoly<T,Storage>& a, 
	  const typename OrthogPoly<T,Storage>::value_type& b)
{
  return a.two_norm() < b;
}

template <typename T, typename Storage>
bool
operator>(const OrthogPoly<T,Storage>& a, 
	  const OrthogPoly<T,Storage>& b)
{
  return a.two_norm() > b.two_norm();
}

template <typename T, typename Storage>
bool
operator>(const typename OrthogPoly<T,Storage>::value_type& a, 
	  const OrthogPoly<T,Storage>& b)
{
  return a > b.two_norm();
}

template <typename T, typename Storage>
bool
operator>(const OrthogPoly<T,Storage>& a, 
	  const typename OrthogPoly<T,Storage>::value_type& b)
{
  return a.two_norm() > b;
}

template <typename T, typename Storage>
bool toBool(const OrthogPoly<T,Storage>& x) {
  bool is_zero = true;
  for (int i=0; i<x.size(); i++)
    is_zero = is_zero && (x.coeff(i) == 0.0);
  return !is_zero;
}

template <typename T, typename Storage>
inline bool
operator && (const OrthogPoly<T,Storage>& x1, const OrthogPoly<T,Storage>& x2)
{
  return toBool(x1) && toBool(x2);
}

template <typename T, typename Storage>
inline bool
operator && (const typename OrthogPoly<T,Storage>::value_type& a, 
	     const OrthogPoly<T,Storage>& x2)
{
  return a && toBool(x2);
}

template <typename T, typename Storage>
inline bool
operator && (const OrthogPoly<T,Storage>& x1, 
	     const typename OrthogPoly<T,Storage>::value_type& b)
{
  return toBool(x1) && b;
}

template <typename T, typename Storage>
inline bool
operator || (const OrthogPoly<T,Storage>& x1, const OrthogPoly<T,Storage>& x2)
{
  return toBool(x1) || toBool(x2);
}

template <typename T, typename Storage>
inline bool
operator || (const typename OrthogPoly<T,Storage>::value_type& a, 
	     const OrthogPoly<T,Storage>& x2)
{
  return a || toBool(x2);
}

template <typename T, typename Storage>
inline bool
operator || (const OrthogPoly<T,Storage>& x1, 
	     const typename OrthogPoly<T,Storage>::value_type& b)
{
  return toBool(x1) || b;
}

template <typename T, typename Storage>
std::ostream& 
operator << (std::ostream& os, const OrthogPoly<T,Storage>& a)
{
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;

  os << "[ ";
      
  for (ordinal_type i=0; i<a.size(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}

template <typename T, typename Storage>
std::istream& 
operator >> (std::istream& is, OrthogPoly<T,Storage>& a)
{
  typedef typename OrthogPoly<T,Storage>::ordinal_type ordinal_type;

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

} // namespace PCE
} // namespace Sacado
