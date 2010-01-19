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
#include "Stokhos_ConstantOrthogPolyExpansion.hpp"

namespace Sacado {
namespace PCE {

template <typename T> 
OrthogPoly<T>::
OrthogPoly() :
  expansion_(),
  th(new Stokhos::OrthogPolyApprox<int,value_type>)
{
  expansion_ = Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<int,T>);
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(const typename OrthogPoly<T>::value_type& x) :
  expansion_(),
  th(new Stokhos::OrthogPolyApprox<int,value_type>)
{
  expansion_ = Teuchos::rcp(new Stokhos::ConstantOrthogPolyExpansion<int,T>);
  (*th)[0] = x;
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(const Teuchos::RCP<expansion_type>& expansion) :
  expansion_(expansion),
  th(new Stokhos::OrthogPolyApprox<int,value_type>(expansion_->getBasis()))
{
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(const Teuchos::RCP<expansion_type>& expansion,
	   ordinal_type sz) :
  expansion_(expansion),
  th(new Stokhos::OrthogPolyApprox<int,value_type>(expansion_->getBasis(), sz))
{
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(const OrthogPoly<T>& x) :
  expansion_(x.expansion_),
  th(x.th)
{
}

template <typename T> 
OrthogPoly<T>::
~OrthogPoly()
{
}

template <typename T> 
void
OrthogPoly<T>::
reset(const Teuchos::RCP<expansion_type>& expansion)
{
  expansion_ = expansion;
  th->reset(expansion_->getBasis());
}

template <typename T> 
void
OrthogPoly<T>::
reset(const Teuchos::RCP<expansion_type>& expansion, ordinal_type sz)
{
  expansion_ = expansion;
  th->reset(expansion_->getBasis(), sz);
}

template <typename T> 
typename OrthogPoly<T>::value_type
OrthogPoly<T>::
evaluate(const Teuchos::Array<typename OrthogPoly<T>::value_type>& point) const
{
  return th->evaluate(point);
}

template <typename T> 
typename OrthogPoly<T>::value_type
OrthogPoly<T>::
evaluate(const Teuchos::Array<typename OrthogPoly<T>::value_type>& point,
         const Teuchos::Array<typename OrthogPoly<T>::value_type>& bvals) const
{
  return th->evaluate(point, bvals);
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator=(const typename OrthogPoly<T>::value_type& v) 
{
  th.makeOwnCopy();

  (*th)[0] = v;
  Sacado::ds_array<value_type>::zero(th->coeff()+1, th->size()-1);

  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator=(const OrthogPoly<T>& x) 
{
  expansion_ = x.expansion_;
  th = x.th;
  return *this;
}

template <typename T> 
OrthogPoly<T>
OrthogPoly<T>::
operator+() const
{
  return *this;
}

template <typename T> 
OrthogPoly<T> 
OrthogPoly<T>::
operator-() const
{
  OrthogPoly<T> x(expansion_);
  expansion_->unaryMinus(*(x.th), *th);
  return x;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator+=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->plusEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator-=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->minusEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator*=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->timesEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator/=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion_->divideEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator+=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->plusEqual(*th, *x.th);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator-=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->minusEqual(*th, *x.th);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator*=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->timesEqual(*th, *x.th);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator/=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = expansion_;
  if (x.size() > size()) {
    e = x.expansion();
    reset(e, size());
  }
  e->divideEqual(*th, *x.th);
  return *this;
}

template <typename T>
OrthogPoly<T>
operator+(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T> c(e, 0);
  e->plus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	  b.getOrthogPolyApprox());

  return c;
}

template <typename T>
OrthogPoly<T>
operator+(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.expansion(), 0);
  b.expansion()->plus(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator+(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->plus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T>
OrthogPoly<T>
operator-(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T> c(e, 0);
  e->minus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	   b.getOrthogPolyApprox());

  return c;
}

template <typename T>
OrthogPoly<T>
operator-(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.expansion(), 0);
  b.expansion()->minus(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator-(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->minus(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T>
OrthogPoly<T>
operator*(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T> c(e, 0);
  e->times(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	   b.getOrthogPolyApprox());

  return c;
}

template <typename T>
OrthogPoly<T>
operator*(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.expansion(), 0);
  b.expansion()->times(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator*(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->times(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T>
OrthogPoly<T>
operator/(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T> c(e, 0);
  e->divide(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	    b.getOrthogPolyApprox());

  return c;
}

template <typename T>
OrthogPoly<T>
operator/(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.expansion(), 0);
  b.expansion()->divide(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator/(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->divide(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T>
OrthogPoly<T>
exp(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->exp(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
log(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->log(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
log10(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->log10(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
sqrt(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->sqrt(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
pow(const OrthogPoly<T>& a, 
    const OrthogPoly<T>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T> c(e, 0);
  e->pow(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	 b.getOrthogPolyApprox());

  return c;
}

template <typename T>
OrthogPoly<T>
pow(const T& a,
    const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.expansion(), 0);
  b.expansion()->pow(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
pow(const OrthogPoly<T>& a,
    const T& b)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->pow(c.getOrthogPolyApprox(),a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T>
OrthogPoly<T>
sin(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->sin(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
cos(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->cos(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
tan(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->tan(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
sinh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->sinh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
cosh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->cosh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
tanh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->tanh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
acos(const OrthogPoly<T>& a)
{
   OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->acos(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
asin(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->asin(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
atan(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->atan(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
acosh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->acosh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
asinh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->asinh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
atanh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->atanh(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
fabs(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->fabs(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
abs(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->abs(c.getOrthogPolyApprox(), a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
max(const OrthogPoly<T>& a,
    const OrthogPoly<T>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T> c(e, 0);
  e->max(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	 b.getOrthogPolyApprox());

  return c;
}

template <typename T>
OrthogPoly<T>
max(const typename OrthogPoly<T>::value_type& a,
    const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.expansion(), 0);
  b.expansion()->max(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
max(const OrthogPoly<T>& a,
    const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->max(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T>
OrthogPoly<T>
min(const OrthogPoly<T>& a,
    const OrthogPoly<T>& b)
{
  // Get expansion
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;
  ordinal_type da = a.size();
  ordinal_type db = b.size();
  Teuchos::RCP<typename OrthogPoly<T>::expansion_type> e = a.expansion();
  if (da == db || da > 1)
    e = a.expansion();
  else
    e = b.expansion();

  OrthogPoly<T> c(e, 0);
  e->min(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), 
	 b.getOrthogPolyApprox());

  return c;
}

template <typename T>
OrthogPoly<T>
min(const typename OrthogPoly<T>::value_type& a,
    const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.expansion(), 0);
  b.expansion()->min(c.getOrthogPolyApprox(), a, b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
min(const OrthogPoly<T>& a,
    const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.expansion(), 0);
  a.expansion()->min(c.getOrthogPolyApprox(), a.getOrthogPolyApprox(), b);
  return c;
}

template <typename T>
bool
operator==(const OrthogPoly<T>& a, 
	   const OrthogPoly<T>& b)
{
  return a.coeff(0) == b.coeff(0);
}

template <typename T>
bool
operator==(const typename OrthogPoly<T>::value_type& a, 
	   const OrthogPoly<T>& b)
{
  return a == b.coeff(0);
}

template <typename T>
bool
operator==(const OrthogPoly<T>& a, 
	   const typename OrthogPoly<T>::value_type& b)
{
  return a.coeff(0) == b;
}

template <typename T>
bool
operator!=(const OrthogPoly<T>& a, 
	   const OrthogPoly<T>& b)
{
  return a.coeff(0) != b.coeff(0);
}

template <typename T>
bool
operator!=(const typename OrthogPoly<T>::value_type& a, 
	   const OrthogPoly<T>& b)
{
  return a != b.coeff(0);
}

template <typename T>
bool
operator!=(const OrthogPoly<T>& a, 
	   const typename OrthogPoly<T>::value_type& b)
{
  return a.coeff(0) != b;
}

template <typename T>
bool
operator<=(const OrthogPoly<T>& a, 
	   const OrthogPoly<T>& b)
{
  return a.coeff(0) <= b.coeff(0);
}

template <typename T>
bool
operator<=(const typename OrthogPoly<T>::value_type& a, 
	   const OrthogPoly<T>& b)
{
  return a <= b.coeff(0);
}

template <typename T>
bool
operator<=(const OrthogPoly<T>& a, 
	   const typename OrthogPoly<T>::value_type& b)
{
  return a.coeff(0) <= b;
}

template <typename T>
bool
operator>=(const OrthogPoly<T>& a, 
	   const OrthogPoly<T>& b)
{
  return a.coeff(0) >= b.coeff(0);
}

template <typename T>
bool
operator>=(const typename OrthogPoly<T>::value_type& a, 
	   const OrthogPoly<T>& b)
{
  return a >= b.coeff(0);
}

template <typename T>
bool
operator>=(const OrthogPoly<T>& a, 
	   const typename OrthogPoly<T>::value_type& b)
{
  return a.coeff(0) >= b;
}

template <typename T>
bool
operator<(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  return a.coeff(0) < b.coeff(0);
}

template <typename T>
bool
operator<(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  return a < b.coeff(0);
}

template <typename T>
bool
operator<(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  return a.coeff(0) < b;
}

template <typename T>
bool
operator>(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  return a.coeff(0) > b.coeff(0);
}

template <typename T>
bool
operator>(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  return a > b.coeff(0);
}

template <typename T>
bool
operator>(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  return a.coeff(0) > b;
}

template <typename T>
std::ostream& 
operator << (std::ostream& os, const OrthogPoly<T>& a)
{
  typedef typename OrthogPoly<T>::ordinal_type ordinal_type;

  os << "[ ";
      
  for (ordinal_type i=0; i<a.size(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}

} // namespace PCE
} // namespace Sacado
