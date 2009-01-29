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

namespace Sacado {
namespace PCE {

// Initialize static data
template <typename T> 
Teuchos::RCP<typename OrthogPoly<T>::expansion_type> 
OrthogPoly<T>::expansion = Teuchos::null;

template <typename T> 
OrthogPoly<T>::
OrthogPoly() :
  th(new Stokhos::OrthogPolyApprox<value_type>)
{
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(const typename OrthogPoly<T>::value_type& x) :
  th(new Stokhos::OrthogPolyApprox<value_type>(x))
{
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(unsigned int sz, const typename OrthogPoly<T>::value_type& x) :
  th(new Stokhos::OrthogPolyApprox<value_type>(sz, x))
{
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(unsigned int sz) :
  th(new Stokhos::OrthogPolyApprox<value_type>(sz))
{
}

template <typename T> 
OrthogPoly<T>::
OrthogPoly(const OrthogPoly<T>& x) :
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
resize(unsigned int sz)
{
  th->resize(sz);
}

template <typename T> 
void
OrthogPoly<T>::
reserve(unsigned int sz)
{
  th->reserve(sz);
}

template <typename T> 
void
OrthogPoly<T>::
initExpansion(const Teuchos::RCP<typename OrthogPoly<T>::expansion_type>& e)
{
  expansion = e;
}

template <typename T> 
Stokhos::Polynomial<typename OrthogPoly<T>::value_type>
OrthogPoly<T>::
toStandardBasis() const
{
  return expansion->getBasis().toStandardBasis(th->coeff(), th->size());
}

template <typename T> 
typename OrthogPoly<T>::value_type
OrthogPoly<T>::
evaluate(const std::vector<typename OrthogPoly<T>::value_type>& point) const
{
  return th->evaluate(expansion->getBasis(), point);
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator=(const typename OrthogPoly<T>::value_type& v) 
{
  th.makeOwnCopy();

  if (th->size() < 1) {
    th->resize(1);
  }

  (*th)[0] = v;
  Sacado::ds_array<value_type>::zero(th->coeff()+1, th->size()-1);

  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator=(const OrthogPoly<T>& x) 
{
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
  OrthogPoly<T> x(th->size());
  expansion->unaryMinus(*(x.th), *th);
  return x;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator+=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion->plusEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator-=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion->minusEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator*=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion->timesEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator/=(const typename OrthogPoly<T>::value_type& v)
{
  th.makeOwnCopy();
  expansion->divideEqual(*th, v);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator+=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  expansion->plusEqual(*th, *x.th);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator-=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  expansion->minusEqual(*th, *x.th);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator*=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  expansion->timesEqual(*th, *x.th);
  return *this;
}

template <typename T> 
OrthogPoly<T>& 
OrthogPoly<T>::
operator/=(const OrthogPoly<T>& x)
{
  th.makeOwnCopy();
  expansion->divideEqual(*th, *x.th);
  return *this;
}

template <typename T>
OrthogPoly<T>
operator+(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<T> c(dc);
  OrthogPoly<T>::expansion->plus(c.getOrthogPolyApprox(), 
				      a.getOrthogPolyApprox(), 
				      b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator+(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.size());
  OrthogPoly<T>::expansion->plus(c.getOrthogPolyApprox(), a, 
				b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator+(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->plus(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox(), 
				b);
  return c;
}

template <typename T>
OrthogPoly<T>
operator-(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<T> c(dc);
  OrthogPoly<T>::expansion->minus(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(),
				 b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator-(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.size());
  OrthogPoly<T>::expansion->minus(c.getOrthogPolyApprox(), a, 
				 b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator-(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->minus(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(),
				 b);
  return c;
}

template <typename T>
OrthogPoly<T>
operator*(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<T> c(dc);
  OrthogPoly<T>::expansion->times(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(), 
				 b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator*(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.size());
  OrthogPoly<T>::expansion->times(c.getOrthogPolyApprox(), a, 
				 b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator*(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->times(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(), 
				 b);
  return c;
}

template <typename T>
OrthogPoly<T>
operator/(const OrthogPoly<T>& a, 
	  const OrthogPoly<T>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<T> c(dc);
  OrthogPoly<T>::expansion->divide(c.getOrthogPolyApprox(), 
				  a.getOrthogPolyApprox(), 
				  b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator/(const typename OrthogPoly<T>::value_type& a, 
	  const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.size());
  OrthogPoly<T>::expansion->divide(c.getOrthogPolyApprox(), a, 
				  b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
operator/(const OrthogPoly<T>& a, 
	  const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->divide(c.getOrthogPolyApprox(), 
				  a.getOrthogPolyApprox(), 
				  b);
  return c;
}

template <typename T>
OrthogPoly<T>
exp(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->exp(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
log(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->log(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
log10(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->log10(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
sqrt(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->sqrt(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
pow(const OrthogPoly<T>& a,
    const OrthogPoly<T>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<T> c(dc);
  OrthogPoly<T>::expansion->pow(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
pow(const T& a,
    const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.size());
  OrthogPoly<T>::expansion->pow(c.getOrthogPolyApprox(), a, 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
pow(const OrthogPoly<T>& a,
    const T& b)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->pow(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b);
  return c;
}

template <typename T>
OrthogPoly<T>
sin(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->sin(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
cos(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->cos(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
tan(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->tan(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
sinh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->sinh(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
cosh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->cosh(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
tanh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->tanh(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
acos(const OrthogPoly<T>& a)
{
   OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->acos(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
asin(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->asin(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
atan(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->atan(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
acosh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->acosh(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
asinh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->asinh(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
atanh(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->atanh(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
fabs(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->fabs(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
abs(const OrthogPoly<T>& a)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->abs(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
max(const OrthogPoly<T>& a,
    const OrthogPoly<T>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<T> c(dc);
  OrthogPoly<T>::expansion->max(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b.getOrthogPolyApprox());
}

template <typename T>
OrthogPoly<T>
max(const typename OrthogPoly<T>::value_type& a,
    const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.size());
  OrthogPoly<T>::expansion->max(c.getOrthogPolyApprox(), a, 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
max(const OrthogPoly<T>& a,
    const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->max(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b);
  return c;
}

template <typename T>
OrthogPoly<T>
min(const OrthogPoly<T>& a,
    const OrthogPoly<T>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<T> c(dc);
  OrthogPoly<T>::expansion->min(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
min(const typename OrthogPoly<T>::value_type& a,
    const OrthogPoly<T>& b)
{
  OrthogPoly<T> c(b.size());
  OrthogPoly<T>::expansion->min(c.getOrthogPolyApprox(), a, 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename T>
OrthogPoly<T>
min(const OrthogPoly<T>& a,
    const typename OrthogPoly<T>::value_type& b)
{
  OrthogPoly<T> c(a.size());
  OrthogPoly<T>::expansion->min(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b);
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
  os << "[ ";
      
  for (unsigned int i=0; i<a.size(); i++) {
    os << a.coeff(i) << " ";
  }

  os << "]\n";
  return os;
}

} // namespace PCE
} // namespace Sacado
