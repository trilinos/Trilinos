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
template <typename BasisT> 
Teuchos::RCP< Stokhos::OrthogPolyExpansion<BasisT> > 
OrthogPoly<BasisT>::expansion = Teuchos::null;

template <typename BasisT> 
OrthogPoly<BasisT>::
OrthogPoly() :
  th(new Stokhos::OrthogPolyApprox<value_type>)
{
}

template <typename BasisT> 
OrthogPoly<BasisT>::
OrthogPoly(const typename OrthogPoly<BasisT>::value_type& x) :
  th(new Stokhos::OrthogPolyApprox<value_type>(x))
{
}

template <typename BasisT> 
OrthogPoly<BasisT>::
OrthogPoly(unsigned int sz, const typename OrthogPoly<BasisT>::value_type& x) :
  th(new Stokhos::OrthogPolyApprox<value_type>(sz, x))
{
}

template <typename BasisT> 
OrthogPoly<BasisT>::
OrthogPoly(unsigned int sz) :
  th(new Stokhos::OrthogPolyApprox<value_type>(sz))
{
}

template <typename BasisT> 
OrthogPoly<BasisT>::
OrthogPoly(const OrthogPoly<BasisT>& x) :
  th(x.th)
{
}

template <typename BasisT> 
OrthogPoly<BasisT>::
~OrthogPoly()
{
}

template <typename BasisT> 
void
OrthogPoly<BasisT>::
resize(unsigned int sz)
{
  th->resize(sz);
}

template <typename BasisT> 
void
OrthogPoly<BasisT>::
reserve(unsigned int sz)
{
  th->reserve(sz);
}

template <typename BasisT> 
void
OrthogPoly<BasisT>::
initExpansion(const Teuchos::RCP< Stokhos::OrthogPolyExpansion<BasisT> >& e)
{
  expansion = e;
}

template <typename BasisT> 
Stokhos::Polynomial<typename OrthogPoly<BasisT>::value_type>
OrthogPoly<BasisT>::
toStandardBasis() const
{
  return expansion->getBasis().toStandardBasis(th->coeff(), th->size());
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator=(const typename OrthogPoly<BasisT>::value_type& val) 
{
  th.makeOwnCopy();

  if (th->length() < 1) {
    th->resize(1);
  }

  (*th)[0] = val;
  Sacado::ds_array<value_type>::zero(th->coeff()+1, th->size()-1);

  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator=(const OrthogPoly<BasisT>& x) 
{
  th = x.th;
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>
OrthogPoly<BasisT>::
operator+() const
{
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT> 
OrthogPoly<BasisT>::
operator-() const
{
  OrthogPoly<BasisT> x(th->size());
  expansion->unaryMinus(*(x.th), *th);
  return x;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator+=(const typename OrthogPoly<BasisT>::value_type& val)
{
  th.makeOwnCopy();
  expansion->plusEqual(*th, val);
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator-=(const typename OrthogPoly<BasisT>::value_type& val)
{
  th.makeOwnCopy();
  expansion->minusEqual(*th, val);
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator*=(const typename OrthogPoly<BasisT>::value_type& val)
{
  th.makeOwnCopy();
  expansion->timesEqual(*th, val);
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator/=(const typename OrthogPoly<BasisT>::value_type& val)
{
  th.makeOwnCopy();
  expansion->divideEqual(*th, val);
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator+=(const OrthogPoly<BasisT>& x)
{
  th.makeOwnCopy();
  expansion->plusEqual(*th, x);
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator-=(const OrthogPoly<BasisT>& x)
{
  th.makeOwnCopy();
  expansion->minusEqual(*th, x);
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator*=(const OrthogPoly<BasisT>& x)
{
  th.makeOwnCopy();
  expansion->timesEqual(*th, x);
  return *this;
}

template <typename BasisT> 
OrthogPoly<BasisT>& 
OrthogPoly<BasisT>::
operator/=(const OrthogPoly<BasisT>& x)
{
  th.makeOwnCopy();
  expansion->divideEqual(*th, x);
  return *this;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator+(const OrthogPoly<BasisT>& a, 
	  const OrthogPoly<BasisT>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<BasisT> c(dc);
  OrthogPoly<BasisT>::expansion->plus(c.getOrthogPolyApprox(), 
				      a.getOrthogPolyApprox(), 
				      b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator+(const typename OrthogPoly<BasisT>::value_type& a, 
	  const OrthogPoly<BasisT>& b)
{
  OrthogPoly<BasisT> c(b.size());
  OrthogPoly<BasisT>::expansion->plus(c.getOrthogPolyApprox(), a, 
				b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator+(const OrthogPoly<BasisT>& a, 
	  const typename OrthogPoly<BasisT>::value_type& b)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->plus(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox(), 
				b);
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator-(const OrthogPoly<BasisT>& a, 
	  const OrthogPoly<BasisT>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<BasisT> c(dc);
  OrthogPoly<BasisT>::expansion->minus(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(),
				 b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator-(const typename OrthogPoly<BasisT>::value_type& a, 
	  const OrthogPoly<BasisT>& b)
{
  OrthogPoly<BasisT> c(b.size());
  OrthogPoly<BasisT>::expansion->minus(c.getOrthogPolyApprox(), a, 
				 b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator-(const OrthogPoly<BasisT>& a, 
	  const typename OrthogPoly<BasisT>::value_type& b)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->minus(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(),
				 b);
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator*(const OrthogPoly<BasisT>& a, 
	  const OrthogPoly<BasisT>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<BasisT> c(dc);
  OrthogPoly<BasisT>::expansion->times(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(), 
				 b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator*(const typename OrthogPoly<BasisT>::value_type& a, 
	  const OrthogPoly<BasisT>& b)
{
  OrthogPoly<BasisT> c(b.size());
  OrthogPoly<BasisT>::expansion->times(c.getOrthogPolyApprox(), a, 
				 b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator*(const OrthogPoly<BasisT>& a, 
	  const typename OrthogPoly<BasisT>::value_type& b)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->times(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox(), 
				 b);
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator/(const OrthogPoly<BasisT>& a, 
	  const OrthogPoly<BasisT>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<BasisT> c(dc);
  OrthogPoly<BasisT>::expansion->divide(c.getOrthogPolyApprox(), 
				  a.getOrthogPolyApprox(), 
				  b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator/(const typename OrthogPoly<BasisT>::value_type& a, 
	  const OrthogPoly<BasisT>& b)
{
  OrthogPoly<BasisT> c(b.size());
  OrthogPoly<BasisT>::expansion->divide(c.getOrthogPolyApprox(), a, 
				  b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
operator/(const OrthogPoly<BasisT>& a, 
	  const typename OrthogPoly<BasisT>::value_type& b)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->divide(c.getOrthogPolyApprox(), 
				  a.getOrthogPolyApprox(), 
				  b);
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
exp(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->exp(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
log(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->log(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
log10(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->log10(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
sqrt(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->sqrt(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
pow(const OrthogPoly<BasisT>& a,
    const OrthogPoly<BasisT>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<BasisT> c(dc);
  OrthogPoly<BasisT>::expansion->pow(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
pow(const typename OrthogPoly<BasisT>::value_type& a,
    const OrthogPoly<BasisT>& b)
{
  OrthogPoly<BasisT> c(b.size());
  OrthogPoly<BasisT>::expansion->pow(c.getOrthogPolyApprox(), a, 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
pow(const OrthogPoly<BasisT>& a,
    const typename OrthogPoly<BasisT>::value_type& b)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->pow(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b);
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
sin(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->sin(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
cos(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->cos(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
tan(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->tan(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
sinh(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->sinh(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
cosh(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->cosh(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
tanh(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->tanh(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
acos(const OrthogPoly<BasisT>& a)
{
   OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->acos(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
asin(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->asin(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
atan(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->atan(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
acosh(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->acosh(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
asinh(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->asinh(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
atanh(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->atanh(c.getOrthogPolyApprox(), 
				 a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
fabs(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->fabs(c.getOrthogPolyApprox(), 
				a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
abs(const OrthogPoly<BasisT>& a)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->abs(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
max(const OrthogPoly<BasisT>& a,
    const OrthogPoly<BasisT>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<BasisT> c(dc);
  OrthogPoly<BasisT>::expansion->max(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b.getOrthogPolyApprox());
}

template <typename BasisT>
OrthogPoly<BasisT>
max(const typename OrthogPoly<BasisT>::value_type& a,
    const OrthogPoly<BasisT>& b)
{
  OrthogPoly<BasisT> c(b.size());
  OrthogPoly<BasisT>::expansion->max(c.getOrthogPolyApprox(), a, 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
max(const OrthogPoly<BasisT>& a,
    const typename OrthogPoly<BasisT>::value_type& b)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->max(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b);
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
min(const OrthogPoly<BasisT>& a,
    const OrthogPoly<BasisT>& b)
{
  unsigned int da = a.size();
  unsigned int db = b.size();
  unsigned int dc = da > db ? da : db;
  OrthogPoly<BasisT> c(dc);
  OrthogPoly<BasisT>::expansion->min(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
min(const typename OrthogPoly<BasisT>::value_type& a,
    const OrthogPoly<BasisT>& b)
{
  OrthogPoly<BasisT> c(b.size());
  OrthogPoly<BasisT>::expansion->min(c.getOrthogPolyApprox(), a, 
			       b.getOrthogPolyApprox());
  return c;
}

template <typename BasisT>
OrthogPoly<BasisT>
min(const OrthogPoly<BasisT>& a,
    const typename OrthogPoly<BasisT>::value_type& b)
{
  OrthogPoly<BasisT> c(a.size());
  OrthogPoly<BasisT>::expansion->min(c.getOrthogPolyApprox(), 
			       a.getOrthogPolyApprox(), 
			       b);
  return c;
}

template <typename BasisT>
bool
operator==(const OrthogPoly<BasisT>& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a.coeff(0) == b.coeff(0);
}

template <typename BasisT>
bool
operator==(const typename OrthogPoly<BasisT>::value_type& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a == b.coeff(0);
}

template <typename BasisT>
bool
operator==(const OrthogPoly<BasisT>& a, 
	   const typename OrthogPoly<BasisT>::value_type& b)
{
  return a.coeff(0) == b;
}

template <typename BasisT>
bool
operator!=(const OrthogPoly<BasisT>& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a.coeff(0) != b.coeff(0);
}

template <typename BasisT>
bool
operator!=(const typename OrthogPoly<BasisT>::value_type& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a != b.coeff(0);
}

template <typename BasisT>
bool
operator!=(const OrthogPoly<BasisT>& a, 
	   const typename OrthogPoly<BasisT>::value_type& b)
{
  return a.coeff(0) != b;
}

template <typename BasisT>
bool
operator<=(const OrthogPoly<BasisT>& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a.coeff(0) <= b.coeff(0);
}

template <typename BasisT>
bool
operator<=(const typename OrthogPoly<BasisT>::value_type& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a <= b.coeff(0);
}

template <typename BasisT>
bool
operator<=(const OrthogPoly<BasisT>& a, 
	   const typename OrthogPoly<BasisT>::value_type& b)
{
  return a.coeff(0) <= b;
}

template <typename BasisT>
bool
operator>=(const OrthogPoly<BasisT>& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a.coeff(0) >= b.coeff(0);
}

template <typename BasisT>
bool
operator>=(const typename OrthogPoly<BasisT>::value_type& a, 
	   const OrthogPoly<BasisT>& b)
{
  return a >= b.coeff(0);
}

template <typename BasisT>
bool
operator>=(const OrthogPoly<BasisT>& a, 
	   const typename OrthogPoly<BasisT>::value_type& b)
{
  return a.coeff(0) >= b;
}

template <typename BasisT>
bool
operator<(const OrthogPoly<BasisT>& a, 
	  const OrthogPoly<BasisT>& b)
{
  return a.coeff(0) < b.coeff(0);
}

template <typename BasisT>
bool
operator<(const typename OrthogPoly<BasisT>::value_type& a, 
	  const OrthogPoly<BasisT>& b)
{
  return a < b.coeff(0);
}

template <typename BasisT>
bool
operator<(const OrthogPoly<BasisT>& a, 
	  const typename OrthogPoly<BasisT>::value_type& b)
{
  return a.coeff(0) < b;
}

template <typename BasisT>
bool
operator>(const OrthogPoly<BasisT>& a, 
	  const OrthogPoly<BasisT>& b)
{
  return a.coeff(0) > b.coeff(0);
}

template <typename BasisT>
bool
operator>(const typename OrthogPoly<BasisT>::value_type& a, 
	  const OrthogPoly<BasisT>& b)
{
  return a > b.coeff(0);
}

template <typename BasisT>
bool
operator>(const OrthogPoly<BasisT>& a, 
	  const typename OrthogPoly<BasisT>::value_type& b)
{
  return a.coeff(0) > b;
}

template <typename BasisT>
std::ostream& 
operator << (std::ostream& os, const OrthogPoly<BasisT>& a)
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
