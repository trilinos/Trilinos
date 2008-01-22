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
template <typename T> typename UnivariateHermite<T>::he_type UnivariateHermite<T>::expansion(1);

template <typename T> 
UnivariateHermite<T>::
UnivariateHermite() :
  th(new Stokhos::HermitePoly<T>)
{
}

template <typename T> 
UnivariateHermite<T>::
UnivariateHermite(const T& x) :
  th(new Stokhos::HermitePoly<T>(x))
{
}

template <typename T> 
UnivariateHermite<T>::
UnivariateHermite(unsigned int d, const T& x) :
  th(new Stokhos::HermitePoly<T>(d, x))
{
}

template <typename T> 
UnivariateHermite<T>::
UnivariateHermite(unsigned int d) :
  th(new Stokhos::HermitePoly<T>(d))
{
}

template <typename T> 
UnivariateHermite<T>::
UnivariateHermite(const UnivariateHermite<T>& x) :
  th(x.th)
{
}

template <typename T> 
UnivariateHermite<T>::
~UnivariateHermite()
{
}

template <typename T> 
void
UnivariateHermite<T>::
resize(unsigned int d, bool keep_coeffs)
{
  th->resize(d,keep_coeffs);
}

template <typename T> 
void
UnivariateHermite<T>::
reserve(unsigned int d)
{
  th->reserve(d);
}

template <typename T> 
void
UnivariateHermite<T>::
initExpansion(unsigned int d)
{
  expansion.resize(d+1);
}

template <typename T> 
Stokhos::StandardPoly<T>
UnivariateHermite<T>::
toStandardBasis() const
{
  return expansion.getBasis().toStandardBasis(th->coeff(), th->degree()+1);
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator=(const T& val) 
{
  th.makeOwnCopy();

  if (th->length() < 1) {
    th->resize(0); // degree 0
  }

  (*th)[0] = val;
  Sacado::ds_array<T>::zero(th->coeff()+1, th->degree());

  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator=(const UnivariateHermite<T>& x) 
{
  th = x.th;
  return *this;
}

template <typename T> 
UnivariateHermite<T>
UnivariateHermite<T>::
operator+() const
{
  return *this;
}

template <typename T> 
UnivariateHermite<T> 
UnivariateHermite<T>::
operator-() const
{
  UnivariateHermite<T> x(th->degree());
  expansion.unaryMinus(*(x.th), *th);
  return x;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator+=(const T& val)
{
  th.makeOwnCopy();
  expansion.plusEqual(*th, val);
  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator-=(const T& val)
{
  th.makeOwnCopy();
  expansion.minusEqual(*th, val);
  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator*=(const T& val)
{
  th.makeOwnCopy();
  expansion.timesEqual(*th, val);
  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator/=(const T& val)
{
  th.makeOwnCopy();
  expansion.divideEqual(*th, val);
  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator+=(const UnivariateHermite<T>& x)
{
  th.makeOwnCopy();
  expansion.plusEqual(*th, x);
  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator-=(const UnivariateHermite<T>& x)
{
  th.makeOwnCopy();
  expansion.minusEqual(*th, x);
  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator*=(const UnivariateHermite<T>& x)
{
  th.makeOwnCopy();
  expansion.timesEqual(*th, x);
  return *this;
}

template <typename T> 
UnivariateHermite<T>& 
UnivariateHermite<T>::
operator/=(const UnivariateHermite<T>& x)
{
  th.makeOwnCopy();
  expansion.divideEqual(*th, x);
  return *this;
}

template <typename T>
UnivariateHermite<T>
operator+(const UnivariateHermite<T>& a, 
	  const UnivariateHermite<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  UnivariateHermite<T> c(dc);
  UnivariateHermite<T>::expansion.plus(c.getHermitePoly(), a.getHermitePoly(), 
				       b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator+(const T& a, const UnivariateHermite<T>& b)
{
  UnivariateHermite<T> c(b.degree());
  UnivariateHermite<T>::expansion.plus(c.getHermitePoly(), a, 
				       b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator+(const UnivariateHermite<T>& a, const T& b)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.plus(c.getHermitePoly(), a.getHermitePoly(), 
				       b);
  return c;
}

template <typename T>
UnivariateHermite<T>
operator-(const UnivariateHermite<T>& a, 
	  const UnivariateHermite<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  UnivariateHermite<T> c(dc);
  UnivariateHermite<T>::expansion.minus(c.getHermitePoly(), a.getHermitePoly(),
					b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator-(const T& a, const UnivariateHermite<T>& b)
{
  UnivariateHermite<T> c(b.degree());
  UnivariateHermite<T>::expansion.minus(c.getHermitePoly(), a, 
					b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator-(const UnivariateHermite<T>& a, const T& b)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.minus(c.getHermitePoly(), a.getHermitePoly(),
					b);
  return c;
}

template <typename T>
UnivariateHermite<T>
operator*(const UnivariateHermite<T>& a, 
	  const UnivariateHermite<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  UnivariateHermite<T> c(dc);
  UnivariateHermite<T>::expansion.times(c.getHermitePoly(), 
					a.getHermitePoly(), 
					b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator*(const T& a, const UnivariateHermite<T>& b)
{
  UnivariateHermite<T> c(b.degree());
  UnivariateHermite<T>::expansion.times(c.getHermitePoly(), a, 
					b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator*(const UnivariateHermite<T>& a, const T& b)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.times(c.getHermitePoly(), 
					a.getHermitePoly(), 
					b);
  return c;
}

template <typename T>
UnivariateHermite<T>
operator/(const UnivariateHermite<T>& a, 
	  const UnivariateHermite<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  UnivariateHermite<T> c(dc);
  UnivariateHermite<T>::expansion.divide(c.getHermitePoly(), 
					 a.getHermitePoly(), 
					 b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator/(const T& a, const UnivariateHermite<T>& b)
{
  UnivariateHermite<T> c(b.degree());
  UnivariateHermite<T>::expansion.divide(c.getHermitePoly(), a, 
					 b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
operator/(const UnivariateHermite<T>& a, const T& b)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.divide(c.getHermitePoly(), 
					 a.getHermitePoly(), 
					 b);
  return c;
}

template <typename T>
UnivariateHermite<T>
exp(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.exp(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
log(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.log(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
log10(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.log10(c.getHermitePoly(), 
					a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
sqrt(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.sqrt(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
pow(const UnivariateHermite<T>& a,
    const UnivariateHermite<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  UnivariateHermite<T> c(dc);
  UnivariateHermite<T>::expansion.pow(c.getHermitePoly(), a.getHermitePoly(), 
				      b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
pow(const T& a,
    const UnivariateHermite<T>& b)
{
  UnivariateHermite<T> c(b.degree());
  UnivariateHermite<T>::expansion.pow(c.getHermitePoly(), a, 
				      b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
pow(const UnivariateHermite<T>& a,
    const T& b)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.pow(c.getHermitePoly(), a.getHermitePoly(), 
				      b);
  return c;
}

template <typename T>
UnivariateHermite<T>
sin(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.sin(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
cos(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.cos(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
tan(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.tan(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
sinh(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.sinh(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
cosh(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.cosh(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
tanh(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.tanh(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
acos(const UnivariateHermite<T>& a)
{
   UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.acos(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
asin(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.asin(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
atan(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.atan(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
acosh(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.acosh(c.getHermitePoly(), 
					a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
asinh(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.asinh(c.getHermitePoly(), 
					a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
atanh(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.atanh(c.getHermitePoly(), 
					a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
fabs(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.fabs(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
abs(const UnivariateHermite<T>& a)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.abs(c.getHermitePoly(), a.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
max(const UnivariateHermite<T>& a,
    const UnivariateHermite<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  UnivariateHermite<T> c(dc);
  UnivariateHermite<T>::expansion.max(c.getHermitePoly(), a.getHermitePoly(), 
				      b.getHermitePoly());
}

template <typename T>
UnivariateHermite<T>
max(const T& a,
    const UnivariateHermite<T>& b)
{
  UnivariateHermite<T> c(b.degree());
  UnivariateHermite<T>::expansion.max(c.getHermitePoly(), a, 
				      b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
max(const UnivariateHermite<T>& a,
    const T& b)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.max(c.getHermitePoly(), a.getHermitePoly(), 
				      b);
  return c;
}

template <typename T>
UnivariateHermite<T>
min(const UnivariateHermite<T>& a,
    const UnivariateHermite<T>& b)
{
  unsigned int da = a.degree();
  unsigned int db = b.degree();
  unsigned int dc = da > db ? da : db;
  UnivariateHermite<T> c(dc);
  UnivariateHermite<T>::expansion.min(c.getHermitePoly(), a.getHermitePoly(), 
				      b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
min(const T& a,
    const UnivariateHermite<T>& b)
{
  UnivariateHermite<T> c(b.degree());
  UnivariateHermite<T>::expansion.min(c.getHermitePoly(), a, 
				      b.getHermitePoly());
  return c;
}

template <typename T>
UnivariateHermite<T>
min(const UnivariateHermite<T>& a,
    const T& b)
{
  UnivariateHermite<T> c(a.degree());
  UnivariateHermite<T>::expansion.min(c.getHermitePoly(), a.getHermitePoly(), 
				      b);
  return c;
}

template <typename T>
bool
operator==(const UnivariateHermite<T>& a, 
	   const UnivariateHermite<T>& b)
{
  return a.coeff(0) == b.coeff(0);
}

template <typename T>
bool
operator==(const T& a, 
	   const UnivariateHermite<T>& b)
{
  return a == b.coeff(0);
}

template <typename T>
bool
operator==(const UnivariateHermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) == b;
}

template <typename T>
bool
operator!=(const UnivariateHermite<T>& a, 
	   const UnivariateHermite<T>& b)
{
  return a.coeff(0) != b.coeff(0);
}

template <typename T>
bool
operator!=(const T& a, 
	   const UnivariateHermite<T>& b)
{
  return a != b.coeff(0);
}

template <typename T>
bool
operator!=(const UnivariateHermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) != b;
}

template <typename T>
bool
operator<=(const UnivariateHermite<T>& a, 
	   const UnivariateHermite<T>& b)
{
  return a.coeff(0) <= b.coeff(0);
}

template <typename T>
bool
operator<=(const T& a, 
	   const UnivariateHermite<T>& b)
{
  return a <= b.coeff(0);
}

template <typename T>
bool
operator<=(const UnivariateHermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) <= b;
}

template <typename T>
bool
operator>=(const UnivariateHermite<T>& a, 
	   const UnivariateHermite<T>& b)
{
  return a.coeff(0) >= b.coeff(0);
}

template <typename T>
bool
operator>=(const T& a, 
	   const UnivariateHermite<T>& b)
{
  return a >= b.coeff(0);
}

template <typename T>
bool
operator>=(const UnivariateHermite<T>& a, 
	   const T& b)
{
  return a.coeff(0) >= b;
}

template <typename T>
bool
operator<(const UnivariateHermite<T>& a, 
	  const UnivariateHermite<T>& b)
{
  return a.coeff(0) < b.coeff(0);
}

template <typename T>
bool
operator<(const T& a, 
	  const UnivariateHermite<T>& b)
{
  return a < b.coeff(0);
}

template <typename T>
bool
operator<(const UnivariateHermite<T>& a, 
	  const T& b)
{
  return a.coeff(0) < b;
}

template <typename T>
bool
operator>(const UnivariateHermite<T>& a, 
	  const UnivariateHermite<T>& b)
{
  return a.coeff(0) > b.coeff(0);
}

template <typename T>
bool
operator>(const T& a, 
	  const UnivariateHermite<T>& b)
{
  return a > b.coeff(0);
}

template <typename T>
bool
operator>(const UnivariateHermite<T>& a, 
	  const T& b)
{
  return a.coeff(0) > b;
}

template <typename T>
std::ostream& 
operator << (std::ostream& os, const UnivariateHermite<T>& a)
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
