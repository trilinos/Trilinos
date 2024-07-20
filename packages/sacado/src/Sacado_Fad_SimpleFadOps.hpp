// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_SIMPLEFADOPS_HPP
#define SACADO_FAD_SIMPLEFADOPS_HPP

#include "Sacado_cmath.hpp"
#include <ostream>	// for std::ostream

namespace Sacado {

  namespace Fad {

    template <typename ValueT>
    SimpleFad<ValueT>
    operator + (const SimpleFad<ValueT>& a) {
      return a;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator - (const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, -a.val(), -1.0);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    exp(const SimpleFad<ValueT>& a) {
      ValueT t1 = std::exp(a.val());
      return SimpleFad<ValueT>(a, t1, t1);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    log(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::log(a.val()), 1.0/a.val());
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    log10(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::log10(a.val()), 
			       1.0/(std::log(10.0)*a.val()));
    }
    
    template <typename ValueT>
    SimpleFad<ValueT>
    sqrt(const SimpleFad<ValueT>& a) {
      ValueT t1 = std::sqrt(a.val());
      ValueT t2 = 1.0/(2.0*t1);
      return SimpleFad<ValueT>(a, t1, t2);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    cbrt(const SimpleFad<ValueT>& a) {
      ValueT t1 = std::cbrt(a.val());
      ValueT t2 = 1.0/(3.0*t1*t1);
      return SimpleFad<ValueT>(a, t1, t2);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    cos(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::cos(a.val()), -std::sin(a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    sin(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::sin(a.val()), std::cos(a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    tan(const SimpleFad<ValueT>& a) {
      ValueT t1 = std::tan(a.val());
      ValueT t2 = 1.0 + t1*t1;
      return SimpleFad<ValueT>(a, t1, t2);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    acos(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::acos(a.val()), 
			       -1.0/std::sqrt(1.0 - a.val()*a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    asin(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::asin(a.val()), 
			       1.0/std::sqrt(1.0 - a.val()*a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    atan(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::atan(a.val()), 
			       1.0/(1.0 + a.val()*a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    cosh(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::cosh(a.val()), std::sinh(a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    sinh(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::sinh(a.val()), std::cosh(a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    tanh(const SimpleFad<ValueT>& a) {
      ValueT t = std::tanh(a.val());
      return SimpleFad<ValueT>(a, t, 1.0-t*t);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    acosh(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::acosh(a.val()), 
			       1.0/std::sqrt(a.val()*a.val()-1.0));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    asinh(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::asinh(a.val()), 
			       1.0/std::sqrt(1.0 + a.val()*a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    atanh(const SimpleFad<ValueT>& a) {
      return SimpleFad<ValueT>(a, std::atanh(a.val()), 
			       1.0 /(1.0 - a.val()*a.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    abs(const SimpleFad<ValueT>& a) {
      ValueT t = 1.0;
      if (a.val() < 0)
	t = -1.0;
      return SimpleFad<ValueT>(a, std::abs(a.val()), t);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    fabs(const SimpleFad<ValueT>& a) {
      ValueT t = 1.0;
      if (a.val() < 0)
	t = -1.0;
      return SimpleFad<ValueT>(a, std::fabs(a.val()), t);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator + (const SimpleFad<ValueT>& a, 
		const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, a.val() + b.val());
      if (a.hasFastAccess() && b.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i) + b.fastAccessDx(i);
      else if (a.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i);
      else if (b.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = b.fastAccessDx(i);

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator + (const typename SimpleFad<ValueT>::value_type& a, 
		const SimpleFad<ValueT>& b) {
      return SimpleFad<ValueT>(b, a+b.val(), 1.0);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator + (const SimpleFad<ValueT>& a, 
		const typename SimpleFad<ValueT>::value_type& b) {
      return SimpleFad<ValueT>(a, a.val()+b, 1.0);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator - (const SimpleFad<ValueT>& a, 
		const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, a.val() - b.val());
      if (a.hasFastAccess() && b.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i) - b.fastAccessDx(i);
      else if (a.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i);
      else if (b.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = -b.fastAccessDx(i);

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator - (const typename SimpleFad<ValueT>::value_type& a, 
		const SimpleFad<ValueT>& b) {
      return SimpleFad<ValueT>(b, a-b.val(), -1.0);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator - (const SimpleFad<ValueT>& a, 
		const typename SimpleFad<ValueT>::value_type& b) {
      return SimpleFad<ValueT>(a, a.val()-b, 1.0);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator * (const SimpleFad<ValueT>& a, 
		const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, a.val() * b.val());
      if (a.hasFastAccess() && b.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = 
	    a.fastAccessDx(i)*b.val() + b.fastAccessDx(i)*a.val();
      else if (a.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i)*b.val();
      else if (b.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = b.fastAccessDx(i)*a.val();

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator * (const typename SimpleFad<ValueT>::value_type& a, 
		const SimpleFad<ValueT>& b) {
      return SimpleFad<ValueT>(b, a*b.val(), a);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator * (const SimpleFad<ValueT>& a, 
		const typename SimpleFad<ValueT>::value_type& b) {
      return SimpleFad<ValueT>(a, a.val()*b, b);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator / (const SimpleFad<ValueT>& a, 
		const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, a.val() / b.val());
      if (a.hasFastAccess() && b.hasFastAccess()) {
	ValueT t = b.val()*b.val();
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = 
	    (a.fastAccessDx(i)*b.val() - b.fastAccessDx(i)*a.val()) / t;
      }
      else if (a.hasFastAccess())
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i)/b.val();
      else if (b.hasFastAccess()) {
	ValueT t = -a.val()/(b.val()*b.val());
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = b.fastAccessDx(i)*t;
      }

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator / (const typename SimpleFad<ValueT>::value_type& a, 
		const SimpleFad<ValueT>& b) {
      return SimpleFad<ValueT>(b, a/b.val(), -a/(b.val()*b.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    operator / (const SimpleFad<ValueT>& a, 
		const typename SimpleFad<ValueT>::value_type& b) {
      return SimpleFad<ValueT>(a, a.val()/b, 1.0/b);
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    pow(const SimpleFad<ValueT>& a, 
	const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, std::pow(a.val(), b.val()));
      typedef typename SimpleFad<ValueT>::value_type value_type;
      if (a.hasFastAccess() && b.hasFastAccess()) {
        if (a.val() != value_type(0)) {
	  ValueT t1 = c.val()*b.val()/a.val();
	  ValueT t2 = c.val()*std::log(a.val());
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = 
	      a.fastAccessDx(i)*t1 + b.fastAccessDx(i)*t2;
	}
      }
      else if (a.hasFastAccess()) {
        if (b.val() == value_type(1)) {
          for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = a.fastAccessDx(i);
        }
        else if (a.val() != value_type(0)) {
	  ValueT t1 = c.val()*b.val()/a.val();
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = a.fastAccessDx(i)*t1;
	}
      }
      else if (b.hasFastAccess()) {
	if (a.val() != value_type(0)) {
	  ValueT t2 = c.val()*std::log(a.val());
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = b.fastAccessDx(i)*t2;
	}
      }

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    pow(const typename SimpleFad<ValueT>::value_type& a, 
	const SimpleFad<ValueT>& b) {
      typedef typename SimpleFad<ValueT>::value_type value_type;
      ValueT t = std::pow(a,b.val());
      if (a != value_type(0))
	return SimpleFad<ValueT>(b, t, t*std::log(a));
      else
	return SimpleFad<ValueT>(b, t, value_type(0));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    pow(const SimpleFad<ValueT>& a, 
	const typename SimpleFad<ValueT>::value_type& b) {
      typedef typename SimpleFad<ValueT>::value_type value_type;
      ValueT t = std::pow(a.val(),b);
      if (b == value_type(1))
        return a;
      else if (a.val() != value_type(0))
	return SimpleFad<ValueT>(a, t, t*b/a.val());
      else
	return SimpleFad<ValueT>(a, t, value_type(0));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    atan2(const SimpleFad<ValueT>& a, 
	  const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, std::atan2(a.val(), b.val()));
      if (a.hasFastAccess() && b.hasFastAccess()) {
	ValueT t = a.val()*a.val() + b.val()*b.val();
	ValueT t1 = b.val()/t;
	ValueT t2 = a.val()/t;
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = 
	    a.fastAccessDx(i)*t1 - b.fastAccessDx(i)*t2;
      }
      else if (a.hasFastAccess()) {
	ValueT t1 = b.val()/(a.val()*a.val() + b.val()*b.val());
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i)*t1;
      }
      else if (b.hasFastAccess()) {
	ValueT t2 = -a.val()/(a.val()*a.val() + b.val()*b.val());
	for (int i=0; i<sz; i++)
	  c.fastAccessDx(i) = b.fastAccessDx(i)*t2;
      }

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    atan2(const typename SimpleFad<ValueT>::value_type& a, 
	  const SimpleFad<ValueT>& b) {
      return SimpleFad<ValueT>(b, std::atan2(a,b.val()), 
			       -a/(a*a + b.val()*b.val()));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    atan2(const SimpleFad<ValueT>& a, 
	  const typename SimpleFad<ValueT>::value_type& b) {
      return SimpleFad<ValueT>(a, std::atan2(a.val(),b), 
			       b/(a.val()*a.val() + b*b));
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    max(const SimpleFad<ValueT>& a, 
	const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, std::max(a.val(), b.val()));
      if (a.hasFastAccess() && b.hasFastAccess()) {
	if (a.val() >= b.val())
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = a.fastAccessDx(i);
	else
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = b.fastAccessDx(i);
      }
      else if (a.hasFastAccess()) {
	if (a.val() >= b.val())
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = a.fastAccessDx(i);
	else
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = 0.0;
      }
      else if (b.hasFastAccess()) {
	if (a.val() >= b.val())
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = 0.0;
	else
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = b.fastAccessDx(i);
      }

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    max(const typename SimpleFad<ValueT>::value_type& a, 
	const SimpleFad<ValueT>& b) {
      SimpleFad<ValueT> c(b.size(), std::max(a, b.val()));
      if (a >= b.val())
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = 0.0;
      else
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = b.fastAccessDx(i);

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    max(const SimpleFad<ValueT>& a, 
	const typename SimpleFad<ValueT>::value_type& b) {
      SimpleFad<ValueT> c(a.size(), std::max(a.val(), b));
      if (a.val() >= b)
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i);
      else
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = 0.0;

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    min(const SimpleFad<ValueT>& a, 
	const SimpleFad<ValueT>& b) {
      int sz = a.size() >= b.size() ? a.size() : b.size();
      SimpleFad<ValueT> c(sz, std::min(a.val(), b.val()));
      if (a.hasFastAccess() && b.hasFastAccess()) {
	if (a.val() <= b.val())
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = a.fastAccessDx(i);
	else
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = b.fastAccessDx(i);
      }
      else if (a.hasFastAccess()) {
	if (a.val() <= b.val())
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = a.fastAccessDx(i);
	else
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = 0.0;
      }
      else if (b.hasFastAccess()) {
	if (a.val() <= b.val())
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = 0.0;
	else
	  for (int i=0; i<sz; i++)
	    c.fastAccessDx(i) = b.fastAccessDx(i);
      }

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    min(const typename SimpleFad<ValueT>::value_type& a, 
	const SimpleFad<ValueT>& b) {
      SimpleFad<ValueT> c(b.size(), std::min(a, b.val()));
      if (a <= b.val())
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = 0.0;
      else
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = b.fastAccessDx(i);

      return c;
    }

    template <typename ValueT>
    SimpleFad<ValueT>
    min(const SimpleFad<ValueT>& a, 
	const typename SimpleFad<ValueT>::value_type& b) {
      SimpleFad<ValueT> c(a.size(), std::min(a.val(), b));
      if (a.val() <= b)
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = a.fastAccessDx(i);
      else
	for (int i=0; i<c.size(); i++)
	  c.fastAccessDx(i) = 0.0;

      return c;
    }

  } // namespace Fad

} // namespace Sacado

    //-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace Fad {							\
    template <typename ValueT>						\
    inline bool								\
    operator OP (const SimpleFad<ValueT>& a,				\
		 const SimpleFad<ValueT>& b)				\
    {									\
      return a.val() OP b.val();					\
    }									\
									\
    template <typename ValueT>						\
    inline bool								\
    operator OP (const ValueT& a,					\
		 const SimpleFad<ValueT>& b)				\
    {									\
      return a OP b.val();						\
    }									\
									\
    template <typename ValueT>						\
    inline bool								\
    operator OP (const SimpleFad<ValueT>& a,				\
		 const ValueT& b)					\
    {									\
      return a.val() OP b;						\
    }									\
  }									\
}

FAD_RELOP_MACRO(==)
FAD_RELOP_MACRO(!=)
FAD_RELOP_MACRO(<)
FAD_RELOP_MACRO(>)
FAD_RELOP_MACRO(<=)
FAD_RELOP_MACRO(>=)
FAD_RELOP_MACRO(<<=)
FAD_RELOP_MACRO(>>=)
FAD_RELOP_MACRO(&)
FAD_RELOP_MACRO(|)

#undef FAD_RELOP_MACRO

namespace Sacado {

  namespace Fad {

    template <typename ValueT>
    inline bool operator ! (const SimpleFad<ValueT>& a) 
    {
      return ! a.val();
    }

  } // namespace Fad

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace Fad {

    template <typename T>
    bool toBool(const SimpleFad<T>& x) {
      bool is_zero = (x.val() == 0.0);
      for (int i=0; i<x.size(); i++)
	is_zero = is_zero && (x.dx(i) == 0.0);
      return !is_zero;
    }

  } // namespace Fad

} // namespace Sacado

#define FAD_BOOL_MACRO(OP)						\
namespace Sacado {							\
  namespace Fad {							\
    template <typename T1, typename T2>					\
    inline bool								\
    operator OP (const SimpleFad<T1>& expr1,				\
		 const SimpleFad<T2>& expr2)				\
    {									\
      return toBool(expr1) OP toBool(expr2);				\
    }									\
									\
    template <typename T2>						\
    inline bool								\
    operator OP (const typename SimpleFad<T2>::value_type& a,		\
		 const SimpleFad<T2>& expr2)				\
    {									\
      return a OP toBool(expr2);					\
    }									\
									\
    template <typename T1>						\
    inline bool								\
    operator OP (const SimpleFad<T1>& expr1,				\
		 const typename SimpleFad<T1>::value_type& b)		\
    {									\
      return toBool(expr1) OP b;					\
    }									\
  }									\
}

FAD_BOOL_MACRO(&&)
FAD_BOOL_MACRO(||)

#undef FAD_BOOL_MACRO

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace Fad {

    template <typename ValueT>
    std::ostream& operator << (std::ostream& os, 
			       const SimpleFad<ValueT>& x) {
      os << x.val() << " [";
      
      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

 } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_SIMPLEFADOPS_HPP
