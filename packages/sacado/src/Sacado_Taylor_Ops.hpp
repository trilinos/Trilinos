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

#ifndef SACADO_TAYLOR_OPS_HPP
#define SACADO_TAYLOR_OPS_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_Taylor_Expression.hpp"

// Import the standard math functions into the Sacado::Taylor namespace
namespace Sacado {
  namespace Taylor {
    using std::exp;
    using std::log;
    using std::log10;
    using std::sqrt;
    using std::cos;
    using std::sin;
    using std::tan;
    using std::acos;
    using std::asin;
    using std::atan;
    using std::cosh;
    using std::sinh;
    using std::tanh;
    using std::abs;
    using std::fabs;
    using std::pow;
    using std::max;
    using std::min;
  }
}

namespace Sacado {

  namespace Taylor {

    // ---------------------- Unary Addition operator ------------------------

    template <typename ExprT>
    class UnaryPlusOp {
    public:

      typedef typename ExprT::value_type value_type;

      UnaryPlusOp(const ExprT& expr) {}

      void allocateCache(unsigned int d) const {}

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	return expr.coeff(i);
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const {
	return expr.fastAccessCoeff(i);
      }

    }; // class UnaryPlusOp

    // ---------------------- Unary Subtraction operator ---------------------

    template <typename ExprT>
    class UnaryMinusOp {
    public:

      typedef typename ExprT::value_type value_type;

      UnaryMinusOp(const ExprT& expr) {}

      void allocateCache(unsigned int d) const {}

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	return -expr.coeff(i);
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const {
	return -expr.fastAccessCoeff(i);
      }

    }; // class UnaryPlusOp

    // -------------------------- exp() function -----------------------------

    template <typename ExprT>
    class ExpOp {
    public:

      typedef typename ExprT::value_type value_type;

      ExpOp(const ExprT& expr) :
	c(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::exp(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*c[k-j]*expr.coeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::exp(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*c[k-j]*expr.fastAccessCoeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ExpOp

    // -------------------------- log() function -----------------------------

    template <typename ExprT>
    class LogOp {
    public:

      typedef typename ExprT::value_type value_type;

      LogOp(const ExprT& expr) :
	c(),
	dc(-1) 
      {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::log(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    c[k] = value_type(k)*expr.coeff(k);
	    for (unsigned int j=1; j<=k-1; j++)
	      c[k] -= value_type(j)*expr.coeff(k-j)*c[j];
	    c[k] /= (value_type(k)*expr.fastAccessCoeff(0));
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::log(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    c[k] = value_type(k)*expr.fastAccessCoeff(k);
	    for (unsigned int j=1; j<=k-1; j++)
	      c[k] -= value_type(j)*expr.fastAccessCoeff(k-j)*c[j];
	    c[k] /= (value_type(k)*expr.fastAccessCoeff(0));
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class LogOp

    // -------------------------- sqrt() function -----------------------------

    template <typename ExprT>
    class SqrtOp {
    public:

      typedef typename ExprT::value_type value_type;

      SqrtOp(const ExprT& expr) :
	c(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::sqrt(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  value_type tmp = value_type(2)*c[0];
	  for (unsigned int k=dc+1; k<=i; k++) {
	    c[k] = expr.coeff(k);
	    for (unsigned int j=1; j<=k-1; j++)
	      c[k] -= c[j]*c[k-j];
	    c[k] /= tmp;
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::sqrt(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  value_type tmp = value_type(2)*c[0];
	  for (unsigned int k=dc+1; k<=i; k++) {
	    c[k] = expr.fastAccessCoeff(k);
	    for (unsigned int j=1; j<=k-1; j++)
	      c[k] -= c[j]*c[k-j];
	    c[k] /= tmp;
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class SqrtOp

    // -------------------------- cos() function -----------------------------

    template <typename ExprT>
    class CosOp {
    public:

      typedef typename ExprT::value_type value_type;

      CosOp(const ExprT& expr) :
	c(),
	s(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
	s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cos(expr.fastAccessCoeff(0));
	    s[0] = std::sin(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] -= value_type(j)*expr.coeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.coeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cos(expr.fastAccessCoeff(0));
	    s[0] = std::sin(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] -= value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class CosOp

    // -------------------------- sin() function -----------------------------

    template <typename ExprT>
    class SinOp {
    public:

      typedef typename ExprT::value_type value_type;

      SinOp(const ExprT& expr) :
	c(),
	s(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
	s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cos(expr.fastAccessCoeff(0));
	    s[0] = std::sin(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] -= value_type(j)*expr.coeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.coeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return s[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cos(expr.fastAccessCoeff(0));
	    s[0] = std::sin(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] -= value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return s[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class SinOp

    // -------------------------- cosh() function -----------------------------

    template <typename ExprT>
    class CoshOp {
    public:

      typedef typename ExprT::value_type value_type;

      CoshOp(const ExprT& expr) :
	c(),
	s(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
	s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cosh(expr.fastAccessCoeff(0));
	    s[0] = std::sinh(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] += value_type(j)*expr.coeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.coeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cosh(expr.fastAccessCoeff(0));
	    s[0] = std::sinh(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] += value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class CoshOp

    // -------------------------- sinh() function -----------------------------

    template <typename ExprT>
    class SinhOp {
    public:

      typedef typename ExprT::value_type value_type;

      SinhOp(const ExprT& expr) :
	c(),
	s(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
	s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cosh(expr.fastAccessCoeff(0));
	    s[0] = std::sinh(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] += value_type(j)*expr.coeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.coeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return s[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::cosh(expr.fastAccessCoeff(0));
	    s[0] = std::sinh(expr.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++) {
	      c[k] += value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
	      s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
	    }
	    c[k] /= value_type(k);
	    s[k] /= value_type(k);
	  }
	  dc = i;
	}
	return s[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class SinhOp

    // -------------------------- fabs() function -----------------------------

    template <typename ExprT>
    class FAbsOp {
    public:

      typedef typename ExprT::value_type value_type;

      FAbsOp(const ExprT& expr) {}

      void allocateCache(unsigned int d) const {}

      value_type computeCoeff(unsigned int i, const ExprT& expr) const {
	if (expr.fastAccessCoeff(0) > 0)
	  return expr.coeff(i);
	else
	  return -expr.coeff(i);
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT& expr) const
      {
	if (expr.fastAccessCoeff(0) > 0)
	  return expr.fastAccessCoeff(i);
	else
	  return -expr.fastAccessCoeff(i);
      }

    }; // class FAbsOp

  } // namespace Taylor

} // namespace Sacado

#define TAYLOR_UNARYOP_MACRO(OPNAME,OP)					\
namespace Sacado {							\
  namespace Taylor {							\
    template <typename T>						\
    inline Expr< UnaryExpr< Expr<T>, OP > >				\
    OPNAME (const Expr<T>& expr)					\
    {									\
      typedef UnaryExpr< Expr<T>, OP > expr_t;				\
      									\
      return Expr<expr_t>(expr_t(expr));				\
    }									\
  }									\
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Taylor::OPNAME;                                         \
}

TAYLOR_UNARYOP_MACRO(operator+, UnaryPlusOp)
TAYLOR_UNARYOP_MACRO(operator-, UnaryMinusOp)
TAYLOR_UNARYOP_MACRO(exp, ExpOp)
TAYLOR_UNARYOP_MACRO(log, LogOp)
TAYLOR_UNARYOP_MACRO(sqrt, SqrtOp)
TAYLOR_UNARYOP_MACRO(cos, CosOp)
TAYLOR_UNARYOP_MACRO(sin, SinOp)
TAYLOR_UNARYOP_MACRO(cosh, CoshOp)
TAYLOR_UNARYOP_MACRO(sinh, SinhOp)
TAYLOR_UNARYOP_MACRO(abs, FAbsOp)
TAYLOR_UNARYOP_MACRO(fabs, FAbsOp)

#undef TAYLOR_UNARYOP_MACRO

namespace Sacado {

  namespace Taylor {

    // ---------------------- Addition operator -----------------------------

    template <typename ExprT1, typename ExprT2>
    class AdditionOp {
    public:
      
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;
    
      AdditionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}
    
      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	return expr1.coeff(i) + expr2.coeff(i);
      }
      
      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1, 
			     const ExprT2& expr2) const {
	return expr1.fastAccessCoeff(i) + expr2.fastAccessCoeff(i);
      }
      
    }; // class AdditionOp

    template <typename ExprT1>
    class AdditionOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      AdditionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.coeff(i) + expr2.coeff(i);
	else
	  return expr1.coeff(i);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.fastAccessCoeff(i) + expr2.fastAccessCoeff(i);
	else
	  return expr1.fastAccessCoeff(i);
      }

    }; // class AdditionOp

    template <typename ExprT2>
    class AdditionOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      AdditionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.coeff(i) + expr2.coeff(i);
	else
	  return expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.fastAccessCoeff(i) + expr2.fastAccessCoeff(i);
	else
	  return expr2.fastAccessCoeff(i);
      }

    }; // class AdditionOp

    // ---------------------- Subtraction operator ---------------------------

    template <typename ExprT1, typename ExprT2>
    class SubtractionOp {
    public:
      
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;
    
      SubtractionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}
    
      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	return expr1.coeff(i) - expr2.coeff(i);
      }
      
      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1, 
			     const ExprT2& expr2) const {
	return expr1.fastAccessCoeff(i) - expr2.fastAccessCoeff(i);
      }
      
    }; // class SubtractionOp

    template <typename ExprT1>
    class SubtractionOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      SubtractionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.coeff(i) - expr2.coeff(i);
	else
	  return expr1.coeff(i);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.fastAccessCoeff(i) - expr2.fastAccessCoeff(i);
	else
	  return expr1.fastAccessCoeff(i);
      }

    }; // class SubtractionOp

    template <typename ExprT2>
    class SubtractionOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      SubtractionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.coeff(i) - expr2.coeff(i);
	else
	  return -expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return expr1.fastAccessCoeff(i) - expr2.fastAccessCoeff(i);
	else
	  return -expr2.fastAccessCoeff(i);
      }

    }; // class SubtractionOp

    // ---------------------- Multiplication operator -------------------------

    template <typename ExprT1, typename ExprT2>
    class MultiplicationOp {
    public:
      
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;
    
      MultiplicationOp(const ExprT1& expr1, const ExprT2 expr2) :
	c(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }
    
      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=0; j<=k; j++)
	      c[k] += expr1.coeff(j)*expr2.coeff(k-j);
	  }
	  dc = i;
	}
	return c[i];
      }
      
      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1, 
			     const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=0; j<=k; j++)
	      c[k] += expr1.fastAccessCoeff(j)*expr2.fastAccessCoeff(k-j);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;
      
    }; // class MultiplicationOp

    template <typename ExprT1>
    class MultiplicationOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      MultiplicationOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	return expr1.coeff(i)*expr2.value();
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	return expr1.fastAccessCoeff(i)*expr2.value();
      }

    }; // class MultiplicationOp

    template <typename ExprT2>
    class MultiplicationOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      MultiplicationOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	return expr1.value()*expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	return expr1.value()*expr2.fastAccessCoeff(i);
      }

    }; // class MultiplicationOp

    // ---------------------- Division operator -------------------------

    template <typename ExprT1, typename ExprT2>
    class DivisionOp {
    public:
      
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;
    
      DivisionOp(const ExprT1& expr1, const ExprT2 expr2) :
	c(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }
    
      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  for (unsigned int k=dc+1; k<=i; k++) {
	    c[k] = expr1.coeff(k);
	    for (unsigned int j=1; j<=k; j++)
	      c[k] -= expr2.coeff(j)*c[k-j];
	    c[k] /= expr2.fastAccessCoeff(0);
	  }
	  dc = i;
	}
	return c[i];
      }
      
      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1, 
			     const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  for (unsigned int k=dc+1; k<=i; k++) {
	    c[k] = expr1.coeff(k);
	    for (unsigned int j=1; j<=k; j++)
	      c[k] -= expr2.fastAccessCoeff(j)*c[k-j];
	    c[k] /= expr2.fastAccessCoeff(0);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;
      
    }; // class DivisionOp

    template <typename ExprT1>
    class DivisionOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      DivisionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	return expr1.coeff(i)/expr2.value();
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	return expr1.fastAccessCoeff(i)/expr2.value();
      }

    }; // class DivisionOp

    template <typename ExprT2>
    class DivisionOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      DivisionOp(const ExprT1& expr1, const ExprT2 expr2) :
	c(),
	dc(-1) {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = expr1.fastAccessCoeff(0) / expr2.fastAccessCoeff(0);
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] -= expr2.coeff(j)*c[k-j];
	    c[k] /= expr2.fastAccessCoeff(0);
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = expr1.fastAccessCoeff(0) / expr2.fastAccessCoeff(0);
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    c[k] = expr1.coeff(k);
	    for (unsigned int j=1; j<=k; j++)
	      c[k] -= expr2.fastAccessCoeff(j)*c[k-j];
	    c[k] /= expr2.fastAccessCoeff(0);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class DivisionOp

    // ---------------------- Max operator -----------------------------

    template <typename ExprT1, typename ExprT2>
    class MaxOp {
    public:
      
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;
    
      MaxOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}
    
      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return max(expr1.coeff(0), expr2.coeff(0));
	else
	  return expr1.coeff(0) >= expr2.coeff(0) ? expr1.coeff(i) : 
	    expr2.coeff(i);
      }
      
      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1, 
			     const ExprT2& expr2) const {
	if (i == 0)
	  return max(expr1.fastAccessCoeff(0), expr2.fastAccessCoeff(0));
	else
	  return expr1.fastAccessCoeff(0) >= expr2.fastAccessCoeff(0) ? 
	    expr1.fastAccessoeff(i) : expr2.fastAccessCoeff(i);
      }
      
    }; // class MaxOp

    template <typename ExprT1>
    class MaxOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      MaxOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return max(expr1.coeff(0), expr2.value());
	else
	  return expr1.coeff(0) >= expr2.value() ? expr1.coeff(i) : 
	    value_type(0);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return max(expr1.fastAccessCoeff(0), expr2.value());
	else
	  return expr1.fastAccessCoeff(0) >= expr2.value() ? 
	    expr1.fastAccessCoeff(i) : value_type(0);
      }

    }; // class MaxOp

    template <typename ExprT2>
    class MaxOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      MaxOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return max(expr1.value(), expr2.coeff(0));
	else
	  return expr1.value() >= expr2.coeff(0) ? value_type(0) : 
	    expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return max(expr1.value(), expr2.fastAccessCoeff(0));
	else
	  return expr1.value() >= expr2.fastAccessCoeff(0) ? value_type(0) : 
	    expr2.fastAccessCoeff(i);
      }

    }; // class MaxOp

    // ---------------------- Min operator -----------------------------

    template <typename ExprT1, typename ExprT2>
    class MinOp {
    public:
      
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;
    
      MinOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}
    
      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return min(expr1.coeff(0), expr2.coeff(0));
	else
	  return expr1.coeff(0) <= expr2.coeff(0) ? expr1.coeff(i) : 
	    expr2.coeff(i);
      }
      
      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1, 
			     const ExprT2& expr2) const {
	if (i == 0)
	  return min(expr1.fastAccessCoeff(0), expr2.fastAccessCoeff(0));
	else
	  return expr1.fastAccessCoeff(0) <= expr2.fastAccessCoeff(0) ? 
	    expr1.fastAccessCoeff(i) : expr2.fastAccessCoeff(i);
      }
      
    }; // class MinOp

    template <typename ExprT1>
    class MinOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      MinOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return min(expr1.coeff(0), expr2.value());
	else
	  return expr1.coeff(0) <= expr2.value() ? expr1.coeff(i) : 
	    value_type(0);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return min(expr1.fastAccessCoeff(0), expr2.value());
	else
	  return expr1.fastAccessCoeff(0) <= expr2.value() ? 
	    expr1.fastAccessCoeff(i) : value_type(0);
      }

    }; // class MinOp

    template <typename ExprT2>
    class MinOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      MinOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(unsigned int d) const {}

      value_type
      computeCoeff(unsigned int i, const ExprT1& expr1, 
		   const ExprT2& expr2) const {
	if (i == 0)
	  return min(expr1.value(), expr2.coeff(0));
	else
	  return expr1.value() <= expr2.coeff(0) ? value_type(0) : 
	    expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(unsigned int i, const ExprT1& expr1,
			     const ExprT2& expr2) const {
	if (i == 0)
	  return min(expr1.value(), expr2.fastAccessCoeff(0));
	else
	  return expr1.value() <= expr2.fastAccessCoeff(0) ? value_type(0) : 
	    expr2.fastAccessCoeff(i);
      }

    }; // class MinOp

    // ------------------------ quadrature function ---------------------------

    template <typename ExprT1, typename ExprT2>
    class ASinQuadOp {
    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      ASinQuadOp(const ExprT1& expr1, const ExprT2& expr2) :
	c(),
	dc(-1)
      {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT1& expr1,
			      const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::asin(expr1.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*expr2.coeff(k-j)*expr1.coeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT1& expr1,
					const ExprT2& expr2) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::asin(expr1.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*expr2.fastAccessCoeff(k-j)*
		expr1.fastAccessCoeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ASinQuadOp

    template <typename ExprT1, typename ExprT2>
    class ACosQuadOp {
    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      ACosQuadOp(const ExprT1& expr1, const ExprT2& expr2) :
	c(),
	dc(-1)
      {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT1& expr1,
			      const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::acos(expr1.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*expr2.coeff(k-j)*expr1.coeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT1& expr1,
					const ExprT2& expr2) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::acos(expr1.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*expr2.fastAccessCoeff(k-j)*
		expr1.fastAccessCoeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ACosQuadOp

    template <typename ExprT1, typename ExprT2>
    class ATanQuadOp {
    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;

      ATanQuadOp(const ExprT1& expr1, const ExprT2& expr2) :
	c(),
	dc(-1)
      {}

      void allocateCache(unsigned int d) const { 
	c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(unsigned int i, const ExprT1& expr1,
			      const ExprT2& expr2) const {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::atan(expr1.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*expr2.coeff(k-j)*expr1.coeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

      value_type computeFastAccessCoeff(unsigned int i, 
					const ExprT1& expr1,
					const ExprT2& expr2) const
      {
	if (static_cast<int>(i) > dc) {
	  if (dc < 0) {
	    c[0] = std::atan(expr1.fastAccessCoeff(0));
	    dc = 0;
	  }
	  for (unsigned int k=dc+1; k<=i; k++) {
	    for (unsigned int j=1; j<=k; j++)
	      c[k] += value_type(j)*expr2.fastAccessCoeff(k-j)*
		expr1.fastAccessCoeff(j);
	    c[k] /= value_type(k);
	  }
	  dc = i;
	}
	return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ATanQuadOp

  } // namespace Taylor

} // namespace Sacado

#define TAYLOR_BINARYOP_MACRO(OPNAME,OP)				\
namespace Sacado {							\
  namespace Taylor {							\
    template <typename T1, typename T2>					\
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, OP > >			\
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)		\
    {									\
      typedef BinaryExpr< Expr<T1>, Expr<T2>, OP > expr_t;		\
    									\
      return Expr<expr_t>(expr_t(expr1, expr2));			\
    }									\
									\
    template <typename T>						\
    inline Expr< BinaryExpr< ConstExpr<typename Expr<T>::value_type>,	\
			     Expr<T>, OP > >				\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef BinaryExpr< ConstT, Expr<T>, OP > expr_t;			\
									\
      return Expr<expr_t>(expr_t(ConstT(c), expr));			\
    }									\
									\
    template <typename T>						\
    inline Expr< BinaryExpr< Expr<T>,					\
			     ConstExpr<typename Expr<T>::value_type>,	\
			     OP > >					\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef BinaryExpr< Expr<T>, ConstT, OP > expr_t;			\
									\
      return Expr<expr_t>(expr_t(expr, ConstT(c)));			\
    }									\
  }									\
}

TAYLOR_BINARYOP_MACRO(operator+, AdditionOp)
TAYLOR_BINARYOP_MACRO(operator-, SubtractionOp)
TAYLOR_BINARYOP_MACRO(operator*, MultiplicationOp)
TAYLOR_BINARYOP_MACRO(operator/, DivisionOp)

#undef TAYLOR_BINARYOP_MACRO

  // The general definition of max/min works for Taylor variables too, except
  // we need to add a case when the argument types are different.  This 
  // can't conflict with the general definition, so we need to use
  // Substitution Failure Is Not An Error
#include "Sacado_mpl_disable_if.hpp"
#include "Sacado_mpl_is_same.hpp"

#define TAYLOR_SFINAE_BINARYOP_MACRO(OPNAME,OP)				\
namespace Sacado {							\
  namespace Taylor {							\
    template <typename T1, typename T2>					\
    inline                                                              \
    typename                                                            \
    mpl::disable_if< mpl::is_same<T1,T2>,                               \
                     Expr<BinaryExpr<Expr<T1>, Expr<T2>, OP> > >::type  \
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)		\
    {									\
      typedef BinaryExpr< Expr<T1>, Expr<T2>, OP > expr_t;		\
    									\
      return Expr<expr_t>(expr_t(expr1, expr2));			\
    }									\
									\
    template <typename T>						\
    inline Expr< BinaryExpr< ConstExpr<typename Expr<T>::value_type>,	\
			     Expr<T>, OP > >				\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef BinaryExpr< ConstT, Expr<T>, OP > expr_t;			\
									\
      return Expr<expr_t>(expr_t(ConstT(c), expr));			\
    }									\
									\
    template <typename T>						\
    inline Expr< BinaryExpr< Expr<T>,					\
			     ConstExpr<typename Expr<T>::value_type>,	\
			     OP > >					\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef BinaryExpr< Expr<T>, ConstT, OP > expr_t;			\
									\
      return Expr<expr_t>(expr_t(expr, ConstT(c)));			\
    }									\
  }									\
}

TAYLOR_SFINAE_BINARYOP_MACRO(max, MaxOp)
TAYLOR_SFINAE_BINARYOP_MACRO(min, MinOp)

#undef TAYLOR_SFINAE_BINARYOP_MACRO

namespace std {
  using Sacado::Taylor::min;
  using Sacado::Taylor::max;
}

namespace Sacado {

  namespace Taylor {

    template <typename T1, typename T2>
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, ASinQuadOp > >
    asin_quad (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef BinaryExpr< Expr<T1>, Expr<T2>, ASinQuadOp > expr_t;
 
      return Expr<expr_t>(expr_t(expr1, expr2));
    }

    template <typename T1, typename T2>
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, ACosQuadOp > >
    acos_quad (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef BinaryExpr< Expr<T1>, Expr<T2>, ACosQuadOp > expr_t;
 
      return Expr<expr_t>(expr_t(expr1, expr2));
    }

    template <typename T1, typename T2>
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, ATanQuadOp > >
    atan_quad (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef BinaryExpr< Expr<T1>, Expr<T2>, ATanQuadOp > expr_t;
 
      return Expr<expr_t>(expr_t(expr1, expr2));
    }

    template <typename ExprT1, typename ExprT2>
    struct PowExprType {
      typedef UnaryExpr< ExprT1, LogOp > T3;
      typedef BinaryExpr< ExprT2, Expr<T3>, MultiplicationOp > T4;
      typedef UnaryExpr< Expr<T4>, ExpOp > T5;
      
      typedef Expr<T5> expr_type;
    };

     template <typename ExprT2>
     struct PowExprType< typename ExprT2::value_type, ExprT2 > {
       typedef typename ExprT2::value_type T1;
      typedef BinaryExpr< ExprT2, ConstExpr<T1>, MultiplicationOp > T4;
      typedef UnaryExpr< Expr<T4>, ExpOp > T5;
      
      typedef Expr<T5> expr_type;
    };

    template <typename ExprT1>
    struct PowExprType<ExprT1,typename ExprT1::value_type> {
      typedef typename ExprT1::value_type T2;
      typedef UnaryExpr< ExprT1, LogOp > T3;
      typedef BinaryExpr< ConstExpr<T2>, Expr<T3>, MultiplicationOp > T4;
      typedef UnaryExpr< Expr<T4>, ExpOp > T5;
      
      typedef Expr<T5> expr_type;
    };

    template <typename T1, typename T2>
    inline typename PowExprType< Expr<T1>, Expr<T2> >::expr_type
    pow (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      // pow(x,y) = exp(y*log(x))
      return exp(expr2*log(expr1));
    }

    template <typename T>
    inline typename PowExprType< typename Expr<T>::value_type, Expr<T> >::expr_type
    pow (const typename Expr<T>::value_type& c, const Expr<T>& expr)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;

      // pow(x,y) = exp(y*log(x))
      return exp(expr*std::log(c));
    }

    template <typename T>
    inline typename PowExprType< Expr<T>, typename Expr<T>::value_type >::expr_type
    pow (const Expr<T>& expr, const typename Expr<T>::value_type& c)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;

      // pow(x,y) = exp(y*log(x))
      return exp(c*log(expr));
    }

    template <typename T>
    struct Log10ExprType {
      typedef UnaryExpr< Expr<T>, LogOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< Expr<T1>, ConstT, DivisionOp > T2;
      typedef Expr<T2> expr_type;
    };

    template <typename T>
    inline typename Log10ExprType<T>::expr_type
    log10 (const Expr<T>& expr)
    {
      // log10(x) = log(x)/log(10)
      return log(expr)/std::log(typename Expr<T>::value_type(10));
    }

    template <typename T>
    struct TanExprType {
      typedef UnaryExpr< Expr<T>, SinOp > T1;
      typedef UnaryExpr< Expr<T>, CosOp > T2;
      typedef BinaryExpr< Expr<T1>, Expr<T2>, DivisionOp > T3;
      typedef Expr<T3> expr_type;
    };

    template <typename T>
    inline typename TanExprType<T>::expr_type
    tan (const Expr<T>& expr)
    {
      // tan(x) = sin(x)/cos(x)
      return sin(expr)/cos(expr);
    }

    template <typename T>
    struct ASinExprType {
      typedef BinaryExpr< Expr<T>, Expr<T>, MultiplicationOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< ConstT, Expr<T1>, SubtractionOp > T2;
      typedef UnaryExpr< Expr<T2>, SqrtOp > T3;
      typedef BinaryExpr< ConstT, Expr<T3>, DivisionOp > T4;
      typedef BinaryExpr< Expr<T>, Expr<T4>, ASinQuadOp > T5;
      typedef Expr<T5> expr_type;
    };

    template <typename T>
    inline typename ASinExprType<T>::expr_type
    asin (const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type value_type;

      // asin(x) = integral of 1/sqrt(1-x*x)
      return asin_quad(expr, value_type(1)/sqrt(value_type(1)-expr*expr));
    }

    template <typename T>
    struct ACosExprType {
      typedef BinaryExpr< Expr<T>, Expr<T>, MultiplicationOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< ConstT, Expr<T1>, SubtractionOp > T2;
      typedef UnaryExpr< Expr<T2>, SqrtOp > T3;
      typedef BinaryExpr< ConstT, Expr<T3>, DivisionOp > T4;
      typedef BinaryExpr< Expr<T>, Expr<T4>, ACosQuadOp > T5;
      typedef Expr<T5> expr_type;
    };

    template <typename T>
    inline typename ACosExprType<T>::expr_type
    acos (const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type value_type;

      // acos(x) = integral of -1/sqrt(1-x*x)
      return acos_quad(expr, value_type(-1)/sqrt(value_type(1)-expr*expr));
    }

    template <typename T>
    struct ATanExprType {
      typedef BinaryExpr< Expr<T>, Expr<T>, MultiplicationOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< ConstT, Expr<T1>, AdditionOp > T2;
      typedef BinaryExpr< ConstT, Expr<T2>, DivisionOp > T3;
      typedef BinaryExpr< Expr<T>, Expr<T3>, ATanQuadOp > T4;
      typedef Expr<T4> expr_type;
    };

    template <typename T>
    inline typename ATanExprType<T>::expr_type
    atan (const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type value_type;

      // atan(x) = integral of 1/(1+x*x)
      return atan_quad(expr, value_type(1)/(value_type(1)+expr*expr));
    }

    template <typename T>
    struct TanhExprType {
      typedef UnaryExpr< Expr<T>, SinhOp > T1;
      typedef UnaryExpr< Expr<T>, CoshOp > T2;
      typedef BinaryExpr< Expr<T1>, Expr<T2>, DivisionOp > T3;
      typedef Expr<T3> expr_type;
    };

    template <typename T>
    inline typename TanhExprType<T>::expr_type
    tanh (const Expr<T>& expr)
    {
      // tanh(x) = sinh(x)/cosh(x)
      return sinh(expr)/cosh(expr);
    }

  } // namespace Taylor

} // namespace Sacado

namespace std {
  using Sacado::Taylor::pow;
  using Sacado::Taylor::log10;
  using Sacado::Taylor::tan;
  using Sacado::Taylor::asin;
  using Sacado::Taylor::acos;
  using Sacado::Taylor::atan;
  using Sacado::Taylor::tanh;
}

//-------------------------- Relational Operators -----------------------

#define TAYLOR_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace Taylor {							\
    template <typename ExprT1, typename ExprT2>				\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return expr1.fastAccessCoeff(0) OP expr2.fastAccessCoeff(0);	\
    }									\
									\
    template <typename ExprT2>						\
    inline bool								\
    operator OP (const typename Expr<ExprT2>::value_type& a,		\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return a OP expr2.fastAccessCoeff(0);				\
    }									\
									\
    template <typename ExprT1>						\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const typename Expr<ExprT1>::value_type& b)		\
    {									\
      return expr1.fastAccessCoeff(0) OP b;				\
    }									\
  }									\
}

TAYLOR_RELOP_MACRO(==)
TAYLOR_RELOP_MACRO(!=)
TAYLOR_RELOP_MACRO(<)
TAYLOR_RELOP_MACRO(>)
TAYLOR_RELOP_MACRO(<=)
TAYLOR_RELOP_MACRO(>=)
TAYLOR_RELOP_MACRO(<<=)
TAYLOR_RELOP_MACRO(>>=)
TAYLOR_RELOP_MACRO(&)
TAYLOR_RELOP_MACRO(|)

#undef TAYLOR_RELOP_MACRO

namespace Sacado {

  namespace Taylor {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.fastAccessCoeff(0);
    }

  } // namespace Taylor

} // namespace Sacado

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace Taylor {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os.setf(std::ios::fixed, std::ios::floatfield);
      os.width(12);
      os << "[";
      
      for (unsigned int i=0; i<=x.degree(); i++) {
	os.width(12);
	os << x.coeff(i);
      }

      os << "]\n";
      return os;
    }

  } // namespace Taylor

} // namespace Sacado


#endif // SACADO_TAYLOR_OPS_HPP
