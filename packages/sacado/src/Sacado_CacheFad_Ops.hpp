// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef SACADO_CACHEFAD_OPS_HPP
#define SACADO_CACHEFAD_OPS_HPP

#include "Sacado_CacheFad_Expression.hpp"

// Import the standard math functions into the Sacado::CacheFad namespace
namespace Sacado {
  namespace CacheFad {
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
  }
}

#define FAD_UNARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE,DX,FASTACCESSDX)	\
namespace Sacado {							\
  namespace CacheFad {							\
									\
    template <typename ExprT>						\
    class OP {								\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
									\
      OP(const ExprT& expr) {}						\
									\
      value_type computeValue(const ExprT& expr) const {		\
	v = expr.val();							\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type computeDx(int i, const ExprT& expr) const {		\
	return DX;							\
      }									\
									\
      value_type computeFastAccessDx(int i, const ExprT& expr) const {	\
	return FASTACCESSDX;						\
      }									\
									\
    protected:								\
									\
      mutable value_type v;						\
      mutable value_type a;						\
    };									\
									\
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
  using Sacado::CacheFad::OPNAME;                                       \
}

FAD_UNARYOP_MACRO(operator+,
		  UnaryPlusOp, 
		  ;,
		  v,
		  expr.dx(i),
		  expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(operator-,
		  UnaryMinusOp,
		  ;, 
		  -v,
		  -expr.dx(i),
		  -expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(exp,
		  ExpOp, 
		  a = exp(v),
		  a,
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  ;,
		  log(v),
		  expr.dx(i)/v,
		  expr.fastAccessDx(i)/v)
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  a = log(value_type(10))*v,
		  log10(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp,
		  a = value_type(2)*sqrt(v),
		  sqrt(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  a = sin(v),
		  cos(v),
		  -expr.dx(i)*a,
		  -expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  a = cos(v),
		  sin(v),
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  value_type t = tan(v); a = value_type(1)+t*t,
		  t,
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  a = -sqrt(value_type(1)-v*v),
		  acos(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  a = sqrt(value_type(1)-v*v),
		  asin(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  a = (value_type(1)+v*v),
		  atan(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  a = sinh(v),
		  cosh(v),
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  a = cosh(v),
		  sinh(v),
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  a = cosh(v); a = a*a,
		  tanh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  a = sqrt((v-value_type(1))/(v+value_type(1))),
		  acosh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  a = sqrt(value_type(1)+v*v),
		  asinh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  a = sqrt(value_type(1)-v*v),
		  atanh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(abs,
		  AbsOp, 
		  ;,
		  abs(v),
		  v >= 0 ? expr.dx(i) : -expr.dx(i),
		  v >= 0 ? expr.fastAccessDx(i) : -expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(fabs,
		  FAbsOp, 
		  ;,
		  fabs(v),
		  v >= 0 ? expr.dx(i) : -expr.dx(i),
		  v >= 0 ? expr.fastAccessDx(i) : -expr.fastAccessDx(i))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE,DX,FASTACCESSDX,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {							\
  namespace CacheFad {							\
									\
    template <typename ExprT1, typename ExprT2>				\
    class OP {								\
									\
    public:								\
									\
      typedef typename ExprT1::value_type value_type_1;			\
      typedef typename ExprT2::value_type value_type_2;			\
      typedef typename Sacado::Promote<value_type_1,			\
				       value_type_1>::type value_type;  \
									\
      OP(const ExprT1& expr1, const ExprT2& expr2) {}			\
									\
      value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2) const {	\
	v1 = expr1.val();						\
	v2 = expr2.val();						\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) const {				\
	return DX;							\
      }									\
									\
      value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) const {			\
	return FASTACCESSDX;						\
      }									\
									\
    protected:								\
									\
      mutable value_type_1 v1;						\
      mutable value_type_2 v2;						\
      mutable value_type a;						\
      mutable value_type b;						\
    };									\
									\
    template <typename ExprT1>						\
    class OP<ExprT1, ConstExpr<typename ExprT1::value_type> > {		\
									\
    public:								\
									\
      typedef typename ExprT1::value_type value_type;			\
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;		\
									\
      OP(const ExprT1& expr1, const ExprT2& expr2) {}			\
									\
      value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2) const {	\
	v1 = expr1.val();						\
	v2 = expr2.val();						\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) const {				\
	return CONST_DX_2;						\
      }									\
									\
      value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) const {			\
	return CONST_FASTACCESSDX_2;					\
      }									\
									\
    protected:								\
									\
      mutable value_type v1;						\
      mutable value_type v2;						\
      mutable value_type a;						\
      mutable value_type b;						\
    };									\
									\
    template <typename ExprT2>						\
    class OP<ConstExpr<typename ExprT2::value_type>, ExprT2 > {		\
									\
    public:								\
									\
      typedef typename ExprT2::value_type value_type;			\
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;		\
									\
      OP(const ExprT1& expr1, const ExprT2& expr2) {}			\
									\
      value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2) const {	\
	v1 = expr1.val();						\
	v2 = expr2.val();						\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) const {				\
	return CONST_DX_1;						\
      }									\
									\
      value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) const {			\
	return CONST_FASTACCESSDX_1;					\
      }									\
									\
    protected:								\
									\
      mutable value_type v1;						\
      mutable value_type v2;						\
      mutable value_type a;						\
      mutable value_type b;						\
    };									\
									\
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
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::CacheFad::OPNAME;                                       \
}

FAD_BINARYOP_MACRO(operator+,
		   AdditionOp, 
		   ;,
		   v1 + v2,
		   expr1.dx(i) + expr2.dx(i),
		   expr1.fastAccessDx(i) + expr2.fastAccessDx(i),
		   expr2.dx(i),
		   expr1.dx(i),
		   expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator-,
		   SubtractionOp, 
		   ;,
		   v1 - v2,
		   expr1.dx(i) - expr2.dx(i),
		   expr1.fastAccessDx(i) - expr2.fastAccessDx(i),
		   -expr2.dx(i),
		   expr1.dx(i),
		   -expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator*,
		   MultiplicationOp, 
		   ;,
		   v1*v2,
		   v1*expr2.dx(i) + expr1.dx(i)*v2,
		   v1*expr2.fastAccessDx(i) + expr1.fastAccessDx(i)*v2,
		   v1*expr2.dx(i),
		   expr1.dx(i)*v2,
		   v1*expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i)*v2)
FAD_BINARYOP_MACRO(operator/,
		   DivisionOp, 
		   value_type c = v1/v2; a = value_type(1)/v2; b = -c/v2,
		   c,
		   expr1.dx(i)*a + expr2.dx(i)*b,
		   expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b,
		   expr2.dx(i)*b,
		   expr1.dx(i)*a,
		   expr2.fastAccessDx(i)*b,
		   expr1.fastAccessDx(i)*a)
FAD_BINARYOP_MACRO(pow,
		   PowerOp,
		   value_type c=pow(v1,v2); a=c*v2/v1; b=c*log(v1),
		   c,
		   expr1.dx(i)*a + expr2.dx(i)*b,
		   expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b,
		   expr2.dx(i)*b,
		   expr1.dx(i)*a,
		   expr2.fastAccessDx(i)*b,
		   expr1.fastAccessDx(i)*a)

#undef FAD_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace CacheFad {							\
    template <typename ExprT1, typename ExprT2>				\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return expr1.val() OP expr2.val();				\
    }									\
									\
    template <typename ExprT2>						\
    inline bool								\
    operator OP (const typename Expr<ExprT2>::value_type& a,		\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return a OP expr2.val();						\
    }									\
									\
    template <typename ExprT1>						\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const typename Expr<ExprT1>::value_type& b)		\
    {									\
      return expr1.val() OP b;						\
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

  namespace CacheFad {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace CacheFad

} // namespace Sacado

//-------------------------- I/O Operators -----------------------

#include <iostream>

namespace Sacado {

  namespace CacheFad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os << x.val() << "  [";
      
      for (int i=0; i< x.size(); i++) {
	os << " " << x.dx(i);
      }

      os << " ]\n";
      return os;
    }

  } // namespace CacheFad

} // namespace Sacado


#endif // SACADO_CACHEFAD_OPS_HPP
