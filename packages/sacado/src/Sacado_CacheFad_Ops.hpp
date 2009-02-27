// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.  
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  A short implementation ( not all operators and 
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_OPS_HPP
#define SACADO_CACHEFAD_OPS_HPP

#include "Sacado_CacheFad_Expression.hpp"
#include "Sacado_cmath.hpp"
#include <ostream>	// for std::ostream

#define FAD_UNARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE,DX,FASTACCESSDX)	\
namespace Sacado {							\
  namespace CacheFad {							\
									\
    template <typename ExprT>						\
    class OP {};							\
									\
    template <typename ExprT>						\
    class Expr< OP<ExprT> > {						\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
									\
      Expr(const ExprT& expr_) : expr(expr_)  {}			\
									\
      int size() const { return expr.size(); }				\
									\
      bool hasFastAccess() const { return expr.hasFastAccess(); }	\
									\
      bool isPassive() const { return expr.isPassive();}		\
									\
      value_type val() const {						\
	v = expr.val();							\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type dx(int i) const {					\
	return DX;							\
      }									\
									\
      value_type fastAccessDx(int i) const {				\
	return FASTACCESSDX;						\
      }									\
									\
    protected:								\
									\
      const ExprT& expr;						\
      mutable value_type v;						\
      mutable value_type a;						\
    };									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T> > >					\
    OPNAME (const Expr<T>& expr)					\
    {									\
      typedef OP< Expr<T> > expr_t;					\
      									\
      return Expr<expr_t>(expr);					\
    }									\
  }									\
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
		  a = std::exp(v),
		  a,
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  ;,
		  std::log(v),
		  expr.dx(i)/v,
		  expr.fastAccessDx(i)/v)
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  a = std::log(value_type(10))*v,
		  std::log10(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp,
		  a = value_type(2)*std::sqrt(v),
		  std::sqrt(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  a = std::sin(v),
		  std::cos(v),
		  -expr.dx(i)*a,
		  -expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  a = std::cos(v),
		  std::sin(v),
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  value_type t = std::tan(v); a = value_type(1)+t*t,
		  t,
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  a = - std::sqrt(value_type(1)-v*v),
		  std::acos(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  a = std::sqrt(value_type(1)-v*v),
		  std::asin(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  a = (value_type(1)+v*v),
		  std::atan(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  a = std::sinh(v),
		  std::cosh(v),
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  a = std::cosh(v),
		  std::sinh(v),
		  expr.dx(i)*a,
		  expr.fastAccessDx(i)*a)
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  a = std::cosh(v); a = a*a,
		  std::tanh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  a = std::sqrt((v-value_type(1))*(v+value_type(1))),
		  std::acosh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  a = std::sqrt(value_type(1)+v*v),
		  std::asinh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  a = value_type(1)-v*v,
		  std::atanh(v),
		  expr.dx(i)/a,
		  expr.fastAccessDx(i)/a)
FAD_UNARYOP_MACRO(abs,
		  AbsOp, 
		  ;,
		  std::abs(v),
		  v >= 0 ? value_type(+expr.dx(i)) : value_type(-expr.dx(i)),
		  v >= 0 ? value_type(+expr.fastAccessDx(i)) : 
		    value_type(-expr.fastAccessDx(i)))
FAD_UNARYOP_MACRO(fabs,
		  FAbsOp, 
		  ;,
		  std::fabs(v),
		  v >= 0 ? value_type(+expr.dx(i)) : value_type(-expr.dx(i)),
		  v >= 0 ? value_type(+expr.fastAccessDx(i)) : 
		    value_type(-expr.fastAccessDx(i)))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE,DX,FASTACCESSDX,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {							\
  namespace CacheFad {							\
									\
    template <typename ExprT1, typename ExprT2>				\
    class OP {};							\
									\
    template <typename ExprT1, typename ExprT2>				\
    class Expr< OP<ExprT1,ExprT2> > {					\
									\
    public:								\
									\
      typedef typename ExprT1::value_type value_type_1;			\
      typedef typename ExprT2::value_type value_type_2;			\
      typedef typename Sacado::Promote<value_type_1,			\
				       value_type_2>::type value_type;  \
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      int size() const {						\
	int sz1 = expr1.size(), sz2 = expr2.size();			\
	return sz1 > sz2 ? sz1 : sz2;					\
      }									\
									\
      bool hasFastAccess() const {					\
	return expr1.hasFastAccess() && expr2.hasFastAccess();		\
      }									\
									\
      bool isPassive() const {						\
	return expr1.isPassive() && expr2.isPassive();			\
      }									\
									\
      value_type val() const {						\
	v1 = expr1.val();						\
	v2 = expr2.val();						\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type dx(int i) const {					\
	return DX;							\
      }									\
									\
      value_type fastAccessDx(int i) const {				\
	return FASTACCESSDX;						\
      }									\
									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ExprT2& expr2;						\
      mutable value_type_1 v1;						\
      mutable value_type_2 v2;						\
      mutable value_type a;						\
      mutable value_type b;						\
    };									\
									\
    template <typename ExprT1>						\
    class Expr< OP<ExprT1, ConstExpr<typename ExprT1::value_type> >  >{ \
									\
    public:								\
									\
      typedef typename ExprT1::value_type value_type;			\
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;		\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      int size() const {						\
	return expr1.size();						\
      }									\
									\
      bool hasFastAccess() const {					\
	return expr1.hasFastAccess();					\
      }									\
									\
      bool isPassive() const {						\
	return expr1.isPassive();					\
      }									\
									\
      value_type val() const {						\
	v1 = expr1.val();						\
	v2 = expr2.val();						\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type dx(int i) const {					\
	return CONST_DX_2;						\
      }									\
									\
      value_type fastAccessDx(int i) const {				\
	return CONST_FASTACCESSDX_2;					\
      }									\
									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ExprT2 expr2;						\
      mutable value_type v1;						\
      mutable value_type v2;						\
      mutable value_type a;						\
      mutable value_type b;						\
    };									\
									\
    template <typename ExprT2>						\
    class Expr< OP<ConstExpr<typename ExprT2::value_type>, ExprT2 >  >{	\
									\
    public:								\
									\
      typedef typename ExprT2::value_type value_type;			\
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;		\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      int size() const {						\
	return expr2.size();						\
      }									\
									\
      bool hasFastAccess() const {					\
	return expr2.hasFastAccess();					\
      }									\
									\
      bool isPassive() const {						\
	return expr2.isPassive();					\
      }									\
									\
      value_type val() const {						\
	v1 = expr1.val();						\
	v2 = expr2.val();						\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      value_type dx(int i) const {					\
	return CONST_DX_1;						\
      }									\
									\
      value_type fastAccessDx(int i) const {				\
	return CONST_FASTACCESSDX_1;					\
      }									\
									\
    protected:								\
									\
      const ExprT1 expr1;						\
      const ExprT2& expr2;						\
      mutable value_type v1;						\
      mutable value_type v2;						\
      mutable value_type a;						\
      mutable value_type b;						\
    };									\
									\
    template <typename T1, typename T2>					\
    inline Expr< OP< Expr<T1>, Expr<T2> > >				\
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)		\
    {									\
      typedef OP< Expr<T1>, Expr<T2> > expr_t;				\
    									\
      return Expr<expr_t>(expr1, expr2);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< ConstExpr<typename Expr<T>::value_type>,		\
			     Expr<T> > >				\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef OP< ConstT, Expr<T> > expr_t;				\
									\
      return Expr<expr_t>(ConstT(c), expr);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>,						\
		     ConstExpr<typename Expr<T>::value_type> > >	\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;		\
      typedef OP< Expr<T>, ConstT > expr_t;				\
									\
      return Expr<expr_t>(expr, ConstT(c));				\
    }									\
  }									\
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
FAD_BINARYOP_MACRO(atan2,
		   Atan2Op,
		   a=v1*v1 + v2*v2,
		   std::atan2(v1,v2),
		   (expr1.dx(i)*v2 - expr2.dx(i)*v1)/a,
		   (expr1.fastAccessDx(i)*v2 - expr2.fastAccessDx(i)*v1)/a,
		   (-expr2.dx(i)*v1)/a,
		   (expr1.dx(i)*v2)/a,
		   (-expr2.fastAccessDx(i)*v1)/a,
		   (expr1.fastAccessDx(i)*v2)/a)
FAD_BINARYOP_MACRO(pow,
		   PowerOp,
		   value_type c= std::pow(v1,v2); a=c*v2/v1; b=c*std::log(v1),
		   c,
		   expr1.dx(i)*a + expr2.dx(i)*b,
		   expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b,
		   expr2.dx(i)*b,
		   expr1.dx(i)*a,
		   expr2.fastAccessDx(i)*b,
		   expr1.fastAccessDx(i)*a)
FAD_BINARYOP_MACRO(max,
		   MaxOp,
		   ;,
		   std::max(expr1.val(), expr2.val()),
		   expr1.val() >= expr2.val() ? expr1.dx(i) : expr2.dx(i),
		   expr1.val() >= expr2.val() ? expr1.fastAccessDx(i) : 
		                                expr2.fastAccessDx(i),
		   expr1.val() >= expr2.val() ? value_type(0) : expr2.dx(i),
		   expr1.val() >= expr2.val() ? expr1.dx(i) : value_type(0),
		   expr1.val() >= expr2.val() ? value_type(0) : 
		                                expr2.fastAccessDx(i),
		   expr1.val() >= expr2.val() ? expr1.fastAccessDx(i) : 
		                                value_type(0))
FAD_BINARYOP_MACRO(min,
		   MinOp,
		   ;,
		   std::min(expr1.val(), expr2.val()),
		   expr1.val() <= expr2.val() ? expr1.dx(i) : expr2.dx(i),
		   expr1.val() <= expr2.val() ? expr1.fastAccessDx(i) : 
		                                expr2.fastAccessDx(i),
		   expr1.val() <= expr2.val() ? value_type(0) : expr2.dx(i),
		   expr1.val() <= expr2.val() ? expr1.dx(i) : value_type(0),
		   expr1.val() <= expr2.val() ? value_type(0) : 
		                                expr2.fastAccessDx(i),
		   expr1.val() <= expr2.val() ? expr1.fastAccessDx(i) : 
		                                value_type(0))

#undef FAD_BINARYOP_MACRO
  
  // The general definition of max/min works for Fad variables too, except
  // we need to add a case when the argument types are different.  This 
  // can't conflict with the general definition, so we need to use
  // Substitution Failure Is Not An Error
#include "Sacado_mpl_disable_if.hpp"
#include "Sacado_mpl_is_same.hpp"

#define FAD_SFINAE_BINARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE,DX,FASTACCESSDX,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
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
				       value_type_2>::type value_type;  \
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

// FAD_SFINAE_BINARYOP_MACRO(max,
// 		   MaxOp,
// 		   ;,
// 		   std::max(expr1.val(), expr2.val()),
// 		   expr1.val() >= expr2.val() ? expr1.dx(i) : expr2.dx(i),
// 		   expr1.val() >= expr2.val() ? expr1.fastAccessDx(i) : 
// 		                                expr2.fastAccessDx(i),
// 		   expr1.val() >= expr2.val() ? value_type(0) : expr2.dx(i),
// 		   expr1.val() >= expr2.val() ? expr1.dx(i) : value_type(0),
// 		   expr1.val() >= expr2.val() ? value_type(0) : 
// 		                                expr2.fastAccessDx(i),
// 		   expr1.val() >= expr2.val() ? expr1.fastAccessDx(i) : 
// 		                                value_type(0))
// FAD_SFINAE_BINARYOP_MACRO(min,
// 		   MinOp,
// 		   ;,
// 		   std::min(expr1.val(), expr2.val()),
// 		   expr1.val() <= expr2.val() ? expr1.dx(i) : expr2.dx(i),
// 		   expr1.val() <= expr2.val() ? expr1.fastAccessDx(i) : 
// 		                                expr2.fastAccessDx(i),
// 		   expr1.val() <= expr2.val() ? value_type(0) : expr2.dx(i),
// 		   expr1.val() <= expr2.val() ? expr1.dx(i) : value_type(0),
// 		   expr1.val() <= expr2.val() ? value_type(0) : 
// 		                                expr2.fastAccessDx(i),
// 		   expr1.val() <= expr2.val() ? expr1.fastAccessDx(i) : 
// 		                                value_type(0))

#undef FAD_SFINAE_BINARYOP_MACRO

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

namespace Sacado {

  namespace CacheFad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os << x.val() << " [";
      
      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

  } // namespace CacheFad

} // namespace Sacado


#endif // SACADO_CACHEFAD_OPS_HPP
