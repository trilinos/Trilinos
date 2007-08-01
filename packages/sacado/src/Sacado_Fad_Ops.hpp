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

#ifndef SACADO_FAD_OPS_HPP
#define SACADO_FAD_OPS_HPP

#include "Sacado_Fad_Expression.hpp"
#include <cmath>
#include <algorithm>	// for std::min and std::max
#include <ostream>	// for std::ostream

#define FAD_UNARYOP_MACRO(OPNAME,OP,VALUE,DX,FASTACCESSDX)		\
namespace Sacado {							\
  namespace Fad {							\
									\
    template <typename ExprT>						\
    class OP {								\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
									\
      OP() {}						                \
									\
      static value_type computeValue(const ExprT& expr) {		\
	return VALUE;							\
      }									\
									\
      static value_type computeDx(int i, const ExprT& expr) {		\
	return DX;							\
      }									\
									\
      static value_type computeFastAccessDx(int i, const ExprT& expr) {	\
	return FASTACCESSDX;						\
      }									\
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
}

FAD_UNARYOP_MACRO(operator+,
		  UnaryPlusOp, 
		  expr.val(),
		  expr.dx(i),
		  expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(operator-,
		  UnaryMinusOp, 
		  -expr.val(),
		  -expr.dx(i),
		  -expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(exp,
		  ExpOp, 
		  std::exp(expr.val()),
		  std::exp(expr.val())*expr.dx(i),
		  std::exp(expr.val())*expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  std::log(expr.val()),
		  expr.dx(i)/expr.val(),
		  expr.fastAccessDx(i)/expr.val())
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  std::log10(expr.val()),
		  expr.dx(i)/( std::log(value_type(10))*expr.val()),
		  expr.fastAccessDx(i) / ( std::log(value_type(10))*expr.val()))
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp, 
		  std::sqrt(expr.val()),
		  expr.dx(i)/(value_type(2)* std::sqrt(expr.val())),
		  expr.fastAccessDx(i)/(value_type(2)* std::sqrt(expr.val())))
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  std::cos(expr.val()),
		  -expr.dx(i)* std::sin(expr.val()),
		  -expr.fastAccessDx(i)* std::sin(expr.val()))
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  std::sin(expr.val()),
		  expr.dx(i)* std::cos(expr.val()),
		  expr.fastAccessDx(i)* std::cos(expr.val()))
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  std::tan(expr.val()),
		  expr.dx(i)*
		    (value_type(1)+ std::tan(expr.val())* std::tan(expr.val())),
		  expr.fastAccessDx(i)*
		    (value_type(1)+ std::tan(expr.val())* std::tan(expr.val())))
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  std::acos(expr.val()),
		  -expr.dx(i)/ std::sqrt(value_type(1)-expr.val()*expr.val()),
		  -expr.fastAccessDx(i) /
		    std::sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  std::asin(expr.val()),
		  expr.dx(i)/ std::sqrt(value_type(1)-expr.val()*expr.val()),
		  expr.fastAccessDx(i) /
		    std::sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  std::atan(expr.val()),
		  expr.dx(i)/(value_type(1)+expr.val()*expr.val()),
		  expr.fastAccessDx(i)/(value_type(1)+expr.val()*expr.val()))
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  std::cosh(expr.val()),
		  expr.dx(i)* std::sinh(expr.val()),
		  expr.fastAccessDx(i)* std::sinh(expr.val()))
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  std::sinh(expr.val()),
		  expr.dx(i)* std::cosh(expr.val()),
		  expr.fastAccessDx(i)* std::cosh(expr.val()))
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  std::tanh(expr.val()),
		  expr.dx(i)/( std::cosh(expr.val())* std::cosh(expr.val())),
		  expr.fastAccessDx(i) / 
		    ( std::cosh(expr.val())* std::cosh(expr.val())))
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  acosh(expr.val()),
		  expr.dx(i)/ std::sqrt((expr.val()-value_type(1)) * 
				       (expr.val()+value_type(1))),
		  expr.fastAccessDx(i)/ std::sqrt((expr.val()-value_type(1)) * 
						 (expr.val()+value_type(1))))
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  asinh(expr.val()),
		  expr.dx(i)/ std::sqrt(value_type(1)+expr.val()*expr.val()),
		  expr.fastAccessDx(i)/ std::sqrt(value_type(1)+
						 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  atanh(expr.val()),
		  expr.dx(i)/(value_type(1)-expr.val()*expr.val()),
		  expr.fastAccessDx(i)/(value_type(1)-
						 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(abs,
		  AbsOp, 
		  std::abs(expr.val()),
		  expr.val() >= 0 ? value_type(+expr.dx(i)) : 
		    value_type(-expr.dx(i)),
		  expr.val() >= 0 ? value_type(+expr.fastAccessDx(i)) : 
		    value_type(-expr.fastAccessDx(i)))
FAD_UNARYOP_MACRO(fabs,
		  FAbsOp, 
		  std::fabs(expr.val()),
		  expr.val() >= 0 ? value_type(+expr.dx(i)) : 
		    value_type(-expr.dx(i)),
		  expr.val() >= 0 ? value_type(+expr.fastAccessDx(i)) : 
		    value_type(-expr.fastAccessDx(i)))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,VALUE,DX,FASTACCESSDX,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {							\
  namespace Fad {							\
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
      OP() {}			                                        \
									\
      static value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2) {	        \
	return VALUE;							\
      }									\
									\
      static value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) {				        \
	return DX;							\
      }									\
									\
      static value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) {			\
	return FASTACCESSDX;						\
      }									\
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
      OP() {}			                                        \
									\
      static value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2) {	        \
	return VALUE;							\
      }									\
									\
      static value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) {				        \
	return CONST_DX_2;						\
      }									\
									\
      static value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) {			\
	return CONST_FASTACCESSDX_2;					\
      }									\
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
      OP() {}			                                        \
									\
      static value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2){	        \
	return VALUE;							\
      }									\
									\
      static value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) {				        \
	return CONST_DX_1;						\
      }									\
									\
      static value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) {			\
	return CONST_FASTACCESSDX_1;					\
      }									\
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
}


FAD_BINARYOP_MACRO(operator+,
		   AdditionOp, 
		   expr1.val() + expr2.val(),
		   expr1.dx(i) + expr2.dx(i),
		   expr1.fastAccessDx(i) + expr2.fastAccessDx(i),
		   expr2.dx(i),
		   expr1.dx(i),
		   expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator-,
		   SubtractionOp, 
		   expr1.val() - expr2.val(),
		   expr1.dx(i) - expr2.dx(i),
		   expr1.fastAccessDx(i) - expr2.fastAccessDx(i),
		   -expr2.dx(i),
		   expr1.dx(i),
		   -expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator*,
		   MultiplicationOp, 
		   expr1.val() * expr2.val(),
		   expr1.val()*expr2.dx(i) + expr1.dx(i)*expr2.val(),
		   expr1.val()*expr2.fastAccessDx(i) + 
		     expr1.fastAccessDx(i)*expr2.val(),
		   expr1.val()*expr2.dx(i),
		   expr1.dx(i)*expr2.val(),
		   expr1.val()*expr2.fastAccessDx(i),
		   expr1.fastAccessDx(i)*expr2.val())
FAD_BINARYOP_MACRO(operator/,
		   DivisionOp, 
		   expr1.val() / expr2.val(),
		   (expr1.dx(i)*expr2.val() - expr2.dx(i)*expr1.val()) /
		     (expr2.val()*expr2.val()),
		   (expr1.fastAccessDx(i)*expr2.val() - 
		      expr2.fastAccessDx(i)*expr1.val()) /
		      (expr2.val()*expr2.val()),
		   -expr2.dx(i)*expr1.val() / (expr2.val()*expr2.val()),
		   expr1.dx(i)/expr2.val(),
		   -expr2.fastAccessDx(i)*expr1.val() / 
		     (expr2.val()*expr2.val()),
		   expr1.fastAccessDx(i)/expr2.val())
FAD_BINARYOP_MACRO(atan2,
		   Atan2Op,
		   std::atan2(expr1.val(), expr2.val()),
		   (expr2.val()*expr1.dx(i) - expr1.val()*expr2.dx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   (expr2.val()*expr1.fastAccessDx(i) - expr1.val()*expr2.fastAccessDx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   (-expr1.val()*expr2.dx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   (expr2.val()*expr1.dx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   (-expr1.val()*expr2.fastAccessDx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   (expr2.val()*expr1.fastAccessDx(i))/
			(expr1.val()*expr1.val() + expr2.val()*expr2.val()))
FAD_BINARYOP_MACRO(pow,
		   PowerOp,
		   std::pow(expr1.val(), expr2.val()),
		   (expr2.dx(i)* std::log(expr1.val())+expr2.val()*expr1.dx(i)/
		    expr1.val())* std::pow(expr1.val(),expr2.val()),
		   (expr2.fastAccessDx(i)* std::log(expr1.val())+
		    expr2.val()*expr1.fastAccessDx(i)/
		    expr1.val())* std::pow(expr1.val(),expr2.val()),
		   expr2.dx(i)* std::log(expr1.val()) *
		     std::pow(expr1.val(),expr2.val()),
		   expr2.val()*expr1.dx(i)/
		   expr1.val()* std::pow(expr1.val(),expr2.val()),
		   expr2.fastAccessDx(i)* std::log(expr1.val()) *
		     std::pow(expr1.val(),expr2.val()),
		   expr2.val()*expr1.fastAccessDx(i)/
		     expr1.val()* std::pow(expr1.val(),expr2.val()))
FAD_BINARYOP_MACRO(max,
                   MaxOp,
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

#define FAD_SFINAE_BINARYOP_MACRO(OPNAME,OP,VALUE,DX,FASTACCESSDX,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {							\
  namespace Fad {							\
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
      OP() {}			                                        \
									\
      static value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2) {	        \
	return VALUE;							\
      }									\
									\
      static value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) {				        \
	return DX;							\
      }									\
									\
      static value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) {			\
	return FASTACCESSDX;						\
      }									\
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
      OP() {}			                                        \
									\
      static value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2) {	        \
	return VALUE;							\
      }									\
									\
      static value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) {				        \
	return CONST_DX_2;						\
      }									\
									\
      static value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) {			\
	return CONST_FASTACCESSDX_2;					\
      }									\
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
      OP() {}			                                        \
									\
      static value_type							\
      computeValue(const ExprT1& expr1, const ExprT2& expr2){	        \
	return VALUE;							\
      }									\
									\
      static value_type							\
      computeDx(int i, const ExprT1& expr1,				\
		const ExprT2& expr2) {				        \
	return CONST_DX_1;						\
      }									\
									\
      static value_type							\
      computeFastAccessDx(int i, const ExprT1& expr1,			\
			  const ExprT2& expr2) {			\
	return CONST_FASTACCESSDX_1;					\
      }									\
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
//                    MaxOp,
//                    std::max(expr1.val(), expr2.val()),
//                    expr1.val() >= expr2.val() ? expr1.dx(i) : expr2.dx(i),
//                    expr1.val() >= expr2.val() ? expr1.fastAccessDx(i) : 
//                                                 expr2.fastAccessDx(i),
//                    expr1.val() >= expr2.val() ? value_type(0) : expr2.dx(i),
//                    expr1.val() >= expr2.val() ? expr1.dx(i) : value_type(0),
//                    expr1.val() >= expr2.val() ? value_type(0) : 
//                                                 expr2.fastAccessDx(i),
//                    expr1.val() >= expr2.val() ? expr1.fastAccessDx(i) : 
//                                                 value_type(0))
// FAD_SFINAE_BINARYOP_MACRO(min,
//                    MinOp,
//                    std::min(expr1.val(), expr2.val()),
//                    expr1.val() <= expr2.val() ? expr1.dx(i) : expr2.dx(i),
//                    expr1.val() <= expr2.val() ? expr1.fastAccessDx(i) : 
//                                                 expr2.fastAccessDx(i),
//                    expr1.val() <= expr2.val() ? value_type(0) : expr2.dx(i),
//                    expr1.val() <= expr2.val() ? expr1.dx(i) : value_type(0),
//                    expr1.val() <= expr2.val() ? value_type(0) : 
//                                                 expr2.fastAccessDx(i),
//                    expr1.val() <= expr2.val() ? expr1.fastAccessDx(i) : 
//                                                 value_type(0))

#undef FAD_SFINAE_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace Fad {							\
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

  namespace Fad {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace Fad

} // namespace Sacado

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace Fad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os << x.val() << "  [";
      
      for (int i=0; i< x.size(); i++) {
	os << " " << x.dx(i);
      }

      os << " ]\n";
      return os;
    }

  } // namespace Fad

} // namespace Sacado


#endif // SACADO_FAD_OPS_HPP
