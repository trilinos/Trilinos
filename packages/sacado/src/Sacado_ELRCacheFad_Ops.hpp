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

#ifndef SACADO_ELRCACHEFAD_OPS_HPP
#define SACADO_ELRCACHEFAD_OPS_HPP

#include "Sacado_ELRCacheFad_Expression.hpp"
#include "Sacado_cmath.hpp"
#include <ostream>	// for std::ostream

#define FAD_UNARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE,ADJOINT)		\
namespace Sacado {							\
  namespace ELRCacheFad {						\
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
      typedef typename ExprT::base_expr_type base_expr_type;		\
									\
      static const int num_args = ExprT::num_args;			\
									\
      Expr(const ExprT& expr_) : expr(expr_)  {}			\
									\
      int size() const { return expr.size(); }				\
									\
      template <int Arg>						\
      bool isActive() const { return expr.template isActive<Arg>(); }	\
									\
      value_type val() const {						\
        v = expr.val();							\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      void computePartials(const value_type& bar,			\
			   value_type partials[]) const {		\
	expr.computePartials(ADJOINT, partials);			\
      }									\
									\
      void getTangents(int i, value_type dots[]) const {		\
	expr.getTangents(i, dots); }					\
									\
      template <int Arg>						\
      value_type getTangent(int i) const {				\
	return expr.template getTangent<Arg>(i);			\
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
		  bar)
FAD_UNARYOP_MACRO(operator-,
		  UnaryMinusOp, 
		  ;, 
		  -v,
		  -bar)
FAD_UNARYOP_MACRO(exp,
		  ExpOp, 
		  a = std::exp(v),
		  a,
		  bar*a)
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  ;,
		  std::log(v),
		  bar/v)
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  a = std::log(value_type(10))*v,
		  std::log10(v),
		  bar/a)
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp, 
		  a = value_type(2)*std::sqrt(v),
		  std::sqrt(v),
		  bar/a)
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  a = std::sin(v),
		  std::cos(v),
		  -bar*a)
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  a = std::cos(v),
		  std::sin(v),
		  bar*a)
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  a = std::tan(v),
		  a,
		  bar*(value_type(1)+a*a))
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  a = - std::sqrt(value_type(1)-v*v),
		  std::acos(v),
		  bar/a)
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  a = std::sqrt(value_type(1)-v*v),
		  std::asin(v),
		  bar/a)
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  a = (value_type(1)+v*v),
		  std::atan(v),
		  bar/a)
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  a = std::sinh(v),
		  std::cosh(v),
		  bar*a)
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  a = std::cosh(v),
		  std::sinh(v),
		  bar*a)
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  a = std::cosh(v),
		  std::tanh(v),
		  bar/(a*a))
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  a = std::sqrt((v-value_type(1))*(v+value_type(1))),
		  std::acosh(v),
		  bar/a)
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  a = std::sqrt(value_type(1)+v*v),
		  std::asinh(v),
		  bar/a)
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  a = value_type(1)-v*v,
		  std::atanh(v),
		  bar/a)
FAD_UNARYOP_MACRO(abs,
		  AbsOp, 
		  ;,
		  std::abs(v),
		  (expr.val() >= value_type(0.)) ? bar : value_type(-bar))
FAD_UNARYOP_MACRO(fabs,
		  FAbsOp, 
		  ;,
		  std::fabs(v),
		  (expr.val() >= value_type(0.)) ? bar : value_type(-bar))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE,LADJOINT,RADJOINT)	\
namespace Sacado {							\
  namespace ELRCacheFad {						\
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
      typedef typename ExprT1::base_expr_type base_expr_type_1;		\
      typedef typename ExprT2::base_expr_type base_expr_type_2;		\
      typedef typename ExprPromote<base_expr_type_1,			\
				   base_expr_type_2>::type base_expr_type; \
									\
      static const int num_args1 = ExprT1::num_args;			\
      static const int num_args2 = ExprT2::num_args;			\
      static const int num_args = num_args1 + num_args2;		\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      int size() const {						\
	int sz1 = expr1.size(), sz2 = expr2.size();			\
	return sz1 > sz2 ? sz1 : sz2;					\
      }									\
									\
      template <int Arg> bool isActive() const {			\
	if (Arg < num_args1)						\
	  return expr1.template isActive<Arg>();			\
	else								\
	  return expr2.template isActive<Arg-num_args1>();		\
      }									\
									\
      value_type val() const {						\
      	v1 = expr1.val();						\
	v2 = expr2.val();						\
	PARTIAL;							\
	return VALUE;							\
      }									\
									\
      void computePartials(const value_type& bar,			\
			   value_type partials[]) const {		\
	if (num_args1 > 0)						\
	  expr1.computePartials(LADJOINT, partials);			\
	if (num_args2 > 0)						\
	  expr2.computePartials(RADJOINT, partials+num_args1);		\
      }									\
									\
      void getTangents(int i, value_type dots[]) const {		\
	expr1.getTangents(i, dots);					\
	expr2.getTangents(i, dots+num_args1);				\
      }									\
									\
      template <int Arg> value_type getTangent(int i) const {		\
	if (Arg < num_args1)						\
	  return expr1.template getTangent<Arg>(i);			\
	else								\
	  return expr2.template getTangent<Arg-num_args1>(i);		\
      }									\
									\
    protected:								\
									\
      typename ExprConstRef<ExprT1>::type expr1;			\
      typename ExprConstRef<ExprT2>::type expr2;			\
      mutable value_type_1 v1;						\
      mutable value_type_2 v2;						\
      mutable value_type a;						\
      mutable value_type b;						\
									\
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
    inline Expr< OP< Expr<T>, Expr<T> > >				\
    OPNAME (const Expr<T>& expr1, const Expr<T>& expr2)			\
    {									\
      typedef OP< Expr<T>, Expr<T> > expr_t;				\
    									\
      return Expr<expr_t>(expr1, expr2);				\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< ConstExpr<typename Expr<T>::value_type>,		\
		     Expr<T> > >					\
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
		   bar,
		   bar)
FAD_BINARYOP_MACRO(operator-,
		   SubtractionOp, 
		   ;,
		   v1 - v2,
		   bar,
		   -bar)
FAD_BINARYOP_MACRO(operator*,
		   MultiplicationOp, 
		   ;,
 		   v1*v2,
		   bar*v2,
		   bar*v1)
FAD_BINARYOP_MACRO(operator/,
		   DivisionOp, 
		   ;,
		   v1/v2,
		   bar/v2,
		   -bar*v1/(v2*v2))
FAD_BINARYOP_MACRO(atan2,
		   Atan2Op,
		   a=v1*v1 + v2*v2,
		   std::atan2(v1,v2),
		   bar*v2/a,
		   -bar*v1/a)
FAD_BINARYOP_MACRO(pow,
		   PowerOp,
		   value_type c=std::pow(v1,v2); a=c*v2/v1; b=c*std::log(v1),
		   c,
		   bar*a,
		   bar*b)
FAD_BINARYOP_MACRO(max,
                   MaxOp,
		   ;,
                   std::max(v1, v2),
                   v1 >= v2 ? bar : value_type(0.),
                   v2 >  v1 ? bar : value_type(0.))
FAD_BINARYOP_MACRO(min,
                   MinOp,
		   ;,
                   std::min(v1, v2),
                   v1 <= v2 ? bar : value_type(0.),
                   v2 <  v1 ? bar : value_type(0.))

#undef FAD_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace ELRCacheFad {						\
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

  namespace ELRCacheFad {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace ELRCacheFad

} // namespace Sacado

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace ELRCacheFad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      typedef typename Expr<ExprT>::base_expr_type base_expr_type;
      return os << base_expr_type(x);
    }

  } // namespace Fad

} // namespace Sacado


#endif // SACADO_CACHEFAD_OPS_HPP
