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

#ifndef SACADO_ELRFAD_OPS_HPP
#define SACADO_ELRFAD_OPS_HPP

#include "Sacado_ELRFad_Expression.hpp"
#include "Sacado_cmath.hpp"
#include <ostream>	// for std::ostream

#define FAD_UNARYOP_MACRO(OPNAME,OP,VALUE,ADJOINT)			\
namespace Sacado {							\
  namespace ELRFad {							\
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
		  expr.val(),
		  bar)
FAD_UNARYOP_MACRO(operator-,
		  UnaryMinusOp, 
		  -expr.val(),
		  -bar)
FAD_UNARYOP_MACRO(exp,
		  ExpOp, 
		  std::exp(expr.val()),
		  bar*std::exp(expr.val()))
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  std::log(expr.val()),
		  bar/expr.val())
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  std::log10(expr.val()),
		  bar/( std::log(value_type(10.))*expr.val() ))
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp, 
		  std::sqrt(expr.val()),
		  value_type(0.5)*bar/std::sqrt(expr.val()))
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  std::cos(expr.val()),
		  -bar*std::sin(expr.val()))
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  std::sin(expr.val()),
		  bar*std::cos(expr.val()))
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  std::tan(expr.val()),
		  bar*(value_type(1.)+ std::tan(expr.val())*std::tan(expr.val())))
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  std::acos(expr.val()),
		  -bar/std::sqrt(value_type(1.)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  std::asin(expr.val()),
		  bar/std::sqrt(value_type(1.)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  std::atan(expr.val()),
		  bar/(value_type(1.)+expr.val()*expr.val()))
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  std::cosh(expr.val()),
		  bar*std::sinh(expr.val()))
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  std::sinh(expr.val()),
		  bar*std::cosh(expr.val()))
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  std::tanh(expr.val()),
		  bar/(std::cosh(expr.val())*std::cosh(expr.val())))
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  std::acosh(expr.val()),
		  bar/std::sqrt((expr.val()-value_type(1.)) * 
				(expr.val()+value_type(1.))))
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  std::asinh(expr.val()),
		  bar/std::sqrt(value_type(1.)+expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  std::atanh(expr.val()),
		  bar/(value_type(1.)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(abs,
		  AbsOp, 
		  std::abs(expr.val()),
		  (expr.val() >= value_type(0.)) ? bar : value_type(-bar))
FAD_UNARYOP_MACRO(fabs,
		  FAbsOp, 
		  std::fabs(expr.val()),
		  (expr.val() >= value_type(0.)) ? bar : value_type(-bar))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,VALUE,LADJOINT,RADJOINT)		\
namespace Sacado {							\
  namespace ELRFad {							\
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
		   expr1.val() + expr2.val(),
		   bar,
		   bar)
FAD_BINARYOP_MACRO(operator-,
		   SubtractionOp, 
		   expr1.val() - expr2.val(),
		   bar,
		   -bar)
FAD_BINARYOP_MACRO(operator*,
		   MultiplicationOp, 
		   expr1.val() * expr2.val(),
		   bar*expr2.val(),
		   bar*expr1.val())
FAD_BINARYOP_MACRO(operator/,
		   DivisionOp, 
		   expr1.val() / expr2.val(),
		   bar/expr2.val(),
		   -bar*expr1.val()/(expr2.val()*expr2.val()))
FAD_BINARYOP_MACRO(atan2,
		   Atan2Op,
		   std::atan2(expr1.val(), expr2.val()),
		   bar*expr2.val()/
		   (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   -bar*expr1.val()/
		   (expr1.val()*expr1.val() + expr2.val()*expr2.val()))
FAD_BINARYOP_MACRO(pow,
		   PowerOp,
		   std::pow(expr1.val(), expr2.val()),
		   bar*std::pow(expr1.val(),expr2.val())*expr2.val()/
		   expr1.val(),
		   bar*std::pow(expr1.val(),expr2.val())*std::log(expr1.val()))
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   std::max(expr1.val(), expr2.val()),
                   expr1.val() >= expr2.val() ? bar : value_type(0.),
                   expr2.val() >= expr1.val() ? bar : value_type(0.))
FAD_BINARYOP_MACRO(min,
                   MinOp,
                   std::min(expr1.val(), expr2.val()),
                   expr1.val() <= expr2.val() ? bar : value_type(0.),
                   expr2.val() <= expr1.val() ? bar : value_type(0.))

#undef FAD_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace ELRFad {							\
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

  namespace ELRFad {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace ELRFad

} // namespace Sacado

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace ELRFad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      typedef typename Expr<ExprT>::base_expr_type base_expr_type;
      return os << base_expr_type(x);
    }

  } // namespace Fad

} // namespace Sacado


#endif // SACADO_FAD_OPS_HPP
