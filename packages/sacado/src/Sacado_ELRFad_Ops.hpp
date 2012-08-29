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

#define FAD_UNARYOP_MACRO(OPNAME,OP,VALUE,ADJOINT,                      \
                          LINEAR,DX,FASTACCESSDX)			\
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
      static const bool is_linear = LINEAR;				\
									\
      Expr(const ExprT& expr_) : expr(expr_)  {}			\
									\
      int size() const { return expr.size(); }				\
									\
      template <int Arg>						\
      bool isActive() const { return expr.template isActive<Arg>(); }	\
									\
      bool isActive2(int j) const { return expr.template isActive2(j); } \
									\
      bool updateValue() const { return expr.updateValue(); }		\
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
      const value_type& getTangent(int i) const {			\
	return expr.template getTangent<Arg>(i);			\
      }									\
									\
      bool isLinear() const {                                           \
        return LINEAR;                                                  \
      }                                                                 \
                                                                        \
      bool hasFastAccess() const {                                      \
        return expr.hasFastAccess();                                    \
      }                                                                 \
                                                                        \
      const value_type dx(int i) const {                                \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      const value_type fastAccessDx(int i) const {                      \
        return FASTACCESSDX;                                            \
      }                                                                 \
									\
      const value_type* getDx(int j) const {				\
	return expr.getDx(j);						\
      }									\
									\
      const base_expr_type& getArg(int j) const {			\
	return expr.getArg(j);						\
      }									\
									\
      int numActiveArgs() const {					\
	return expr.numActiveArgs();					\
      }									\
									\
      void computeActivePartials(const value_type& bar,			\
				 value_type *partials) const {		\
        expr.computePartials(ADJOINT, partials);			\
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
		  bar,
                  true,
                  expr.dx(i),
                  expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(operator-,
		  UnaryMinusOp, 
		  -expr.val(),
		  -bar,
                  true,
                  -expr.dx(i),
                  -expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(exp,
		  ExpOp, 
		  std::exp(expr.val()),
		  bar*std::exp(expr.val()),
                  false,
                  std::exp(expr.val())*expr.dx(i),
                  std::exp(expr.val())*expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(log,
		  LogOp, 
		  std::log(expr.val()),
		  bar/expr.val(),
                  false,
                  expr.dx(i)/expr.val(),
                  expr.fastAccessDx(i)/expr.val())
FAD_UNARYOP_MACRO(log10,
		  Log10Op, 
		  std::log10(expr.val()),
		  bar/( std::log(value_type(10.))*expr.val() ),
                  false,
                  expr.dx(i)/( std::log(value_type(10))*expr.val()),
                  expr.fastAccessDx(i) / ( std::log(value_type(10))*expr.val()))
FAD_UNARYOP_MACRO(sqrt,
		  SqrtOp, 
		  std::sqrt(expr.val()),
		  value_type(0.5)*bar/std::sqrt(expr.val()),
                  false,
                  expr.dx(i)/(value_type(2)* std::sqrt(expr.val())),
                  expr.fastAccessDx(i)/(value_type(2)* std::sqrt(expr.val())))
FAD_UNARYOP_MACRO(cos,
		  CosOp, 
		  std::cos(expr.val()),
		  -bar*std::sin(expr.val()),
                  false,
                  -expr.dx(i)* std::sin(expr.val()),
                  -expr.fastAccessDx(i)* std::sin(expr.val()))
FAD_UNARYOP_MACRO(sin,
		  SinOp, 
		  std::sin(expr.val()),
		  bar*std::cos(expr.val()),
                  false,
                  expr.dx(i)* std::cos(expr.val()),
                  expr.fastAccessDx(i)* std::cos(expr.val()))
FAD_UNARYOP_MACRO(tan,
		  TanOp, 
		  std::tan(expr.val()),
		  bar*(value_type(1.)+ std::tan(expr.val())*std::tan(expr.val())),
                  false,
                  expr.dx(i)*
                    (value_type(1)+ std::tan(expr.val())* std::tan(expr.val())),
                  expr.fastAccessDx(i)*
                    (value_type(1)+ std::tan(expr.val())* std::tan(expr.val())))
FAD_UNARYOP_MACRO(acos,
		  ACosOp, 
		  std::acos(expr.val()),
		  -bar/std::sqrt(value_type(1.)-expr.val()*expr.val()),
                  false,
                  -expr.dx(i)/ std::sqrt(value_type(1)-expr.val()*expr.val()),
                  -expr.fastAccessDx(i) /
                    std::sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(asin,
		  ASinOp, 
		  std::asin(expr.val()),
		  bar/std::sqrt(value_type(1.)-expr.val()*expr.val()),
                  false,
                  expr.dx(i)/ std::sqrt(value_type(1)-expr.val()*expr.val()),
                  expr.fastAccessDx(i) /
                    std::sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atan,
		  ATanOp, 
		  std::atan(expr.val()),
		  bar/(value_type(1.)+expr.val()*expr.val()),
                  false,
                  expr.dx(i)/(value_type(1)+expr.val()*expr.val()),
                  expr.fastAccessDx(i)/(value_type(1)+expr.val()*expr.val()))
FAD_UNARYOP_MACRO(cosh,
		  CoshOp, 
		  std::cosh(expr.val()),
		  bar*std::sinh(expr.val()),
                  false,
                  expr.dx(i)* std::sinh(expr.val()),
                  expr.fastAccessDx(i)* std::sinh(expr.val()))
FAD_UNARYOP_MACRO(sinh,
		  SinhOp, 
		  std::sinh(expr.val()),
		  bar*std::cosh(expr.val()),
                  false,
                  expr.dx(i)* std::cosh(expr.val()),
                  expr.fastAccessDx(i)* std::cosh(expr.val()))
FAD_UNARYOP_MACRO(tanh,
		  TanhOp, 
		  std::tanh(expr.val()),
		  bar/(std::cosh(expr.val())*std::cosh(expr.val())),
                  false,
                  expr.dx(i)/( std::cosh(expr.val())* std::cosh(expr.val())),
                  expr.fastAccessDx(i) /
                    ( std::cosh(expr.val())* std::cosh(expr.val())))
FAD_UNARYOP_MACRO(acosh,
		  ACoshOp, 
		  std::acosh(expr.val()),
		  bar/std::sqrt((expr.val()-value_type(1.)) * 
				(expr.val()+value_type(1.))),
                  false,
                  expr.dx(i)/ std::sqrt((expr.val()-value_type(1)) *
                                       (expr.val()+value_type(1))),
                  expr.fastAccessDx(i)/ std::sqrt((expr.val()-value_type(1)) *
                                                 (expr.val()+value_type(1))))
FAD_UNARYOP_MACRO(asinh,
		  ASinhOp, 
		  std::asinh(expr.val()),
		  bar/std::sqrt(value_type(1.)+expr.val()*expr.val()),
                  false,
                  expr.dx(i)/ std::sqrt(value_type(1)+expr.val()*expr.val()),
                  expr.fastAccessDx(i)/ std::sqrt(value_type(1)+
                                                 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atanh,
		  ATanhOp, 
		  std::atanh(expr.val()),
		  bar/(value_type(1.)-expr.val()*expr.val()),
                  false,
                  expr.dx(i)/(value_type(1)-expr.val()*expr.val()),
                  expr.fastAccessDx(i)/(value_type(1)-
                                                 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(abs,
		  AbsOp, 
		  std::abs(expr.val()),
		  (expr.val() >= value_type(0.)) ? bar : value_type(-bar),
                  false,
                  expr.val() >= 0 ? value_type(+expr.dx(i)) :
                    value_type(-expr.dx(i)),
                  expr.val() >= 0 ? value_type(+expr.fastAccessDx(i)) :
                    value_type(-expr.fastAccessDx(i)))
FAD_UNARYOP_MACRO(fabs,
		  FAbsOp, 
		  std::fabs(expr.val()),
		  (expr.val() >= value_type(0.)) ? bar : value_type(-bar),
                  false,
                  expr.val() >= 0 ? value_type(+expr.dx(i)) :
                    value_type(-expr.dx(i)),
                  expr.val() >= 0 ? value_type(+expr.fastAccessDx(i)) :
                    value_type(-expr.fastAccessDx(i)))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(						\
  OPNAME,OP,VALUE,LADJOINT,RADJOINT,					\
  LINEAR,CONST_LINEAR_1, CONST_LINEAR_2,				\
  LINEAR_2,CONST_LINEAR_1_2, CONST_LINEAR_2_2,				\
  DX,FASTACCESSDX,CONST_DX_1,CONST_DX_2,				\
  CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2)		                \
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
      static const bool is_linear = LINEAR_2;				\
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
      bool isActive2(int j) const {					\
	if (j < num_args1)						\
	  return expr1.template isActive2(j);				\
	else								\
	  return expr2.template isActive2(j);				\
      }									\
									\
      bool updateValue() const {					\
	return expr1.updateValue() && expr2.updateValue();		\
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
      template <int Arg> const value_type& getTangent(int i) const {	\
	if (Arg < num_args1)						\
	  return expr1.template getTangent<Arg>(i);			\
	else								\
	  return expr2.template getTangent<Arg-num_args1>(i);		\
      }									\
									\
      bool isLinear() const {                                           \
        return LINEAR;                                                  \
      }                                                                 \
                                                                        \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess() && expr2.hasFastAccess();          \
      }                                                                 \
                                                                        \
      const value_type dx(int i) const {                                \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      const value_type fastAccessDx(int i) const {                      \
        return FASTACCESSDX;                                            \
      }                                                                 \
                                                                        \
      const value_type* getDx(int j) const {				\
	if (j < num_args1)						\
	  return expr1.getDx(j);					\
	else								\
	  return expr2.getDx(j-num_args1);				\
      }									\
									\
      const base_expr_type& getArg(int j) const {			\
	if (j < num_args1)						\
	  return expr1.getArg(j);					\
	else								\
	  return expr2.getArg(j-num_args1);				\
      }									\
									\
      int numActiveArgs() const {					\
	return expr1.numActiveArgs() + expr2.numActiveArgs();		\
      }									\
									\
      void computeActivePartials(const value_type& bar,			\
				 value_type *partials) const {		\
	if (expr1.numActiveArgs() > 0)					\
	  expr1.computePartials(LADJOINT, partials);			\
	if (expr2.numActiveArgs() > 0)					\
	  expr2.computePartials(RADJOINT, partials+expr2.numActiveArgs()); \
      }									\
    protected:								\
									\
      typename ExprConstRef<ExprT1>::type expr1;			\
      typename ExprConstRef<ExprT2>::type expr2;			\
									\
    };									\
                                                                        \
    template <typename ExprT1, typename T2>				\
    class Expr< OP<ExprT1, ConstExpr<T2> > > {				\
									\
    public:								\
									\
      typedef ConstExpr<T2> ExprT2;                                     \
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
      static const int num_args = ExprT1::num_args;			\
									\
      static const bool is_linear = CONST_LINEAR_2_2;			\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      int size() const {						\
	return expr1.size();						\
      }									\
									\
      template <int Arg> bool isActive() const {			\
	return expr1.template isActive<Arg>();				\
      }									\
									\
      bool isActive2(int j) const { return expr1.template isActive2(j); } \
									\
      bool updateValue() const {					\
	return expr1.updateValue();					\
      }									\
									\
      value_type val() const {						\
	return VALUE;							\
      }									\
									\
      void computePartials(const value_type& bar,			\
			   value_type partials[]) const {		\
	expr1.computePartials(LADJOINT, partials);			\
      }									\
									\
      void getTangents(int i, value_type dots[]) const {		\
	expr1.getTangents(i, dots);					\
      }									\
									\
      template <int Arg> const value_type& getTangent(int i) const {	\
	return expr1.template getTangent<Arg>(i);			\
      }									\
									\
      bool isLinear() const {                                           \
        return CONST_LINEAR_2;                                          \
      }                                                                 \
                                                                        \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      const value_type dx(int i) const {                                \
        return CONST_DX_2;                                              \
      }                                                                 \
                                                                        \
      const value_type fastAccessDx(int i) const {                      \
        return CONST_FASTACCESSDX_2;                                    \
      }                                                                 \
									\
      const value_type* getDx(int j) const {				\
	return expr1.getDx(j);						\
      }									\
									\
      const base_expr_type& getArg(int j) const {			\
	return expr1.getArg(j);						\
      }									\
									\
      int numActiveArgs() const {					\
	return expr1.numActiveArgs();					\
      }									\
									\
      void computeActivePartials(const value_type& bar,			\
				 value_type *partials) const {		\
        expr1.computePartials(LADJOINT, partials);			\
      }									\
                                                                        \
    protected:								\
									\
      typename ExprConstRef<ExprT1>::type expr1;			\
      typename ExprConstRef<ExprT2>::type expr2;			\
									\
    };									\
                                                                        \
    template <typename T1, typename ExprT2>				\
    class Expr< OP<ConstExpr<T1>,ExprT2> > {				\
									\
    public:								\
									\
      typedef ConstExpr<T1> ExprT1;                                     \
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
      static const int num_args = ExprT2::num_args;			\
									\
      static const bool is_linear = CONST_LINEAR_1_2;			\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      int size() const {						\
	return expr2.size();						\
      }									\
									\
      template <int Arg> bool isActive() const {			\
	return expr2.template isActive<Arg>();				\
      }									\
									\
      bool isActive2(int j) const { return expr2.template isActive2(j); } \
									\
      bool updateValue() const {					\
	return expr2.updateValue();					\
      }									\
									\
      value_type val() const {						\
	return VALUE;							\
      }									\
									\
      void computePartials(const value_type& bar,			\
			   value_type partials[]) const {		\
	expr2.computePartials(RADJOINT, partials);			\
      }									\
									\
      void getTangents(int i, value_type dots[]) const {		\
	expr2.getTangents(i, dots);					\
      }									\
									\
      template <int Arg> const value_type& getTangent(int i) const {	\
	return expr2.template getTangent<Arg>(i);			\
      }									\
									\
      bool isLinear() const {                                           \
        return CONST_LINEAR_1;                                          \
      }                                                                 \
                                                                        \
      bool hasFastAccess() const {                                      \
        return expr2.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      const value_type dx(int i) const {                                \
        return CONST_DX_1;                                              \
      }                                                                 \
                                                                        \
      const value_type fastAccessDx(int i) const {                      \
        return CONST_FASTACCESSDX_1;                                    \
      }                                                                 \
									\
      const value_type* getDx(int j) const {				\
	return expr2.getDx(j);						\
      }									\
                                                                        \
									\
      const base_expr_type& getArg(int j) const {			\
	return expr2.getArg(j);						\
      }									\
									\
      int numActiveArgs() const {					\
	return expr2.numActiveArgs();					\
      }									\
									\
      void computeActivePartials(const value_type& bar,			\
				 value_type *partials) const {		\
        expr2.computePartials(RADJOINT, partials);			\
      }									\
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
		   expr1.val() + expr2.val(),
		   bar,
		   bar,
                   expr1.isLinear() && expr2.isLinear(),
                   expr2.isLinear(),
                   expr1.isLinear(),
		   ExprT1::is_linear && ExprT2::is_linear,
		   ExprT2::is_linear,
		   ExprT1::is_linear,
                   expr1.dx(i) + expr2.dx(i),
                   expr1.fastAccessDx(i) + expr2.fastAccessDx(i),
                   expr2.dx(i),
                   expr1.dx(i),
                   expr2.fastAccessDx(i),
                   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator-,
		   SubtractionOp, 
		   expr1.val() - expr2.val(),
		   bar,
		   -bar,
                   expr1.isLinear() && expr2.isLinear(),
                   expr2.isLinear(),
                   expr1.isLinear(),
		   ExprT1::is_linear && ExprT2::is_linear,
		   ExprT2::is_linear,
		   ExprT1::is_linear,
                   expr1.dx(i) - expr2.dx(i),
                   expr1.fastAccessDx(i) - expr2.fastAccessDx(i),
                   -expr2.dx(i),
                   expr1.dx(i),
                   -expr2.fastAccessDx(i),
                   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator*,
		   MultiplicationOp, 
		   expr1.val() * expr2.val(),
		   bar*expr2.val(),
		   bar*expr1.val(),
                   false,
                   expr2.isLinear(),
                   expr1.isLinear(),
		   false,
		   ExprT2::is_linear,
		   ExprT1::is_linear,
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
		   bar/expr2.val(),
		   -bar*expr1.val()/(expr2.val()*expr2.val()),
                   false,
                   false,
                   expr1.isLinear(),
		   false,
		   false,
		   ExprT1::is_linear,
                   (expr1.dx(i)*expr2.val() - expr2.dx(i)*expr1.val()) /
                     (expr2.val()*expr2.val()),
                   (expr1.fastAccessDx(i)*expr2.val() -
                      expr2.fastAccessDx(i)*expr1.val()) /
                      (expr2.val()*expr2.val()),
                   -expr2.dx(i)*expr1.val() / (expr2.val()*expr2.val()),
                   expr1.dx(i)/expr2.val(),
                   -expr2.fastAccessDx(i)*expr1.val() / (expr2.val()*expr2.val()),
                   expr1.fastAccessDx(i)/expr2.val())
FAD_BINARYOP_MACRO(atan2,
		   Atan2Op,
		   std::atan2(expr1.val(), expr2.val()),
		   bar*expr2.val()/
		   (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
		   -bar*expr1.val()/
		   (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   false,
                   false,
                   false,
		   false,
                   false,
                   false,
                   (expr2.val()*expr1.dx(i) - expr1.val()*expr2.dx(i))/                              (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   (expr2.val()*expr1.fastAccessDx(i) - expr1.val()*expr2.fastAccessDx(i))/
                     (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   (-expr1.val()*expr2.dx(i)) / (expr1.val()*expr1.val() + expr2.val()*expr2.val()),                   
                   (expr2.val()*expr1.dx(i))/ (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   (-expr1.val()*expr2.fastAccessDx(i))/ (expr1.val()*expr1.val() + expr2.val()*expr2.val()),                   
                   (expr2.val()*expr1.fastAccessDx(i))/ (expr1.val()*expr1.val() + expr2.val()*expr2.val()))
FAD_BINARYOP_MACRO(pow,
		   PowerOp,
		   std::pow(expr1.val(), expr2.val()),
		   expr1.val() == value_type(0) ? value_type(0) : value_type(bar*std::pow(expr1.val(),expr2.val())*expr2.val()/expr1.val()),
                   expr1.val() == value_type(0) ? value_type(0) : value_type(bar*std::pow(expr1.val(),expr2.val())*std::log(expr1.val())),
                   false,
                   false,
                   false,
		   false,
                   false,
                   false,
                   expr1.val() == value_type(0) ? value_type(0) : value_type((expr2.dx(i)*std::log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0.0) : value_type((expr2.fastAccessDx(i)*std::log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0) : value_type(expr2.dx(i)*std::log(expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0.0) : value_type(expr2.val()*expr1.dx(i)/expr1.val()*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0) : value_type(expr2.fastAccessDx(i)*std::log(expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0.0) : value_type(expr2.val()*expr1.fastAccessDx(i)/expr1.val()*std::pow(expr1.val(),expr2.val())))
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   std::max(expr1.val(), expr2.val()),
                   expr1.val() >= expr2.val() ? bar : value_type(0.),
                   expr2.val() > expr1.val() ? bar : value_type(0.),
                   expr1.isLinear() && expr2.isLinear(),
                   expr2.isLinear(),
                   expr1.isLinear(),
		   ExprT1::is_linear && ExprT2::is_linear,
		   ExprT2::is_linear,
		   ExprT1::is_linear,
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
                   expr1.val() <= expr2.val() ? bar : value_type(0.),
                   expr2.val() < expr1.val() ? bar : value_type(0.),
                   expr1.isLinear() && expr2.isLinear(),
                   expr2.isLinear(),
                   expr1.isLinear(),
		   ExprT1::is_linear && ExprT2::is_linear,
		   ExprT2::is_linear,
		   ExprT1::is_linear,
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

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace ELRFad {

    template <typename ExprT>
    bool toBool(const Expr<ExprT>& x) {
      bool is_zero = (x.val() == 0.0);
      for (int i=0; i<x.size(); i++)
	is_zero = is_zero && (x.dx(i) == 0.0);
      return !is_zero;
    }

  } // namespace Fad

} // namespace Sacado

#define FAD_BOOL_MACRO(OP)						\
namespace Sacado {							\
  namespace ELRFad {							\
    template <typename ExprT1, typename ExprT2>				\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return toBool(expr1) OP toBool(expr2);				\
    }									\
									\
    template <typename ExprT2>						\
    inline bool								\
    operator OP (const typename Expr<ExprT2>::value_type& a,		\
		 const Expr<ExprT2>& expr2)				\
    {									\
      return a OP toBool(expr2);					\
    }									\
									\
    template <typename ExprT1>						\
    inline bool								\
    operator OP (const Expr<ExprT1>& expr1,				\
		 const typename Expr<ExprT1>::value_type& b)		\
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

  namespace ELRFad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      typedef typename Expr<ExprT>::base_expr_type base_expr_type;
      return os << base_expr_type(x);
    }

  } // namespace Fad

} // namespace Sacado


#endif // SACADO_FAD_OPS_HPP
