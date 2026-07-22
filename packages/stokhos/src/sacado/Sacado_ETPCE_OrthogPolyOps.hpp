// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_ETPCE_ORTHOGPOLY_OPS_HPP
#define SACADO_ETPCE_ORTHOGPOLY_OPS_HPP

#include "Sacado_cmath.hpp"
#include <ostream>	// for std::ostream

#define LINEAR_PCE_UNARYOP_MACRO(OPNAME,OP,OPER)			\
namespace Sacado {							\
  namespace ETPCE {							\
									\
    template <typename ExprT>						\
    class OP {};							\
									\
    template <typename ExprT>						\
    class Expr< OP<ExprT> > {						\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
      typedef typename ExprT::approx_type approx_type;			\
      typedef typename ExprT::expansion_type expansion_type;		\
      typedef typename ExprT::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT::storage_type storage_type;		\
      typedef typename ExprT::base_expr_type base_expr_type;            \
									\
      static const int num_args = ExprT::num_args;			\
									\
      Expr(const ExprT& expr_) : expr(expr_)  {}			\
									\
      std::string name() const {					\
	return std::string(#OPER) + expr.name();			\
      }									\
									\
      int size() const { return expr.size(); }				\
									\
      const approx_type& getArg(int i) const { return expr.getArg(i); }	\
      									\
      bool has_nonconst_expansion() const {				\
	return expr.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
	return expr.expansion();					\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
	return expr.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const { return expr.has_fast_access(sz); }	\
      									\
      int order() const { return expr.order(); }			\
									\
      value_type val() const { return OPER (expr.val()); }		\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return OPER (expr.fast_higher_order_coeff(i));			\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return OPER (expr.higher_order_coeff(i));			\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return OPER (expr.template eval_sample<offset,tuple_type>(x));	\
      }									\
									\
    protected:								\
									\
      const ExprT& expr;						\
									\
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

#define NONLINEAR_PCE_UNARYOP_MACRO(OPNAME,OP,OPER)			\
namespace Sacado {							\
  namespace ETPCE {							\
									\
    template <typename ExprT>						\
    class OP {};							\
									\
    template <typename ExprT>						\
    class Expr< OP<ExprT> > {						\
    public:								\
									\
      typedef typename ExprT::value_type value_type;			\
      typedef typename ExprT::approx_type approx_type;			\
      typedef typename ExprT::expansion_type expansion_type;		\
      typedef typename ExprT::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT::storage_type storage_type;		\
      typedef typename ExprT::base_expr_type base_expr_type;            \
									\
      static const int num_args = ExprT::num_args;			\
									\
      Expr(const ExprT& expr_) : expr(expr_)  {}			\
									\
      std::string name() const {					\
	return std::string(#OPER) + std::string("(") + expr.name() +	\
	  std::string(")");						\
      }									\
									\
      int size() const {						\
	if (expr.size() == 1)						\
	  return 1;							\
	else								\
	  return expansion()->size();					\
      }									\
     									\
      const approx_type& getArg(int i) const { return expr.getArg(i); }	\
									\
      bool has_nonconst_expansion() const {				\
	return expr.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
	return expr.expansion();					\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
	return expr.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const { return false; }			\
      									\
      int order() const { return size() == 0 ? 0 : 100; }		\
									\
      value_type val() const { return OPER (expr.val()); }		\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return OPER (expr.template eval_sample<offset,tuple_type>(x));	\
      }									\
									\
    protected:								\
									\
      const ExprT& expr;						\
									\
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

LINEAR_PCE_UNARYOP_MACRO(operator+, UnaryPlusOp, +)
LINEAR_PCE_UNARYOP_MACRO(operator-, UnaryMinusOp, -)

NONLINEAR_PCE_UNARYOP_MACRO(exp, ExpOp, std::exp)
NONLINEAR_PCE_UNARYOP_MACRO(log, LogOp, std::log)
NONLINEAR_PCE_UNARYOP_MACRO(log10, Log10Op, std::log10)
NONLINEAR_PCE_UNARYOP_MACRO(sqrt, SqrtOp, std::sqrt)
NONLINEAR_PCE_UNARYOP_MACRO(cbrt, CbrtOp, std::cbrt)
NONLINEAR_PCE_UNARYOP_MACRO(cos, CosOp, std::cos)
NONLINEAR_PCE_UNARYOP_MACRO(sin, SinOp, std::sin)
NONLINEAR_PCE_UNARYOP_MACRO(tan, TanOp, std::tan)
NONLINEAR_PCE_UNARYOP_MACRO(acos, ACosOp, std::acos)
NONLINEAR_PCE_UNARYOP_MACRO(asin, ASinOp, std::asin)
NONLINEAR_PCE_UNARYOP_MACRO(atan, ATanOp, std::atan)
NONLINEAR_PCE_UNARYOP_MACRO(cosh, CoshOp, std::cosh)
NONLINEAR_PCE_UNARYOP_MACRO(sinh, SinhOp, std::sinh)
NONLINEAR_PCE_UNARYOP_MACRO(tanh, TanhOp, std::tanh)
NONLINEAR_PCE_UNARYOP_MACRO(acosh, ACoshOp, std::acosh)
NONLINEAR_PCE_UNARYOP_MACRO(asinh, ASinhOp, std::asinh)
NONLINEAR_PCE_UNARYOP_MACRO(atanh, ATanhOp, std::atanh)
NONLINEAR_PCE_UNARYOP_MACRO(abs, AbsOp, std::abs)
NONLINEAR_PCE_UNARYOP_MACRO(fabs, FAbsOp, std::fabs)

#undef LINEAR_PCE_UNARYOP_MACRO
#undef NONLINEAR_PCE_UNARYOP_MACRO

#define LINEAR_PCE_BINARYOP_MACRO(OPNAME,OP,OPER)			\
namespace Sacado {							\
  namespace ETPCE {							\
									\
    template <typename ExprT1, typename ExprT2>				\
    class OP {};							\
									\
    template <typename T1, typename T2>					\
    class Expr< OP< Expr<T1>, Expr<T2> > > {				\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT1::value_type value_type_1;			\
      typedef typename ExprT2::value_type value_type_2;			\
      typedef typename Sacado::Promote<value_type_1,			\
				       value_type_2>::type value_type;  \
									\
      typedef typename ExprT1::approx_type approx_type;			\
      typedef typename ExprT1::expansion_type expansion_type;		\
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
      static const int num_args1 = ExprT1::num_args;			\
      static const int num_args2 = ExprT2::num_args;			\
      static const int num_args = num_args1 + num_args2;		\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
      									\
      std::string name() const {					\
	return expr1.name() + std::string(#OPER) + expr2.name();	\
      }									\
									\
      int size() const {						\
	int sz1 = expr1.size(), sz2 = expr2.size();			\
	return sz1 > sz2 ? sz1 : sz2;					\
      }									\
									\
      const approx_type& getArg(int i) const {				\
	if (i < num_args1)						\
	  return expr1.getArg(i);					\
	else								\
	  return expr2.getArg(i-num_args1);				\
      }									\
									\
      bool has_nonconst_expansion() const {				\
	return expr1.has_nonconst_expansion() ||			\
	  expr2.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
        return expr1.has_nonconst_expansion() ? expr1.expansion() :	\
          expr2.expansion();						\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
        return expr1.quad_expansion() != Teuchos::null ?		\
	  expr1.quad_expansion() :					\
          expr2.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const {					\
	return expr1.has_fast_access(sz) && expr2.has_fast_access(sz);	\
      }									\
									\
      int order() const {						\
	int o1 = expr1.order(), o2 = expr2.order();			\
	return o1 > o2 ? o1 : o2;					\
      }									\
									\
      value_type val() const {						\
	return expr1.val() OPER expr2.val();				\
      }									\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return expr1.fast_higher_order_coeff(i) OPER			\
	  expr2.fast_higher_order_coeff(i);				\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return expr1.higher_order_coeff(i) OPER				\
	  expr2.higher_order_coeff(i);					\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return expr1.template eval_sample<offset,tuple_type>(x) OPER	\
	  expr2.template eval_sample<offset+num_args1,tuple_type>(x);	\
      }									\
      									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ExprT2& expr2;						\
									\
    };									\
									\
    template <typename T1>						\
    class Expr< OP< Expr<T1>, typename Expr<T1>::value_type> > {	\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef typename ExprT1::value_type value_type;			\
      typedef typename ExprT1::value_type ConstT;			\
									\
      typedef typename ExprT1::approx_type approx_type;			\
      typedef typename ExprT1::expansion_type expansion_type;		\
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
      static const int num_args = ExprT1::num_args;			\
									\
      Expr(const ExprT1& expr1_, const ConstT& c_) :			\
	expr1(expr1_), c(c_) {}						\
									\
      std::string name() const {					\
	return expr1.name() + std::string(#OPER) + std::string("c");	\
      }									\
									\
      int size() const { return expr1.size(); }				\
									\
      const approx_type& getArg(int i) const {				\
	return expr1.getArg(i);						\
      }									\
									\
      bool has_nonconst_expansion() const {				\
	return expr1.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
	return expr1.expansion();					\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
	return expr1.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const { return expr1.has_fast_access(sz); }	\
									\
      int order() const { return expr1.order(); }			\
									\
      value_type val() const {						\
	return expr1.val() OPER c;					\
      }									\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return expr1.fast_higher_order_coeff(i);			\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return expr1.higher_order_coeff(i);				\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return expr1.template eval_sample<offset,tuple_type>(x) OPER c; \
      }									\
									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ConstT& c;							\
    };									\
									\
    template <typename T2>						\
    class Expr< OP< typename Expr<T2>::value_type, Expr<T2> > > {	\
									\
    public:								\
									\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT2::value_type value_type;			\
      typedef typename ExprT2::value_type ConstT;			\
									\
      typedef typename ExprT2::approx_type approx_type;			\
      typedef typename ExprT2::expansion_type expansion_type;		\
      typedef typename ExprT2::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT2::storage_type storage_type;		\
      typedef typename ExprT2::base_expr_type base_expr_type;           \
									\
      static const int num_args = ExprT2::num_args;			\
									\
      Expr(const ConstT& c_, const ExprT2& expr2_) :			\
	c(c_), expr2(expr2_) {}						\
									\
      std::string name() const {					\
	return std::string("c") + std::string(#OPER) + expr2.name();	\
      }									\
									\
      int size() const { return expr2.size(); }				\
									\
     const approx_type& getArg(int i) const { return expr2.getArg(i); } \
									\
      bool has_nonconst_expansion() const {				\
	return expr2.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
	return expr2.expansion();					\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
	return expr2.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const { return expr2.has_fast_access(sz); }	\
									\
      int order() const { return expr2.order(); }			\
									\
      value_type val() const {						\
	return c OPER expr2.val();					\
      }									\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return  OPER expr2.fast_higher_order_coeff(i);			\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return  OPER expr2.higher_order_coeff(i);			\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return c OPER expr2.template eval_sample<offset,tuple_type>(x); \
      }									\
      									\
    protected:								\
									\
      const ConstT& c;							\
      const ExprT2& expr2;						\
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
    inline Expr< OP< typename Expr<T>::value_type, Expr<T> > >		\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< ConstT, Expr<T> > expr_t;				\
									\
      return Expr<expr_t>(c, expr);					\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>, typename Expr<T>::value_type > >		\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< Expr<T>, ConstT > expr_t;				\
									\
      return Expr<expr_t>(expr, c);					\
    }									\
  }									\
}

#define NONLINEAR_PCE_BINARYOP_MACRO(OPNAME,OP,OPER)			\
namespace Sacado {							\
  namespace ETPCE {							\
									\
    template <typename ExprT1, typename ExprT2>				\
    class OP {};							\
									\
    template <typename T1, typename T2>					\
    class Expr< OP< Expr<T1>, Expr<T2> > > {				\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT1::value_type value_type_1;			\
      typedef typename ExprT2::value_type value_type_2;			\
      typedef typename Sacado::Promote<value_type_1,			\
				       value_type_2>::type value_type;  \
									\
      typedef typename ExprT1::approx_type approx_type;			\
      typedef typename ExprT1::expansion_type expansion_type;		\
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
      static const int num_args1 = ExprT1::num_args;			\
      static const int num_args2 = ExprT2::num_args;			\
      static const int num_args = num_args1 + num_args2;		\
									\
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :		\
	expr1(expr1_), expr2(expr2_) {}					\
									\
      std::string name() const {					\
	return std::string(#OPER) + std::string("(") + expr1.name() +	\
	  std::string(",") + expr2.name() + std::string(")");		\
      }									\
									\
      int size() const {						\
	if (expr1.size() == 1 && expr2.size() == 1)			\
	  return 1;							\
	else								\
	  return expansion()->size();					\
      }									\
									\
      const approx_type& getArg(int i) const {				\
	if (i < num_args1)						\
	  return expr1.getArg(i);					\
	else								\
	  return expr2.getArg(i-num_args1);				\
      }									\
									\
      bool has_nonconst_expansion() const {				\
	return expr1.has_nonconst_expansion() ||			\
	  expr2.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
	return expr1.has_nonconst_expansion() ? expr1.expansion() :	\
	  expr2.expansion();						\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
        return expr1.quad_expansion() != Teuchos::null ?		\
	  expr1.quad_expansion() :					\
          expr2.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const { return false; }			\
									\
      int order() const { return size() == 0 ? 0 : 100; }		\
									\
      value_type val() const {						\
	return OPER (expr1.val(), expr2.val());				\
      }									\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return OPER (expr1.template eval_sample<offset,tuple_type>(x),	\
		     expr2.template eval_sample<offset+num_args1,tuple_type>(x)); \
      }									\
      									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ExprT2& expr2;						\
									\
    };									\
									\
    template <typename T1>						\
    class Expr< OP< Expr<T1>, typename Expr<T1>::value_type> > {	\
									\
    public:								\
									\
      typedef Expr<T1> ExprT1;						\
      typedef typename ExprT1::value_type value_type;			\
      typedef typename ExprT1::value_type ConstT;			\
									\
      typedef typename ExprT1::approx_type approx_type;			\
      typedef typename ExprT1::expansion_type expansion_type;		\
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT1::storage_type storage_type;		\
      typedef typename ExprT1::base_expr_type base_expr_type;           \
									\
      static const int num_args = ExprT1::num_args;			\
									\
      Expr(const ExprT1& expr1_, const ConstT& c_) :			\
	expr1(expr1_), c(c_) {}						\
									\
      std::string name() const {					\
	return std::string(#OPER) + std::string("(") + expr1.name() +	\
	  std::string(",c)");						\
      }									\
									\
      int size() const {						\
	if (expr1.size() == 1)						\
	  return 1;							\
	else								\
	  return expansion()->size();					\
      }									\
									\
      const approx_type& getArg(int i) const {				\
	return expr1.getArg(i);						\
      }									\
									\
      bool has_nonconst_expansion() const {				\
	return expr1.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
	return expr1.expansion();					\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
	return expr1.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const { return false; }			\
									\
      int order() const { return size() == 0 ? 0 : 100; }		\
									\
      value_type val() const {						\
	return OPER (expr1.val(), c);					\
      }									\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return OPER (expr1.template eval_sample<offset,tuple_type>(x), c); \
      }									\
									\
    protected:								\
									\
      const ExprT1& expr1;						\
      const ConstT& c;							\
    };									\
									\
    template <typename T2>						\
    class Expr< OP< typename Expr<T2>::value_type, Expr<T2> > > {	\
									\
    public:								\
									\
      typedef Expr<T2> ExprT2;						\
      typedef typename ExprT2::value_type value_type;			\
      typedef typename ExprT2::value_type ConstT;			\
									\
      typedef typename ExprT2::approx_type approx_type;			\
      typedef typename ExprT2::expansion_type expansion_type;		\
      typedef typename ExprT2::quad_expansion_type quad_expansion_type;	\
      typedef typename ExprT2::storage_type storage_type;		\
      typedef typename ExprT2::base_expr_type base_expr_type;           \
									\
      static const int num_args = ExprT2::num_args;			\
									\
      Expr(const ConstT& c_, const ExprT2& expr2_) :			\
	c(c_), expr2(expr2_) {}						\
									\
      std::string name() const {					\
	return std::string(#OPER) + std::string("(c,") +		\
	  expr2.name() + std::string(")");				\
      }									\
									\
      int size() const {						\
	if (expr2.size() == 1)						\
	  return 1;							\
	else								\
	  return expansion()->size();					\
      }									\
									\
     const approx_type& getArg(int i) const { return expr2.getArg(i); } \
									\
      Teuchos::RCP<expansion_type> expansion() const {			\
	return expr2.expansion();					\
      }									\
									\
      bool has_nonconst_expansion() const {				\
	return expr2.has_nonconst_expansion();				\
      }									\
									\
      Teuchos::RCP<quad_expansion_type> quad_expansion() const {	\
	return expr2.quad_expansion();					\
      }									\
									\
      bool has_fast_access(int sz) const { return false; }			\
									\
      int order() const { return size() == 0 ? 0 : 100; }		\
									\
      value_type val() const {						\
	return OPER (c, expr2.val());					\
      }									\
									\
      value_type fast_higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      value_type higher_order_coeff(int i) const {			\
	return value_type(0);						\
      }									\
									\
      template <int offset, typename tuple_type>			\
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {	\
	return OPER (c, expr2.template eval_sample<offset,tuple_type>(x)); \
      }									\
      									\
    protected:								\
									\
      const ConstT& c;							\
      const ExprT2& expr2;						\
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
    inline Expr< OP< typename Expr<T>::value_type, Expr<T> > >		\
    OPNAME (const typename Expr<T>::value_type& c,			\
	    const Expr<T>& expr)					\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< ConstT, Expr<T> > expr_t;				\
									\
      return Expr<expr_t>(c, expr);					\
    }									\
									\
    template <typename T>						\
    inline Expr< OP< Expr<T>, typename Expr<T>::value_type > >		\
    OPNAME (const Expr<T>& expr,					\
	    const typename Expr<T>::value_type& c)			\
    {									\
      typedef typename Expr<T>::value_type ConstT;			\
      typedef OP< Expr<T>, ConstT > expr_t;				\
									\
      return Expr<expr_t>(expr, c);					\
    }									\
  }									\
}

LINEAR_PCE_BINARYOP_MACRO(operator+, AdditionOp, +)
LINEAR_PCE_BINARYOP_MACRO(operator-, SubtractionOp, -)

NONLINEAR_PCE_BINARYOP_MACRO(atan2, Atan2Op, std::atan2)
NONLINEAR_PCE_BINARYOP_MACRO(pow, PowerOp, std::pow)
NONLINEAR_PCE_BINARYOP_MACRO(max, MaxOp, std::max)
NONLINEAR_PCE_BINARYOP_MACRO(min, MinOp, std::min)

#undef LINEAR_PCE_BINARYOP_MACRO
#undef NONLINEAR_PCE_BINARYOP_MACRO

//-------------------------- Multiplication Operator -----------------------

namespace Sacado {
  namespace ETPCE {

    template <typename ExprT1, typename ExprT2>
    class MultiplicationOp {};

    template <typename T1, typename T2>
    class Expr< MultiplicationOp< Expr<T1>, Expr<T2> > > {

    public:

      typedef Expr<T1> ExprT1;
      typedef Expr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;  

      typedef typename ExprT1::approx_type approx_type;
      typedef typename ExprT1::expansion_type expansion_type;
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;
      typedef typename ExprT1::storage_type storage_type;
      typedef typename ExprT1::base_expr_type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}
      
      std::string name() const {
	return expr1.name() + std::string("*") + expr2.name();
      }

      int size() const {
	int sz1 = expr1.size();
	int sz2 = expr2.size();
	if (sz1 > 1 && sz2 > 1)
	  return expansion()->size();
	else
	  return sz1*sz2;
      }

      const approx_type& getArg(int i) const {
	if (i < num_args1)
	  return expr1.getArg(i);
	else
	  return expr2.getArg(i-num_args1);
      }

      bool has_nonconst_expansion() const {
	return expr1.has_nonconst_expansion() || expr2.has_nonconst_expansion();
      }

      Teuchos::RCP<expansion_type> expansion() const {
	return expr1.has_nonconst_expansion() ? expr1.expansion() : 
	  expr2.expansion();
      }

      Teuchos::RCP<quad_expansion_type> quad_expansion() const {
        return expr1.quad_expansion() != Teuchos::null ?
	  expr1.quad_expansion() :
          expr2.quad_expansion();
      }

      bool has_fast_access(int sz) const {
	return expr1.has_fast_access(sz) && expr2.has_fast_access(sz);
      }

      int order() const { return expr1.order() + expr2.order(); }

      value_type val() const {
	if (order() == 0)
	  return expr1.val() * expr2.val();
	else
	  return quad_expansion()->compute_times_coeff(0,expr1,expr2);
      }

      value_type fast_higher_order_coeff(int i) const {
	return quad_expansion()->fast_compute_times_coeff(i,expr1,expr2);
      }

      value_type higher_order_coeff(int i) const {
	return quad_expansion()->compute_times_coeff(i,expr1,expr2);
      }

      template <int offset, typename tuple_type>
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {
	return expr1.template eval_sample<offset,tuple_type>(x) * 
	  expr2.template eval_sample<offset+num_args1,tuple_type>(x);
      }
      
    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename T1>
    class Expr< MultiplicationOp< Expr<T1>, typename Expr<T1>::value_type> > {

    public:

      typedef Expr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::value_type ConstT;

      typedef typename ExprT1::approx_type approx_type;
      typedef typename ExprT1::expansion_type expansion_type;
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;
      typedef typename ExprT1::storage_type storage_type;
      typedef typename ExprT1::base_expr_type base_expr_type;

      static const int num_args = ExprT1::num_args;

      Expr(const ExprT1& expr1_, const ConstT& c_) :
	expr1(expr1_), c(c_) {}

      std::string name() const {
	return expr1.name() + std::string("*c");
      }

      int size() const { return expr1.size(); }

      const approx_type& getArg(int i) const {
	return expr1.getArg(i);
      }
      
      bool has_nonconst_expansion() const {
	return expr1.has_nonconst_expansion();
      }

      Teuchos::RCP<expansion_type> expansion() const {
	return expr1.expansion();
      }

      Teuchos::RCP<quad_expansion_type> quad_expansion() const {
	return expr1.quad_expansion();
      }

      bool has_fast_access(int sz) const { return expr1.has_fast_access(sz); }

      int order() const { return expr1.order(); }

      value_type val() const {
	return expr1.val() * c;
      }

      value_type fast_higher_order_coeff(int i) const {
	return expr1.fast_higher_order_coeff(i) * c;
      }

      value_type higher_order_coeff(int i) const {
	return expr1.higher_order_coeff(i) * c;
      }

      template <int offset, typename tuple_type>
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {
	return expr1.template eval_sample<offset,tuple_type>(x) * c;
      }

    protected:

      const ExprT1& expr1;
      const ConstT& c;
    };

    template <typename T2>
    class Expr< MultiplicationOp< typename Expr<T2>::value_type, Expr<T2> > > {

    public:

      typedef Expr<T2> ExprT2;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::value_type ConstT;

      typedef typename ExprT2::approx_type approx_type;
      typedef typename ExprT2::expansion_type expansion_type;
      typedef typename ExprT2::quad_expansion_type quad_expansion_type;
      typedef typename ExprT2::storage_type storage_type;
      typedef typename ExprT2::base_expr_type base_expr_type;

      static const int num_args = ExprT2::num_args;

      Expr(const ConstT& c_, const ExprT2& expr2_) :
	c(c_), expr2(expr2_) {}

      std::string name() const {
	return std::string("c*") + expr2.name();
      }

      int size() const { return expr2.size(); }

      const approx_type& getArg(int i) const { return expr2.getArg(i); } 

      bool has_nonconst_expansion() const {
	return expr2.has_nonconst_expansion();
      }

      Teuchos::RCP<expansion_type> expansion() const {
	return expr2.expansion();
      }

      Teuchos::RCP<quad_expansion_type> quad_expansion() const {
	return expr2.quad_expansion();
      }

      bool has_fast_access(int sz) const { return expr2.has_fast_access(sz); }

      int order() const { return expr2.order(); }

      value_type val() const {
	return c * expr2.val();
      }

      value_type fast_higher_order_coeff(int i) const {
	return c * expr2.fast_higher_order_coeff(i);
      }

      value_type higher_order_coeff(int i) const {
	return c * expr2.higher_order_coeff(i);
      }

      template <int offset, typename tuple_type>
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {
	return c * expr2.template eval_sample<offset,tuple_type>(x);
      }
      
    protected:

      const ConstT& c;
      const ExprT2& expr2;
    };

    template <typename T1, typename T2>
    inline Expr< MultiplicationOp< Expr<T1>, Expr<T2> > >
    operator* (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef MultiplicationOp< Expr<T1>, Expr<T2> > expr_t;
    
      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    inline Expr< MultiplicationOp< Expr<T>, Expr<T> > >
    operator* (const Expr<T>& expr1, const Expr<T>& expr2)
    {
      typedef MultiplicationOp< Expr<T>, Expr<T> > expr_t;
    
      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    inline Expr< MultiplicationOp< typename Expr<T>::value_type, Expr<T> > >
    operator* (const typename Expr<T>::value_type& c,
	    const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type ConstT;
      typedef MultiplicationOp< ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(c, expr);
    }

    template <typename T>
    inline Expr< MultiplicationOp< Expr<T>, typename Expr<T>::value_type > >
    operator* (const Expr<T>& expr,
	    const typename Expr<T>::value_type& c)
    {
      typedef typename Expr<T>::value_type ConstT;
      typedef MultiplicationOp< Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(expr, c);
    }
  }
}

//-------------------------- Division Operator -----------------------

namespace Sacado {
  namespace ETPCE {

    template <typename ExprT1, typename ExprT2>
    class DivisionOp {};

    template <typename T1, typename T2>
    class Expr< DivisionOp< Expr<T1>, Expr<T2> > > {

    public:

      typedef Expr<T1> ExprT1;
      typedef Expr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
				       value_type_2>::type value_type;  

      typedef typename ExprT1::approx_type approx_type;
      typedef typename ExprT1::expansion_type expansion_type;
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;
      typedef typename ExprT1::storage_type storage_type;
      typedef typename ExprT1::base_expr_type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
	expr1(expr1_), expr2(expr2_) {}

      std::string name() const {
	return expr1.name() + std::string("/") + expr2.name();
      }

      int size() const {
	if (expr2.size() == 1)
	  return expr1.size();
	else
	  return expansion()->size();
      }

      const approx_type& getArg(int i) const {
	if (i < num_args1)
	  return expr1.getArg(i);
	else
	  return expr2.getArg(i-num_args1);
      }

      bool has_nonconst_expansion() const {
	return expr1.has_nonconst_expansion() || expr2.has_nonconst_expansion();
      }

      Teuchos::RCP<expansion_type> expansion() const {
	return expr1.has_nonconst_expansion() ? expr1.expansion() : 
	  expr2.expansion();
      }

      Teuchos::RCP<quad_expansion_type> quad_expansion() const {
        return expr1.quad_expansion() != Teuchos::null ?
	  expr1.quad_expansion() :
          expr2.quad_expansion();
      }

      bool has_fast_access(int sz) const {
	return expr1.has_fast_access(sz) && (expr2.order() == 0); 
      }

      int order() const { return expr2.order() == 0 ? expr1.order() : 100; }

      value_type val() const {
	return expr1.val() / expr2.val();
      }

      value_type fast_higher_order_coeff(int i) const {
	return expr1.fast_higher_order_coeff(i) / expr2.val();
      }

      value_type higher_order_coeff(int i) const {
	return expr1.higher_order_coeff(i) / expr2.val();
      }

      template <int offset, typename tuple_type>
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {
	return expr1.template eval_sample<offset,tuple_type>(x) / 
	  expr2.template eval_sample<offset+num_args1,tuple_type>(x);
      }
      
    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename T1>
    class Expr< DivisionOp< Expr<T1>, typename Expr<T1>::value_type> > {

    public:

      typedef Expr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::value_type ConstT;

      typedef typename ExprT1::approx_type approx_type;
      typedef typename ExprT1::expansion_type expansion_type;
      typedef typename ExprT1::quad_expansion_type quad_expansion_type;
      typedef typename ExprT1::storage_type storage_type;
      typedef typename ExprT1::base_expr_type base_expr_type;

      static const int num_args = ExprT1::num_args;

      Expr(const ExprT1& expr1_, const ConstT& c_) :
	expr1(expr1_), c(c_) {}

      std::string name() const {
	return expr1.name() + std::string("/c");
      }

      int size() const { return expr1.size(); }

      const approx_type& getArg(int i) const {
	return expr1.getArg(i);
      }

      bool has_nonconst_expansion() const {
	return expr1.has_nonconst_expansion();
      }

      Teuchos::RCP<expansion_type> expansion() const {
	return expr1.expansion();
      }

      Teuchos::RCP<quad_expansion_type> quad_expansion() const {
	return expr1.quad_expansion();
      }

      bool has_fast_access(int sz) const { return expr1.has_fast_access(sz); }

      int order() const { return expr1.order(); }

      value_type val() const {
	return expr1.val() / c;
      }

      value_type fast_higher_order_coeff(int i) const {
	return expr1.fast_higher_order_coeff(i) / c;
      }

      value_type higher_order_coeff(int i) const {
	return expr1.higher_order_coeff(i) / c;
      }

      template <int offset, typename tuple_type>
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {
	return expr1.template eval_sample<offset,tuple_type>(x) / c;
      }

    protected:

      const ExprT1& expr1;
      const ConstT& c;
    };

    template <typename T2>
    class Expr< DivisionOp< typename Expr<T2>::value_type, Expr<T2> > > {

    public:

      typedef Expr<T2> ExprT2;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::value_type ConstT;

      typedef typename ExprT2::approx_type approx_type;
      typedef typename ExprT2::expansion_type expansion_type;
      typedef typename ExprT2::quad_expansion_type quad_expansion_type;
      typedef typename ExprT2::storage_type storage_type;
      typedef typename ExprT2::base_expr_type base_expr_type;

      static const int num_args = ExprT2::num_args;

      Expr(const ConstT& c_, const ExprT2& expr2_) :
	c(c_), expr2(expr2_) {}

      std::string name() const {
	return std::string("c/") + expr2.name();
      }

      int size() const {
	if (expr2.size() == 1)
	  return 1;
	else
	  return expansion()->size();
      }

     const approx_type& getArg(int i) const { return expr2.getArg(i); } 

      bool has_nonconst_expansion() const {
	return expr2.has_nonconst_expansion();
      }

      Teuchos::RCP<expansion_type> expansion() const {
	return expr2.expansion();
      }

      Teuchos::RCP<quad_expansion_type> quad_expansion() const {
	return expr2.quad_expansion();
      }

      bool has_fast_access(int sz) const { return false; }

      int order() const { return expr2.order() == 0 ? 0 : 100; }

      value_type val() const {
	return c / expr2.val();
      }

      value_type fast_higher_order_coeff(int i) const {
	return value_type(0);
      }

      value_type higher_order_coeff(int i) const {
	return value_type(0);
      }

      template <int offset, typename tuple_type>
      KERNEL_PREFIX value_type eval_sample(tuple_type x) const {
	return c / expr2.template eval_sample<offset,tuple_type>(x);
      }
      
    protected:

      const ConstT& c;
      const ExprT2& expr2;
    };

    template <typename T1, typename T2>
    inline Expr< DivisionOp< Expr<T1>, Expr<T2> > >
    operator/ (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef DivisionOp< Expr<T1>, Expr<T2> > expr_t;
    
      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    inline Expr< DivisionOp< Expr<T>, Expr<T> > >
    operator/ (const Expr<T>& expr1, const Expr<T>& expr2)
    {
      typedef DivisionOp< Expr<T>, Expr<T> > expr_t;
    
      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    inline Expr< DivisionOp< typename Expr<T>::value_type, Expr<T> > >
    operator/ (const typename Expr<T>::value_type& c,
	    const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type ConstT;
      typedef DivisionOp< ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(c, expr);
    }

    template <typename T>
    inline Expr< DivisionOp< Expr<T>, typename Expr<T>::value_type > >
    operator/ (const Expr<T>& expr,
	    const typename Expr<T>::value_type& c)
    {
      typedef typename Expr<T>::value_type ConstT;
      typedef DivisionOp< Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(expr, c);
    }
  }
}

//-------------------------- Relational Operators -----------------------

#define PCE_RELOP_MACRO(OP)						\
namespace Sacado {							\
  namespace ETPCE {							\
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

PCE_RELOP_MACRO(==)
PCE_RELOP_MACRO(!=)
PCE_RELOP_MACRO(<)
PCE_RELOP_MACRO(>)
PCE_RELOP_MACRO(<=)
PCE_RELOP_MACRO(>=)
PCE_RELOP_MACRO(<<=)
PCE_RELOP_MACRO(>>=)
PCE_RELOP_MACRO(&)
PCE_RELOP_MACRO(|)

#undef PCE_RELOP_MACRO

namespace Sacado {

  namespace ETPCE {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr) 
    {
      return ! expr.val();
    }

  } // namespace ETPCE

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace ETPCE {

    template <typename ExprT>
    bool toBool(const Expr<ExprT>& x) {
      bool is_zero = true;
      for (int i=0; i<x.size(); i++)
	is_zero = is_zero && (x.coeff(i) == 0.0);
      return !is_zero;
    }

  } // namespace ETPCE

} // namespace Sacado

#define PCE_BOOL_MACRO(OP)						\
namespace Sacado {							\
  namespace ETPCE {							\
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

PCE_BOOL_MACRO(&&)
PCE_BOOL_MACRO(||)

#undef PCE_BOOL_MACRO

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace ETPCE {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::storage_type storage_type;
      OrthogPoly<value_type, storage_type> a(x);
      os << a;
      return os;
    }

  } // namespace ETPCE

} // namespace Sacado




#endif // SACADO_ETPCE_ORTHOGPOLY_OPS_HPP
