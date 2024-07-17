// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
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
#include <ostream>      // for std::ostream

#define FAD_UNARYOP_MACRO(OPNAME,OP,VALUE,ADJOINT,                      \
                          LINEAR,DX,FASTACCESSDX)                       \
namespace Sacado {                                                      \
  namespace ELRFad {                                                    \
                                                                        \
    template <typename ExprT>                                           \
    class OP {};                                                        \
                                                                        \
    template <typename ExprT>                                           \
    class Expr< OP<ExprT> > {                                           \
    public:                                                             \
                                                                        \
      typedef typename ExprT::value_type value_type;                    \
      typedef typename ExprT::scalar_type scalar_type;                  \
      typedef typename ExprT::base_expr_type base_expr_type;            \
                                                                        \
      static const int num_args = ExprT::num_args;                      \
                                                                        \
      static const bool is_linear = LINEAR;                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}               \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const { return expr.size(); }                          \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive() const { return expr.template isActive<Arg>(); }   \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive2(int j) const { return expr.isActive2(j); }         \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool updateValue() const { return expr.updateValue(); }           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {}                                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computePartials(const value_type& bar,                       \
                           value_type partials[]) const {               \
        expr.computePartials(ADJOINT, partials);                        \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void getTangents(int i, value_type dots[]) const {                \
        expr.getTangents(i, dots); }                                    \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      const value_type& getTangent(int i) const {                       \
        return expr.template getTangent<Arg>(i);                        \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isLinear() const {                                           \
        return LINEAR;                                                  \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr.hasFastAccess();                                    \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        return FASTACCESSDX;                                            \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type* getDx(int j) const {                            \
        return expr.getDx(j);                                           \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int numActiveArgs() const {                                       \
        return expr.numActiveArgs();                                    \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computeActivePartials(const value_type& bar,                 \
                                 value_type *partials) const {          \
        expr.computePartials(ADJOINT, partials);                        \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ExprT& expr;                                                \
    };                                                                  \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< OP< Expr<T> > >                                               \
    OPNAME (const Expr<T>& expr)                                        \
    {                                                                   \
      typedef OP< Expr<T> > expr_t;                                     \
                                                                        \
      return Expr<expr_t>(expr);                                        \
    }                                                                   \
  }                                                                     \
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
FAD_UNARYOP_MACRO(safe_sqrt,
                  SafeSqrtOp,
                  std::sqrt(expr.val()),
                  expr.val() == value_type(0.0) ? value_type(0.0) : value_type(value_type(0.5)*bar/std::sqrt(expr.val())),
                  false,
                  expr.val() == value_type(0.0) ? value_type(0.0) : value_type(expr.dx(i)/(value_type(2)*std::sqrt(expr.val()))),
                  expr.val() == value_type(0.0) ? value_type(0.0) : value_type(expr.fastAccessDx(i)/(value_type(2)*std::sqrt(expr.val()))))
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
                  bar*(value_type(1)-std::tanh(expr.val())*std::tanh(expr.val())),
                  false,
                  expr.dx(i)*(value_type(1)-std::tanh(expr.val())*std::tanh(expr.val())),
                  expr.fastAccessDx(i)*(value_type(1)-std::tanh(expr.val())*std::tanh(expr.val())))
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
FAD_UNARYOP_MACRO(cbrt,
                  CbrtOp,
                  std::cbrt(expr.val()),
                  bar/(value_type(3)*std::cbrt(expr.val()*expr.val())),
                  false,
                  expr.dx(i)/(value_type(3)*std::cbrt(expr.val()*expr.val())),
                  expr.fastAccessDx(i)/(value_type(3)*std::cbrt(expr.val()*expr.val())))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(                                             \
  OPNAME,OP,VALUE,LADJOINT,RADJOINT,                                    \
  LINEAR,CONST_LINEAR_1, CONST_LINEAR_2,                                \
  LINEAR_2,CONST_LINEAR_1_2, CONST_LINEAR_2_2,                          \
  DX,FASTACCESSDX,CONST_DX_1,CONST_DX_2,                                \
  CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2)                            \
namespace Sacado {                                                      \
  namespace ELRFad {                                                    \
                                                                        \
    template <typename ExprT1, typename ExprT2>                         \
    class OP {};                                                        \
                                                                        \
    template <typename ExprT1, typename ExprT2>                         \
    class Expr< OP<ExprT1,ExprT2> > {                                   \
                                                                        \
    public:                                                             \
                                                                        \
      typedef typename ExprT1::value_type value_type_1;                 \
      typedef typename ExprT2::value_type value_type_2;                 \
      typedef typename Sacado::Promote<value_type_1,                    \
                                       value_type_2>::type value_type;  \
                                                                        \
      typedef typename ExprT1::scalar_type scalar_type_1;               \
      typedef typename ExprT2::scalar_type scalar_type_2;               \
      typedef typename Sacado::Promote<scalar_type_1,                   \
                                       scalar_type_2>::type scalar_type; \
                                                                        \
      typedef typename ExprT1::base_expr_type base_expr_type_1;         \
      typedef typename ExprT2::base_expr_type base_expr_type_2;         \
      typedef typename Sacado::Promote<base_expr_type_1,                \
                                       base_expr_type_2>::type base_expr_type; \
                                                                        \
      static const int num_args1 = ExprT1::num_args;                    \
      static const int num_args2 = ExprT2::num_args;                    \
      static const int num_args = num_args1 + num_args2;                \
                                                                        \
      static const bool is_linear = LINEAR_2;                           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :                \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const {                                                \
        int sz1 = expr1.size(), sz2 = expr2.size();                     \
        return sz1 > sz2 ? sz1 : sz2;                                   \
      }                                                                 \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive() const {                                           \
        if (Arg < num_args1)                                            \
          return expr1.template isActive<Arg>();                        \
        else                                                            \
          return expr2.template isActive<Arg-num_args1>();              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive2(int j) const {                                     \
        if (j < num_args1)                                              \
          return expr1.isActive2(j);                                    \
        else                                                            \
          return expr2.isActive2(j);                                    \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool updateValue() const {                                        \
        return expr1.updateValue() && expr2.updateValue();              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {}                                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computePartials(const value_type& bar,                       \
                           value_type partials[]) const {               \
        if (num_args1 > 0)                                              \
          expr1.computePartials(LADJOINT, partials);                    \
        if (num_args2 > 0)                                              \
          expr2.computePartials(RADJOINT, partials+num_args1);          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void getTangents(int i, value_type dots[]) const {                \
        expr1.getTangents(i, dots);                                     \
        expr2.getTangents(i, dots+num_args1);                           \
      }                                                                 \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      const value_type& getTangent(int i) const {                       \
        if (Arg < num_args1)                                            \
          return expr1.template getTangent<Arg>(i);                     \
        else                                                            \
          return expr2.template getTangent<Arg-num_args1>(i);           \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isLinear() const {                                           \
        return LINEAR;                                                  \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess() && expr2.hasFastAccess();          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        return FASTACCESSDX;                                            \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type* getDx(int j) const {                            \
        if (j < num_args1)                                              \
          return expr1.getDx(j);                                        \
        else                                                            \
          return expr2.getDx(j-num_args1);                              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int numActiveArgs() const {                                       \
        return expr1.numActiveArgs() + expr2.numActiveArgs();           \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computeActivePartials(const value_type& bar,                 \
                                 value_type *partials) const {          \
        if (expr1.numActiveArgs() > 0)                                  \
          expr1.computePartials(LADJOINT, partials);                    \
        if (expr2.numActiveArgs() > 0)                                  \
          expr2.computePartials(RADJOINT, partials+expr2.numActiveArgs()); \
      }                                                                 \
    protected:                                                          \
                                                                        \
      typename ExprConstRef<ExprT1>::type expr1;                        \
      typename ExprConstRef<ExprT2>::type expr2;                        \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename ExprT1, typename T2>                             \
    class Expr< OP<ExprT1, ConstExpr<T2> > > {                          \
                                                                        \
    public:                                                             \
                                                                        \
      typedef ConstExpr<T2> ExprT2;                                     \
      typedef typename ExprT1::value_type value_type_1;                 \
      typedef typename ExprT2::value_type value_type_2;                 \
      typedef typename Sacado::Promote<value_type_1,                    \
                                       value_type_2>::type value_type;  \
                                                                        \
      typedef typename ExprT1::scalar_type scalar_type_1;               \
      typedef typename ExprT2::scalar_type scalar_type_2;               \
      typedef typename Sacado::Promote<scalar_type_1,                   \
                                       scalar_type_2>::type scalar_type; \
                                                                        \
      typedef typename ExprT1::base_expr_type base_expr_type_1;         \
      typedef typename ExprT2::base_expr_type base_expr_type_2;         \
      typedef typename Sacado::Promote<base_expr_type_1,                \
                                       base_expr_type_2>::type base_expr_type; \
                                                                        \
      static const int num_args = ExprT1::num_args;                     \
                                                                        \
      static const bool is_linear = CONST_LINEAR_2_2;                   \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :                \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr1.size();                                            \
      }                                                                 \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive() const {                                           \
        return expr1.template isActive<Arg>();                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive2(int j) const { return expr1.isActive2(j); }        \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool updateValue() const {                                        \
        return expr1.updateValue();                                     \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {}                                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computePartials(const value_type& bar,                       \
                           value_type partials[]) const {               \
        expr1.computePartials(LADJOINT, partials);                      \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void getTangents(int i, value_type dots[]) const {                \
        expr1.getTangents(i, dots);                                     \
      }                                                                 \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      const value_type& getTangent(int i) const {                       \
        return expr1.template getTangent<Arg>(i);                       \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isLinear() const {                                           \
        return CONST_LINEAR_2;                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        return CONST_DX_2;                                              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        return CONST_FASTACCESSDX_2;                                    \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type* getDx(int j) const {                            \
        return expr1.getDx(j);                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int numActiveArgs() const {                                       \
        return expr1.numActiveArgs();                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computeActivePartials(const value_type& bar,                 \
                                 value_type *partials) const {          \
        expr1.computePartials(LADJOINT, partials);                      \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      typename ExprConstRef<ExprT1>::type expr1;                        \
      typename ExprConstRef<ExprT2>::type expr2;                        \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename T1, typename ExprT2>                             \
    class Expr< OP<ConstExpr<T1>,ExprT2> > {                            \
                                                                        \
    public:                                                             \
                                                                        \
      typedef ConstExpr<T1> ExprT1;                                     \
      typedef typename ExprT1::value_type value_type_1;                 \
      typedef typename ExprT2::value_type value_type_2;                 \
      typedef typename Sacado::Promote<value_type_1,                    \
                                       value_type_2>::type value_type;  \
                                                                        \
      typedef typename ExprT1::scalar_type scalar_type_1;               \
      typedef typename ExprT2::scalar_type scalar_type_2;               \
      typedef typename Sacado::Promote<scalar_type_1,                   \
                                       scalar_type_2>::type scalar_type; \
                                                                        \
      typedef typename ExprT1::base_expr_type base_expr_type_1;         \
      typedef typename ExprT2::base_expr_type base_expr_type_2;         \
      typedef typename Sacado::Promote<base_expr_type_1,                \
                                       base_expr_type_2>::type base_expr_type; \
                                                                        \
      static const int num_args = ExprT2::num_args;                     \
                                                                        \
      static const bool is_linear = CONST_LINEAR_1_2;                   \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :                \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr2.size();                                            \
      }                                                                 \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive() const {                                           \
        return expr2.template isActive<Arg>();                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isActive2(int j) const { return expr2.isActive2(j); }        \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool updateValue() const {                                        \
        return expr2.updateValue();                                     \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {}                                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computePartials(const value_type& bar,                       \
                           value_type partials[]) const {               \
        expr2.computePartials(RADJOINT, partials);                      \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void getTangents(int i, value_type dots[]) const {                \
        expr2.getTangents(i, dots);                                     \
      }                                                                 \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      const value_type& getTangent(int i) const {                       \
        return expr2.template getTangent<Arg>(i);                       \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isLinear() const {                                           \
        return CONST_LINEAR_1;                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr2.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        return CONST_DX_1;                                              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        return CONST_FASTACCESSDX_1;                                    \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type* getDx(int j) const {                            \
        return expr2.getDx(j);                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int numActiveArgs() const {                                       \
        return expr2.numActiveArgs();                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computeActivePartials(const value_type& bar,                 \
                                 value_type *partials) const {          \
        expr2.computePartials(RADJOINT, partials);                      \
      }                                                                 \
    protected:                                                          \
                                                                        \
      typename ExprConstRef<ExprT1>::type expr1;                        \
      typename ExprConstRef<ExprT2>::type expr2;                        \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_EXPR(OP)                                  \
    OPNAME (const T1& expr1, const T2& expr2)                           \
    {                                                                   \
      typedef OP< T1, T2 > expr_t;                                      \
                                                                        \
      return Expr<expr_t>(expr1, expr2);                                \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< OP< Expr<T>, Expr<T> > >                                      \
    OPNAME (const Expr<T>& expr1, const Expr<T>& expr2)                 \
    {                                                                   \
      typedef OP< Expr<T>, Expr<T> > expr_t;                            \
                                                                        \
      return Expr<expr_t>(expr1, expr2);                                \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< OP< ConstExpr<typename Expr<T>::value_type>,                  \
              Expr<T> > >                                               \
    OPNAME (const typename Expr<T>::value_type& c,                      \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;           \
      typedef OP< ConstT, Expr<T> > expr_t;                             \
                                                                        \
      return Expr<expr_t>(ConstT(c), expr);                             \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< OP< Expr<T>,                                                  \
              ConstExpr<typename Expr<T>::value_type> > >               \
    OPNAME (const Expr<T>& expr,                                        \
            const typename Expr<T>::value_type& c)                      \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;           \
      typedef OP< Expr<T>, ConstT > expr_t;                             \
                                                                        \
      return Expr<expr_t>(expr, ConstT(c));                             \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_SCALAR_EXPR(OP)                                \
    OPNAME (const typename Expr<T>::scalar_type& c,                     \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;          \
      typedef OP< ConstT, Expr<T> > expr_t;                             \
                                                                        \
      return Expr<expr_t>(ConstT(c), expr);                             \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_SCALAR(OP)                                \
    OPNAME (const Expr<T>& expr,                                        \
            const typename Expr<T>::scalar_type& c)                     \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;          \
      typedef OP< Expr<T>, ConstT > expr_t;                             \
                                                                        \
      return Expr<expr_t>(expr, ConstT(c));                             \
    }                                                                   \
  }                                                                     \
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
                   expr2.size() == 0 && expr2.val() == value_type(1) ? bar : expr1.val() == value_type(0) ? value_type(0) : value_type(bar*std::pow(expr1.val(),expr2.val())*expr2.val()/expr1.val()),
                   expr1.val() == value_type(0) ? value_type(0) : value_type(bar*std::pow(expr1.val(),expr2.val())*std::log(expr1.val())),
                   false,
                   false,
                   false,
                   false,
                   false,
                   false,
                   expr1.val() == value_type(0) ? value_type(0) : value_type((expr2.dx(i)*std::log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0) : value_type((expr2.fastAccessDx(i)*std::log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0) : value_type(expr2.dx(i)*std::log(expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr2.val() == value_type(1) ? expr1.dx(i) : expr1.val() == value_type(0) ? value_type(0) : value_type(expr2.val()*expr1.dx(i)/expr1.val()*std::pow(expr1.val(),expr2.val())),
                   expr1.val() == value_type(0) ? value_type(0) : value_type(expr2.fastAccessDx(i)*std::log(expr1.val())*std::pow(expr1.val(),expr2.val())),
                   expr2.val() == value_type(1) ? expr1.fastAccessDx(i) : expr1.val() == value_type(0) ? value_type(0) : value_type(expr2.val()*expr1.fastAccessDx(i)/expr1.val()*std::pow(expr1.val(),expr2.val())))
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   expr1.val() >= expr2.val() ? expr1.val() : expr2.val(),
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
                   expr1.val() <= expr2.val() ? expr1.val() : expr2.val(),
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

#define FAD_RELOP_MACRO(OP)                                             \
namespace Sacado {                                                      \
  namespace ELRFad {                                                    \
    template <typename ExprT1, typename ExprT2>                         \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return expr1.val() OP expr2.val();                                \
    }                                                                   \
                                                                        \
    template <typename ExprT2>                                          \
    SACADO_INLINE_FUNCTION                                            \
    bool                                                                \
    operator OP (const typename Expr<ExprT2>::value_type& a,            \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return a OP expr2.val();                                          \
    }                                                                   \
                                                                        \
    template <typename ExprT1>                                          \
    SACADO_INLINE_FUNCTION                                            \
    bool                                                                \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const typename Expr<ExprT1>::value_type& b)            \
    {                                                                   \
      return expr1.val() OP b;                                          \
    }                                                                   \
  }                                                                     \
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
    SACADO_INLINE_FUNCTION
    bool operator ! (const Expr<ExprT>& expr)
    {
      return ! expr.val();
    }

  } // namespace ELRFad

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace ELRFad {

    template <typename ExprT>
    SACADO_INLINE_FUNCTION
    bool toBool(const Expr<ExprT>& x) {
      bool is_zero = (x.val() == 0.0);
      for (int i=0; i<x.size(); i++)
        is_zero = is_zero && (x.dx(i) == 0.0);
      return !is_zero;
    }

  } // namespace Fad

} // namespace Sacado

#define FAD_BOOL_MACRO(OP)                                              \
namespace Sacado {                                                      \
  namespace ELRFad {                                                    \
    template <typename ExprT1, typename ExprT2>                         \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return toBool(expr1) OP toBool(expr2);                            \
    }                                                                   \
                                                                        \
    template <typename ExprT2>                                          \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename Expr<ExprT2>::value_type& a,            \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return a OP toBool(expr2);                                        \
    }                                                                   \
                                                                        \
    template <typename ExprT1>                                          \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const typename Expr<ExprT1>::value_type& b)            \
    {                                                                   \
      return toBool(expr1) OP b;                                        \
    }                                                                   \
  }                                                                     \
}

FAD_BOOL_MACRO(&&)
FAD_BOOL_MACRO(||)

#undef FAD_BOOL_MACRO

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace ELRFad {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os << x.val() << " [";

      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_OPS_HPP
