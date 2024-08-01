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

#ifndef SACADO_ELRCACHEFAD_OPS_HPP
#define SACADO_ELRCACHEFAD_OPS_HPP

#include "Sacado_ELRCacheFad_Expression.hpp"
#include "Sacado_cmath.hpp"
#include "Sacado_mpl_disable_if.hpp"
#include <ostream>      // for std::ostream

namespace Sacado {
  namespace ELRCacheFad {

    //
    // UnaryPlusOp
    //

    template <typename ExprT>
    class UnaryPlusOp {};

    template <typename ExprT>
    class Expr< UnaryPlusOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;

      typedef typename ExprT::base_expr_type base_expr_type;

      static const int num_args = ExprT::num_args;

      static const bool is_linear = true;

      SACADO_INLINE_FUNCTION
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const { return expr.template isActive<Arg>(); }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr.computePartials(bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr.getTangents(i, dots); }

      template <int Arg>
      value_type getTangent(int i) const {
        return expr.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr.getDx(j);
      }

    protected:

      const ExprT& expr;
    };

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< UnaryPlusOp< Expr<T> > >
    operator+ (const Expr<T>& expr)
    {
      typedef UnaryPlusOp< Expr<T> > expr_t;

      return Expr<expr_t>(expr);
    }

    //
    // UnaryMinusOp
    //
    template <typename ExprT>
    class UnaryMinusOp {};

    template <typename ExprT>
    class Expr< UnaryMinusOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;

      typedef typename ExprT::base_expr_type base_expr_type;

      static const int num_args = ExprT::num_args;

      static const bool is_linear = true;

      SACADO_INLINE_FUNCTION
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const { return expr.template isActive<Arg>(); }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return -expr.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr.computePartials(-bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr.getTangents(i, dots); }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return -expr.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return -expr.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr.getDx(j);
      }

    protected:

      const ExprT& expr;
    };

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< UnaryMinusOp< Expr<T> > >
    operator- (const Expr<T>& expr)
    {
      typedef UnaryMinusOp< Expr<T> > expr_t;

      return Expr<expr_t>(expr);
    }

    //
    // AbsOp
    //

    template <typename ExprT>
    class AbsOp {};

    template <typename ExprT>
    class Expr< AbsOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;

      typedef typename ExprT::base_expr_type base_expr_type;

      static const int num_args = ExprT::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const { return expr.template isActive<Arg>(); }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr.cache();
        v = expr.val();
        v_pos = (v >= 0);
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return std::abs(v);
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (v_pos)
          expr.computePartials(bar, partials);
        else
          expr.computePartials(-bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr.getTangents(i, dots); }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        if (v_pos) return expr.dx(i);
        else return -expr.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        if (v_pos) return expr.fastAccessDx(i);
        else return -expr.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr.getDx(j);
      }

    protected:

      const ExprT& expr;
      mutable value_type v;
      mutable bool v_pos;
    };

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< AbsOp< Expr<T> > >
    abs (const Expr<T>& expr)
    {
      typedef AbsOp< Expr<T> > expr_t;

      return Expr<expr_t>(expr);
    }

    //
    // FAbsOp
    //

    template <typename ExprT>
    class FAbsOp {};

    template <typename ExprT>
    class Expr< FAbsOp<ExprT> > {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;

      typedef typename ExprT::base_expr_type base_expr_type;

      static const int num_args = ExprT::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const { return expr.template isActive<Arg>(); }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr.cache();
        v = expr.val();
        v_pos = (v >= 0);
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return std::fabs(v);
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (v_pos)
          expr.computePartials(bar, partials);
        else
          expr.computePartials(-bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr.getTangents(i, dots); }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        if (v_pos) return expr.dx(i);
        else return -expr.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        if (v_pos) return expr.fastAccessDx(i);
        else return -expr.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr.getDx(j);
      }

    protected:

      const ExprT& expr;
      mutable value_type v;
      mutable bool v_pos;
    };

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< FAbsOp< Expr<T> > >
    fabs (const Expr<T>& expr)
    {
      typedef FAbsOp< Expr<T> > expr_t;

      return Expr<expr_t>(expr);
    }

  }
}

#define FAD_UNARYOP_MACRO(OPNAME,OP,PARTIAL,VALUE)                      \
namespace Sacado {                                                      \
  namespace ELRCacheFad {                                               \
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
                                                                        \
      typedef typename ExprT::base_expr_type base_expr_type;            \
                                                                        \
      static const int num_args = ExprT::num_args;                      \
                                                                        \
      static const bool is_linear = false;                              \
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
      bool updateValue() const { return expr.updateValue(); }           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {                                              \
        expr.cache();                                                   \
        v = expr.val();                                                 \
        PARTIAL;                                                        \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void computePartials(const value_type& bar,                       \
                           value_type partials[]) const {               \
        expr.computePartials(bar*a, partials);                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void getTangents(int i, value_type dots[]) const {                \
        expr.getTangents(i, dots); }                                    \
                                                                        \
      template <int Arg>                                                \
      SACADO_INLINE_FUNCTION                                            \
      value_type getTangent(int i) const {                              \
        return expr.template getTangent<Arg>(i);                        \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isLinear() const {                                           \
        return false;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr.hasFastAccess();                                    \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        return expr.dx(i)*a;                                            \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        return expr.fastAccessDx(i)*a;                                  \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type* getDx(int j) const {                            \
        return expr.getDx(j);                                           \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ExprT& expr;                                                \
      mutable value_type v;                                             \
      mutable value_type a;                                             \
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

FAD_UNARYOP_MACRO(exp,
                  ExpOp,
                  a = std::exp(v),
                  a)
FAD_UNARYOP_MACRO(log,
                  LogOp,
                  a=scalar_type(1.0)/v,
                  std::log(v))
FAD_UNARYOP_MACRO(log10,
                  Log10Op,
                  a = scalar_type(1.0)/(std::log(scalar_type(10.0))*v),
                  std::log10(v))
FAD_UNARYOP_MACRO(sqrt,
                  SqrtOp,
                  a = scalar_type(1.0)/(scalar_type(2.0)*std::sqrt(v)),
                  std::sqrt(v))
FAD_UNARYOP_MACRO(safe_sqrt,
                  SafeSqrtOp,
                  a = (v == value_type(0.0) ? value_type(0.0) : value_type(scalar_type(1.0)/(scalar_type(2.0)*std::sqrt(v)))),
                  std::sqrt(v))
FAD_UNARYOP_MACRO(cos,
                  CosOp,
                  a = -std::sin(v),
                  std::cos(v))
FAD_UNARYOP_MACRO(sin,
                  SinOp,
                  a = std::cos(v),
                  std::sin(v))
FAD_UNARYOP_MACRO(tan,
                  TanOp,
                  a = scalar_type(1.0)+std::tan(v)*std::tan(v),
                  std::tan(v))
FAD_UNARYOP_MACRO(acos,
                  ACosOp,
                  a = scalar_type(-1.0)/std::sqrt(scalar_type(1.0)-v*v),
                  std::acos(v))
FAD_UNARYOP_MACRO(asin,
                  ASinOp,
                  a = scalar_type(1.0)/std::sqrt(scalar_type(1.0)-v*v),
                  std::asin(v))
FAD_UNARYOP_MACRO(atan,
                  ATanOp,
                  a = scalar_type(1.0)/(scalar_type(1.0)+v*v),
                  std::atan(v))
FAD_UNARYOP_MACRO(cosh,
                  CoshOp,
                  a = std::sinh(v),
                  std::cosh(v))
FAD_UNARYOP_MACRO(sinh,
                  SinhOp,
                  a = std::cosh(v),
                  std::sinh(v))
FAD_UNARYOP_MACRO(tanh,
                  TanhOp,
                  a = scalar_type(1.0)-std::tanh(v)*std::tanh(v),
                  std::tanh(v))
FAD_UNARYOP_MACRO(acosh,
                  ACoshOp,
                  a = scalar_type(1.0)/std::sqrt((v-scalar_type(1.0))*(v+scalar_type(1.0))),
                  std::acosh(v))
FAD_UNARYOP_MACRO(asinh,
                  ASinhOp,
                  a = scalar_type(1.0)/std::sqrt(scalar_type(1.0)+v*v),
                  std::asinh(v))
FAD_UNARYOP_MACRO(atanh,
                  ATanhOp,
                  a = scalar_type(1.0)/(scalar_type(1.0)-v*v),
                  std::atanh(v))
FAD_UNARYOP_MACRO(cbrt,
                  CbrtOp,
                  a = scalar_type(1.0)/(scalar_type(3.0)*std::cbrt(v*v)),
                  std::cbrt(v))

#undef FAD_UNARYOP_MACRO

//
// Binary operators
//
namespace Sacado {
  namespace ELRCacheFad {

    //
    // AdditionOp
    //

    template <typename ExprT1, typename ExprT2>
    class AdditionOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< AdditionOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = ExprT1::is_linear && ExprT2::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val()+expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0)
          expr1.computePartials(bar, partials);
        if (num_args2 > 0)
          expr2.computePartials(bar, partials+num_args1);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr1.isLinear() && expr2.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i) + expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i) + expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< AdditionOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = ExprT1::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val() + expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr1.computePartials(bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr1.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;

    };

    template <typename T1, typename ExprT2>
    class Expr< AdditionOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = ExprT2::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val() + expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr2.computePartials(bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr2.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;

    };

    //
    // SubtractionOp
    //

    template <typename ExprT1, typename ExprT2>
    class SubtractionOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< SubtractionOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = ExprT1::is_linear && ExprT2::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val()-expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0)
          expr1.computePartials(bar, partials);
        if (num_args2 > 0)
          expr2.computePartials(-bar, partials+num_args1);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr1.isLinear() && expr2.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i) - expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i) - expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< SubtractionOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = ExprT1::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val() - expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr1.computePartials(bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr1.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;

    };

    template <typename T1, typename ExprT2>
    class Expr< SubtractionOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = ExprT2::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val() - expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr2.computePartials(-bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr2.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return -expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return -expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;

    };

    //
    // MultiplicationOp
    //

    template <typename ExprT1, typename ExprT2>
    class MultiplicationOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< MultiplicationOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
        v1 = expr1.val();
        v2 = expr2.val();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v1*v2;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0)
          expr1.computePartials(bar*v2, partials);
        if (num_args2 > 0)
          expr2.computePartials(bar*v1, partials+num_args1);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
          return v1*expr2.dx(i) + expr1.dx(i)*v2;
        else if (expr1.size() > 0)
          return expr1.dx(i)*v2;
        else
          return v1*expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return v1*expr2.fastAccessDx(i) + expr1.fastAccessDx(i)*v2;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type_1 v1;
      mutable value_type_2 v2;

    };

    template <typename ExprT1, typename T2>
    class Expr< MultiplicationOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = ExprT1::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val()*expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr1.computePartials(bar*expr2.val(), partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr1.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i)*expr2.val();
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*expr2.val();
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;

    };

    template <typename T1, typename ExprT2>
    class Expr< MultiplicationOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = ExprT2::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return expr1.val()*expr2.val();
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr2.computePartials(bar*expr1.val(), partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr2.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.val()*expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.val()*expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;

    };

    //
    // DivisionOp
    //

    template <typename ExprT1, typename ExprT2>
    class DivisionOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< DivisionOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        a = scalar_type(1.0)/v2;
        v = v1*a;
        b = -v/v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0)
          expr1.computePartials(bar*a, partials);
        if (num_args2 > 0)
          expr2.computePartials(bar*b, partials+num_args1);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
          return expr1.dx(i)*a + expr2.dx(i)*b;
        else if (expr1.size() > 0)
          return expr1.dx(i)*a;
        else
          return expr1.val()*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return  expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type a;
      mutable value_type b;

    };

    template <typename ExprT1, typename T2>
    class Expr< DivisionOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = ExprT1::is_linear;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        const value_type_1 v1 = expr1.val();
        a = scalar_type(1.0)/expr2.val();
        v = v1*a;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr1.computePartials(bar*a, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return expr1.isLinear();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i)*a;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable value_type a;

    };

    template <typename T1, typename ExprT2>
    class Expr< DivisionOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
        const value_type_2 v2 = expr2.val();
        v = expr1.val()/v2;
        b = -v/v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr2.computePartials(bar*b, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr2.dx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type b;

    };

    //
    // Atan2Op
    //

    template <typename ExprT1, typename ExprT2>
    class Atan2Op {};

    template <typename ExprT1, typename ExprT2>
    class Expr< Atan2Op<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        a = scalar_type(1.0)/(v1*v1 + v2*v2);
        b = -v1*a;
        a = v2*a;
        v = std::atan2(v1,v2);
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0)
          expr1.computePartials(bar*a, partials);
        if (num_args2 > 0)
          expr2.computePartials(bar*b, partials+num_args1);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
          return expr1.dx(i)*a + expr2.dx(i)*b;
        else if (expr1.size() > 0)
          return expr1.dx(i)*a;
        else
          return expr1.val()*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type a;
      mutable value_type b;

    };

    template <typename ExprT1, typename T2>
    class Expr< Atan2Op<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        a = v2/(v1*v1 + v2*v2);
        v = std::atan2(v1,v2);
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr1.computePartials(bar*a, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i)*a;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable value_type a;

    };

    template <typename T1, typename ExprT2>
    class Expr< Atan2Op< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        b = -v1/(v1*v1 + v2*v2);
        v = std::atan2(v1,v2);
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr2.computePartials(bar*b, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr2.dx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type b;

    };

    //
    // PowerOp
    //

    template <typename ExprT1, typename ExprT2>
    class PowerOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< PowerOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        v = std::pow(v1,v2);
        if (expr2.size() == 0 && v2 == scalar_type(1.0)) {
          a = scalar_type(1.0);
          b = scalar_type(0.0);
        }
        else if (v1 == scalar_type(0.0)) {
          a = scalar_type(0.0);
          b = scalar_type(0.0);
        }
        else {
          a = v*v2/v1;
          b = v*std::log(v1);
        }
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0)
          expr1.computePartials(bar*a, partials);
        if (num_args2 > 0)
          expr2.computePartials(bar*b, partials+num_args1);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
          return expr1.dx(i)*a + expr2.dx(i)*b;
        else if (expr1.size() > 0)
          return expr1.dx(i)*a;
        else
          return expr1.val()*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a + expr2.fastAccessDx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type a;
      mutable value_type b;

    };

    template <typename ExprT1, typename T2>
    class Expr< PowerOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        v = std::pow(v1,v2);
        if (v2 == scalar_type(1.0)) {
          a = scalar_type(1.0);
        }
        else if (v1 == scalar_type(0.0)) {
          a = scalar_type(0.0);
        }
        else {
          a = v*v2/v1;
        }
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr1.computePartials(bar*a, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i)*a;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*a;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable value_type a;

    };

    template <typename T1, typename ExprT2>
    class Expr< PowerOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        v = std::pow(v1,v2);
        if (v1 == scalar_type(0.0)) {
          b = scalar_type(0.0);
        }
        else {
          b = v*std::log(v1);
        }
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        expr2.computePartials(bar*b, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr2.dx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr2.fastAccessDx(i)*b;
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable value_type b;

    };

    //
    // MaxOp
    //

    template <typename ExprT1, typename ExprT2>
    class MaxOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< MaxOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        max_v1 = (v1 >= v2);
        v = max_v1 ? v1 : v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0) {
          if (max_v1)
            expr1.computePartials(bar, partials);
          else
            expr1.computePartials(value_type(0.0), partials);
        }
        if (num_args2 > 0) {
          if (max_v1)
            expr2.computePartials(value_type(0.0), partials+num_args1);
          else
            expr2.computePartials(bar, partials+num_args1);
        }
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return max_v1 ? expr1.dx(i) : expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return max_v1 ? expr1.fastAccessDx(i) : expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool max_v1;

    };

    template <typename ExprT1, typename T2>
    class Expr< MaxOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        max_v1 = (v1 >= v2);
        v = max_v1 ? v1 : v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (max_v1)
          expr1.computePartials(bar, partials);
        else
          expr1.computePartials(value_type(0.0), partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return max_v1 ? expr1.dx(i) : value_type(0.0);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return max_v1 ? expr1.fastAccessDx(i) : value_type(0.0);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable bool max_v1;

    };

    template <typename T1, typename ExprT2>
    class Expr< MaxOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        max_v1 = (v1 >= v2);
        v = max_v1 ? v1 : v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (max_v1)
          expr2.computePartials(value_type(0.0), partials);
        else
          expr2.computePartials(bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return max_v1 ? value_type(0.0) : expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return max_v1 ? value_type(0.0) : expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool max_v1;

    };

    //
    // MinOp
    //

    template <typename ExprT1, typename ExprT2>
    class MinOp {};

    template <typename ExprT1, typename ExprT2>
    class Expr< MinOp<ExprT1,ExprT2> > {

    public:

      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args1 = ExprT1::num_args;
      static const int num_args2 = ExprT2::num_args;
      static const int num_args = num_args1 + num_args2;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        if (Arg < num_args1)
          return expr1.template isActive<Arg>();
        else
          return expr2.template isActive<Arg-num_args1>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        min_v1 = (v1 <= v2);
        v = min_v1 ? v1 : v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (num_args1 > 0) {
          if (min_v1)
            expr1.computePartials(bar, partials);
          else
            expr1.computePartials(value_type(0.0), partials);
        }
        if (num_args2 > 0) {
          if (min_v1)
            expr2.computePartials(value_type(0.0), partials+num_args1);
          else
            expr2.computePartials(bar, partials+num_args1);
        }
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
        expr2.getTangents(i, dots+num_args1);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        if (Arg < num_args1)
          return expr1.template getTangent<Arg>(i);
        else
          return expr2.template getTangent<Arg-num_args1>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return min_v1 ? expr1.dx(i) : expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return min_v1 ? expr1.fastAccessDx(i) : expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        if (j < num_args1)
          return expr1.getDx(j);
        else
          return expr2.getDx(j-num_args1);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool min_v1;

    };

    template <typename ExprT1, typename T2>
    class Expr< MinOp<ExprT1, ConstExpr<T2> > > {

    public:

      typedef ConstExpr<T2> ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT1::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr1.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr1.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        min_v1 = (v1 <= v2);
        v = min_v1 ? v1 : v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (min_v1)
          expr1.computePartials(bar, partials);
        else
          expr1.computePartials(value_type(0.0), partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr1.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr1.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return min_v1 ? expr1.dx(i) : value_type(0.0);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return min_v1 ? expr1.fastAccessDx(i) : value_type(0.0);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr1.getDx(j);
      }

    protected:

      const ExprT1& expr1;
      ExprT2 expr2;
      mutable value_type v;
      mutable bool min_v1;

    };

    template <typename T1, typename ExprT2>
    class Expr< MinOp< ConstExpr<T1>,ExprT2> > {

    public:

      typedef ConstExpr<T1> ExprT1;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;
      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename ExprT1::base_expr_type base_expr_type_1;
      typedef typename ExprT2::base_expr_type base_expr_type_2;
      typedef typename Sacado::Promote<base_expr_type_1,
                                       base_expr_type_2>::type base_expr_type;

      static const int num_args = ExprT2::num_args;

      static const bool is_linear = false;

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const {
        return expr2.template isActive<Arg>();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {
        expr2.cache();
        const value_type_1 v1 = expr1.val();
        const value_type_2 v2 = expr2.val();
        min_v1 = (v1 <= v2);
        v = min_v1 ? v1 : v2;
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        return v;
      }

      SACADO_INLINE_FUNCTION
      void computePartials(const value_type& bar,
                           value_type partials[]) const {
        if (min_v1)
          expr2.computePartials(value_type(0.0), partials);
        else
          expr2.computePartials(bar, partials);
      }

      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        expr2.getTangents(i, dots);
      }

      template <int Arg>
      SACADO_INLINE_FUNCTION
      value_type getTangent(int i) const {
        return expr2.template getTangent<Arg>(i);
      }

      SACADO_INLINE_FUNCTION
      bool isLinear() const {
        return false;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return min_v1 ? value_type(0.0) : expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return min_v1 ? value_type(0.0) : expr2.fastAccessDx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const {
        return expr2.getDx(j);
      }

    protected:

      ExprT1 expr1;
      const ExprT2& expr2;
      mutable value_type v;
      mutable bool min_v1;

    };

  }

}

#define FAD_BINARYOP_MACRO(OPNAME,OP)                                   \
namespace Sacado {                                                      \
  namespace ELRCacheFad {                                               \
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


FAD_BINARYOP_MACRO(operator+, AdditionOp)
FAD_BINARYOP_MACRO(operator-, SubtractionOp)
FAD_BINARYOP_MACRO(operator*, MultiplicationOp)
FAD_BINARYOP_MACRO(operator/, DivisionOp)
FAD_BINARYOP_MACRO(atan2, Atan2Op)
FAD_BINARYOP_MACRO(pow, PowerOp)
FAD_BINARYOP_MACRO(max, MaxOp)
FAD_BINARYOP_MACRO(min, MinOp)

#undef FAD_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define FAD_RELOP_MACRO(OP)                                             \
namespace Sacado {                                                      \
  namespace ELRCacheFad {                                               \
    template <typename ExprT1, typename ExprT2>                         \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      expr1.cache();                                                    \
      expr2.cache();                                                    \
      return expr1.val() OP expr2.val();                                \
    }                                                                   \
                                                                        \
    template <typename ExprT2>                                            \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename Expr<ExprT2>::value_type& a,            \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      expr2.cache();                                                    \
      return a OP expr2.val();                                          \
    }                                                                   \
                                                                        \
    template <typename ExprT1>                                          \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const typename Expr<ExprT1>::value_type& b)            \
    {                                                                   \
      expr1.cache();                                                    \
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

  namespace ELRCacheFad {

    template <typename ExprT>
    SACADO_INLINE_FUNCTION
    bool operator ! (const Expr<ExprT>& expr)
    {
      expr.cache();
      return ! expr.val();
    }

  } // namespace ELRCacheFad

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace ELRCacheFad {

    template <typename ExprT>
    SACADO_INLINE_FUNCTION
    bool toBool(const Expr<ExprT>& x) {
      x.cache();
      bool is_zero = (x.val() == 0.0);
      for (int i=0; i<x.size(); i++)
        is_zero = is_zero && (x.dx(i) == 0.0);
      return !is_zero;
    }

  } // namespace Fad

} // namespace Sacado

#define FAD_BOOL_MACRO(OP)                                              \
namespace Sacado {                                                      \
  namespace ELRCacheFad {                                               \
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

  namespace ELRCacheFad {

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

#endif // SACADO_CACHEFAD_OPS_HPP
