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

#ifndef SACADO_FAD_OPS_HPP
#define SACADO_FAD_OPS_HPP

#include "Sacado_Fad_Expression.hpp"
#include "Sacado_Fad_Ops_Fwd.hpp"
#include "Sacado_cmath.hpp"
#include <ostream>      // for std::ostream

#define FAD_UNARYOP_MACRO(OPNAME,OP,USING,VALUE,DX,FASTACCESSDX)        \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
                                                                        \
    template <typename ExprT>                                           \
    class OP {};                                                        \
                                                                        \
    template <typename ExprT>                                           \
    struct ExprSpec< OP<ExprT> > {                                      \
      typedef typename ExprSpec<ExprT>::type type;                      \
    };                                                                  \
                                                                        \
    template <typename ExprT>                                           \
    class Expr< OP<ExprT>,ExprSpecDefault > {                           \
    public:                                                             \
                                                                        \
      typedef typename ExprT::value_type value_type;                    \
      typedef typename ExprT::scalar_type scalar_type;                  \
      typedef typename ExprT::base_expr_type base_expr_type;            \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}               \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const { return expr.size(); }                          \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const { return expr.hasFastAccess(); }       \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isPassive() const { return expr.isPassive();}                \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool updateValue() const { return expr.updateValue(); }           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {}                                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        USING                                                           \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type dx(int i) const {                                      \
        USING                                                           \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type fastAccessDx(int i) const {                            \
        USING                                                           \
        return FASTACCESSDX;                                            \
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
                                                                        \
}

FAD_UNARYOP_MACRO(operator+,
                  UnaryPlusOp,
                  ;,
                  expr.val(),
                  expr.dx(i),
                  expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(operator-,
                  UnaryMinusOp,
                  ;,
                  -expr.val(),
                  -expr.dx(i),
                  -expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(exp,
                  ExpOp,
                  using std::exp;,
                  exp(expr.val()),
                  exp(expr.val())*expr.dx(i),
                  exp(expr.val())*expr.fastAccessDx(i))
FAD_UNARYOP_MACRO(log,
                  LogOp,
                  using std::log;,
                  log(expr.val()),
                  expr.dx(i)/expr.val(),
                  expr.fastAccessDx(i)/expr.val())
FAD_UNARYOP_MACRO(log10,
                  Log10Op,
                  using std::log10; using std::log;,
                  log10(expr.val()),
                  expr.dx(i)/( log(value_type(10))*expr.val()),
                  expr.fastAccessDx(i) / ( log(value_type(10))*expr.val()))
FAD_UNARYOP_MACRO(sqrt,
                  SqrtOp,
                  using std::sqrt;,
                  sqrt(expr.val()),
                  expr.dx(i)/(value_type(2)* sqrt(expr.val())),
                  expr.fastAccessDx(i)/(value_type(2)* sqrt(expr.val())))
FAD_UNARYOP_MACRO(cos,
                  CosOp,
                  using std::cos; using std::sin;,
                  cos(expr.val()),
                  -expr.dx(i)* sin(expr.val()),
                  -expr.fastAccessDx(i)* sin(expr.val()))
FAD_UNARYOP_MACRO(sin,
                  SinOp,
                  using std::cos; using std::sin;,
                  sin(expr.val()),
                  expr.dx(i)* cos(expr.val()),
                  expr.fastAccessDx(i)* cos(expr.val()))
FAD_UNARYOP_MACRO(tan,
                  TanOp,
                  using std::tan;,
                  std::tan(expr.val()),
                  expr.dx(i)*
                    (value_type(1)+ tan(expr.val())* tan(expr.val())),
                  expr.fastAccessDx(i)*
                    (value_type(1)+ tan(expr.val())* tan(expr.val())))
FAD_UNARYOP_MACRO(acos,
                  ACosOp,
                  using std::acos; using std::sqrt;,
                  acos(expr.val()),
                  -expr.dx(i)/ sqrt(value_type(1)-expr.val()*expr.val()),
                  -expr.fastAccessDx(i) /
                    sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(asin,
                  ASinOp,
                  using std::asin; using std::sqrt;,
                  asin(expr.val()),
                  expr.dx(i)/ sqrt(value_type(1)-expr.val()*expr.val()),
                  expr.fastAccessDx(i) /
                    sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atan,
                  ATanOp,
                  using std::atan;,
                  atan(expr.val()),
                  expr.dx(i)/(value_type(1)+expr.val()*expr.val()),
                  expr.fastAccessDx(i)/(value_type(1)+expr.val()*expr.val()))
FAD_UNARYOP_MACRO(cosh,
                  CoshOp,
                  using std::cosh; using std::sinh;,
                  cosh(expr.val()),
                  expr.dx(i)* sinh(expr.val()),
                  expr.fastAccessDx(i)* sinh(expr.val()))
FAD_UNARYOP_MACRO(sinh,
                  SinhOp,
                  using std::cosh; using std::sinh;,
                  sinh(expr.val()),
                  expr.dx(i)* cosh(expr.val()),
                  expr.fastAccessDx(i)* cosh(expr.val()))
FAD_UNARYOP_MACRO(tanh,
                  TanhOp,
                  using std::tanh;,
                  tanh(expr.val()),
                  expr.dx(i)*(value_type(1)-tanh(expr.val())*tanh(expr.val())),
                  expr.fastAccessDx(i)*(value_type(1)-tanh(expr.val())*tanh(expr.val())))
FAD_UNARYOP_MACRO(acosh,
                  ACoshOp,
                  using std::acosh; using std::sqrt;,
                  acosh(expr.val()),
                  expr.dx(i)/ sqrt((expr.val()-value_type(1)) *
                                       (expr.val()+value_type(1))),
                  expr.fastAccessDx(i)/ sqrt((expr.val()-value_type(1)) *
                                                 (expr.val()+value_type(1))))
FAD_UNARYOP_MACRO(asinh,
                  ASinhOp,
                  using std::asinh; using std::sqrt;,
                  asinh(expr.val()),
                  expr.dx(i)/ sqrt(value_type(1)+expr.val()*expr.val()),
                  expr.fastAccessDx(i)/ sqrt(value_type(1)+
                                                 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atanh,
                  ATanhOp,
                  using std::atanh;,
                  atanh(expr.val()),
                  expr.dx(i)/(value_type(1)-expr.val()*expr.val()),
                  expr.fastAccessDx(i)/(value_type(1)-
                                                 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(abs,
                  AbsOp,
                  using std::abs; using Sacado::if_then_else;,
                  abs(expr.val()),
                  if_then_else( expr.val() >= 0, expr.dx(i), value_type(-expr.dx(i)) ),
                  if_then_else( expr.val() >= 0, expr.fastAccessDx(i), value_type(-expr.fastAccessDx(i)) ) )
FAD_UNARYOP_MACRO(fabs,
                  FAbsOp,
                  using std::fabs; using Sacado::if_then_else;,
                  fabs(expr.val()),
                  if_then_else( expr.val() >= 0, expr.dx(i), value_type(-expr.dx(i)) ),
                  if_then_else( expr.val() >= 0, expr.fastAccessDx(i), value_type(-expr.fastAccessDx(i)) ) )
FAD_UNARYOP_MACRO(cbrt,
                  CbrtOp,
                  using std::cbrt;,
                  cbrt(expr.val()),
                  expr.dx(i)/(value_type(3)*cbrt(expr.val()*expr.val())),
                  expr.fastAccessDx(i)/(value_type(3)*cbrt(expr.val()*expr.val())))

#undef FAD_UNARYOP_MACRO

// Special handling for safe_sqrt() to provide specializations of SafeSqrtOp for
// "simd" value types that use if_then_else(). The only reason for not using
// if_then_else() always is to avoid evaluating the derivative if the value is
// zero to avoid throwing FPEs.
namespace Sacado {
  namespace Fad {

    template <typename ExprT, bool is_simd>
    class SafeSqrtOp {};

    template <typename ExprT>
    struct ExprSpec< SafeSqrtOp<ExprT> > {
      typedef typename ExprSpec<ExprT>::type type;
    };

    //
    // Implementation for simd type using if_then_else()
    //
    template <typename ExprT>
    class Expr< SafeSqrtOp<ExprT,true>,ExprSpecDefault > {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      SACADO_INLINE_FUNCTION
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const { return expr.hasFastAccess(); }

      SACADO_INLINE_FUNCTION
      bool isPassive() const { return expr.isPassive();}

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      value_type val() const {
        using std::sqrt;
        return sqrt(expr.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::sqrt; using Sacado::if_then_else;
        return if_then_else(
          expr.val() == value_type(0.0), value_type(0.0),
          value_type(expr.dx(i)/(value_type(2)*sqrt(expr.val()))));
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::sqrt; using Sacado::if_then_else;
        return if_then_else(
          expr.val() == value_type(0.0), value_type(0.0),
          value_type(expr.fastAccessDx(i)/(value_type(2)*sqrt(expr.val()))));
      }

    protected:

      const ExprT& expr;
    };

    //
    // Specialization for scalar types using ternary operator
    //
    template <typename ExprT>
    class Expr< SafeSqrtOp<ExprT,false>,ExprSpecDefault > {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      SACADO_INLINE_FUNCTION
      explicit Expr(const ExprT& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const { return expr.hasFastAccess(); }

      SACADO_INLINE_FUNCTION
      bool isPassive() const { return expr.isPassive();}

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      value_type val() const {
        using std::sqrt;
        return sqrt(expr.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::sqrt;
        return expr.val() == value_type(0.0) ? value_type(0.0) :
          value_type(expr.dx(i)/(value_type(2)*sqrt(expr.val())));
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::sqrt;
        return expr.val() == value_type(0.0) ? value_type(0.0) :
          value_type(expr.fastAccessDx(i)/(value_type(2)*sqrt(expr.val())));
      }

    protected:

      const ExprT& expr;
    };

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< SafeSqrtOp< Expr<T> > >
    safe_sqrt (const Expr<T>& expr)
    {
      typedef SafeSqrtOp< Expr<T> > expr_t;

      return Expr<expr_t>(expr);
    }
  }

}

#define FAD_BINARYOP_MACRO(OPNAME,OP,USING,VALUE,DX,FASTACCESSDX,VAL_CONST_DX_1,VAL_CONST_DX_2,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
                                                                        \
    template <typename ExprT1, typename ExprT2>                         \
    class OP {};                                                        \
                                                                        \
    template <typename ExprT1, typename ExprT2>                         \
    struct ExprSpec< OP< ExprT1, ExprT2 > > {                           \
      typedef typename ExprSpec<ExprT1>::type type;                     \
    };                                                                  \
                                                                        \
    template <typename ExprT1, typename ExprT2>                         \
    class Expr< OP< ExprT1, ExprT2 >,ExprSpecDefault > {                \
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
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess() && expr2.hasFastAccess();          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isPassive() const {                                          \
        return expr1.isPassive() && expr2.isPassive();                  \
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
      const value_type val() const {                                    \
        USING                                                           \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        USING                                                           \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        USING                                                           \
        return FASTACCESSDX;                                            \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ExprT1& expr1;                                              \
      const ExprT2& expr2;                                              \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename ExprT1, typename T2>                             \
    struct ExprSpec< OP< ExprT1, ConstExpr<T2> > > {                    \
      typedef typename ExprSpec<ExprT1>::type type;                     \
    };                                                                  \
                                                                        \
    template <typename ExprT1, typename T2>                             \
    class Expr< OP< ExprT1, ConstExpr<T2> >,ExprSpecDefault > {         \
                                                                        \
    public:                                                             \
                                                                        \
      typedef ConstExpr<T2> ConstT;                                     \
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
      SACADO_INLINE_FUNCTION                                            \
      Expr(const ExprT1& expr1_, const ConstT& c_) :                    \
        expr1(expr1_), c(c_) {}                                         \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr1.size();                                            \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isPassive() const {                                          \
        return expr1.isPassive();                                       \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool updateValue() const { return expr1.updateValue(); }          \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {}                                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type val() const {                                    \
        USING                                                           \
        return VAL_CONST_DX_2;                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        USING                                                           \
        return CONST_DX_2;                                              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        USING                                                           \
        return CONST_FASTACCESSDX_2;                                    \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ExprT1& expr1;                                              \
      ConstT c;                                                         \
    };                                                                  \
                                                                        \
    template <typename T1, typename ExprT2>                             \
    struct ExprSpec< OP< ConstExpr<T1>, ExprT2 > > {                    \
      typedef typename ExprSpec<ExprT2>::type type;                     \
    };                                                                  \
                                                                        \
    template <typename T1, typename ExprT2>                             \
    class Expr< OP< ConstExpr<T1>, ExprT2 >,ExprSpecDefault > {         \
                                                                        \
    public:                                                             \
                                                                        \
      typedef ConstExpr<T1> ConstT;                                     \
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
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      Expr(const ConstT& c_, const ExprT2& expr2_) :                    \
        c(c_), expr2(expr2_) {}                                         \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr2.size();                                            \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr2.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool isPassive() const {                                          \
        return expr2.isPassive();                                       \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool updateValue() const { return expr2.updateValue(); }          \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      void cache() const {}                                             \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type val() const {                                    \
        USING                                                           \
        return VAL_CONST_DX_1;                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type dx(int i) const {                                \
        USING                                                           \
        return CONST_DX_1;                                              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      const value_type fastAccessDx(int i) const {                      \
        USING                                                           \
        return CONST_FASTACCESSDX_1;                                    \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      ConstT c;                                                         \
      const ExprT2& expr2;                                              \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    typename mpl::enable_if_c<                                          \
       ExprLevel< Expr<T1> >::value == ExprLevel< Expr<T2> >::value,    \
       Expr< OP< Expr<T1>, Expr<T2> > >                                 \
     >::type                                                            \
    /*SACADO_FAD_OP_ENABLE_EXPR_EXPR(OP)*/                              \
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)               \
    {                                                                   \
      typedef OP< Expr<T1>, Expr<T2> > expr_t;                          \
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
                                                                        \
}


FAD_BINARYOP_MACRO(operator+,
                   AdditionOp,
                   ;,
                   expr1.val() + expr2.val(),
                   expr1.dx(i) + expr2.dx(i),
                   expr1.fastAccessDx(i) + expr2.fastAccessDx(i),
                   c.val() + expr2.val(),
                   expr1.val() + c.val(),
                   expr2.dx(i),
                   expr1.dx(i),
                   expr2.fastAccessDx(i),
                   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator-,
                   SubtractionOp,
                   ;,
                   expr1.val() - expr2.val(),
                   expr1.dx(i) - expr2.dx(i),
                   expr1.fastAccessDx(i) - expr2.fastAccessDx(i),
                   c.val() - expr2.val(),
                   expr1.val() - c.val(),
                   -expr2.dx(i),
                   expr1.dx(i),
                   -expr2.fastAccessDx(i),
                   expr1.fastAccessDx(i))
// FAD_BINARYOP_MACRO(operator*,
//                 MultiplicationOp,
//                 ;,
//                 expr1.val() * expr2.val(),
//                 expr1.val()*expr2.dx(i) + expr1.dx(i)*expr2.val(),
//                 expr1.val()*expr2.fastAccessDx(i) +
//                   expr1.fastAccessDx(i)*expr2.val(),
//                 c.val() * expr2.val(),
//                 expr1.val() * c.val(),
//                 c.val()*expr2.dx(i),
//                 expr1.dx(i)*c.val(),
//                 c.val()*expr2.fastAccessDx(i),
//                 expr1.fastAccessDx(i)*c.val())
FAD_BINARYOP_MACRO(operator/,
                   DivisionOp,
                   ;,
                   expr1.val() / expr2.val(),
                   (expr1.dx(i)*expr2.val() - expr2.dx(i)*expr1.val()) /
                     (expr2.val()*expr2.val()),
                   (expr1.fastAccessDx(i)*expr2.val() -
                      expr2.fastAccessDx(i)*expr1.val()) /
                      (expr2.val()*expr2.val()),
                   c.val() / expr2.val(),
                   expr1.val() / c.val(),
                   -expr2.dx(i)*c.val() / (expr2.val()*expr2.val()),
                   expr1.dx(i)/c.val(),
                   -expr2.fastAccessDx(i)*c.val() / (expr2.val()*expr2.val()),
                   expr1.fastAccessDx(i)/c.val())
FAD_BINARYOP_MACRO(atan2,
                   Atan2Op,
                   using std::atan2;,
                   atan2(expr1.val(), expr2.val()),
                   (expr2.val()*expr1.dx(i) - expr1.val()*expr2.dx(i))/
                        (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   (expr2.val()*expr1.fastAccessDx(i) - expr1.val()*expr2.fastAccessDx(i))/
                        (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   atan2(c.val(), expr2.val()),
                   atan2(expr1.val(), c.val()),
                   (-c.val()*expr2.dx(i)) / (c.val()*c.val() + expr2.val()*expr2.val()),
                   (c.val()*expr1.dx(i))/ (expr1.val()*expr1.val() + c.val()*c.val()),
                   (-c.val()*expr2.fastAccessDx(i))/ (c.val()*c.val() + expr2.val()*expr2.val()),
                   (c.val()*expr1.fastAccessDx(i))/ (expr1.val()*expr1.val() + c.val()*c.val()))
// FAD_BINARYOP_MACRO(pow,
//                    PowerOp,
//                    using std::pow; using std::log; using Sacado::if_then_else;,
//                    pow(expr1.val(), expr2.val()),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type((expr2.dx(i)*log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*pow(expr1.val(),expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type((expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val())) ),
//                    pow(c.val(), expr2.val()),
//                    pow(expr1.val(), c.val()),
//                    if_then_else( c.val() == value_type(0.0), value_type(0.0), value_type(expr2.dx(i)*log(c.val())*pow(c.val(),expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type(c.val()*expr1.dx(i)/expr1.val()*pow(expr1.val(),c.val())) ),
//                    if_then_else( c.val() == value_type(0.0), value_type(0.0), value_type(expr2.fastAccessDx(i)*log(c.val())*pow(c.val(),expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type(c.val()*expr1.fastAccessDx(i)/expr1.val()*pow(expr1.val(),c.val()))) )
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   using Sacado::if_then_else;,
                   if_then_else( expr1.val() >= expr2.val(),  expr1.val(), expr2.val() ),
                   if_then_else( expr1.val() >= expr2.val(), expr1.dx(i), expr2.dx(i) ),
                   if_then_else( expr1.val() >= expr2.val(), expr1.fastAccessDx(i), expr2.fastAccessDx(i) ),
                   if_then_else( c.val() >= expr2.val(), value_type(c.val()),  expr2.val() ),
                   if_then_else( expr1.val() >= c.val(), expr1.val(), value_type(c.val()) ),
                   if_then_else( c.val() >= expr2.val(), value_type(0.0),  expr2.dx(i) ),
                   if_then_else( expr1.val() >= c.val(), expr1.dx(i), value_type(0.0) ),
                   if_then_else( c.val() >= expr2.val(), value_type(0.0), expr2.fastAccessDx(i) ),
                   if_then_else( expr1.val() >= c.val(), expr1.fastAccessDx(i), value_type(0.0) ) )
FAD_BINARYOP_MACRO(min,
                   MinOp,
                   using Sacado::if_then_else;,
                   if_then_else( expr1.val() <= expr2.val(), expr1.val(), expr2.val() ),
                   if_then_else( expr1.val() <= expr2.val(), expr1.dx(i), expr2.dx(i) ),
                   if_then_else( expr1.val() <= expr2.val(), expr1.fastAccessDx(i), expr2.fastAccessDx(i) ),
                   if_then_else( c.val() <= expr2.val(), value_type(c.val()), expr2.val() ),
                   if_then_else( expr1.val() <= c.val(), expr1.val(), value_type(c.val()) ),
                   if_then_else( c.val() <= expr2.val(), value_type(0), expr2.dx(i) ),
                   if_then_else( expr1.val() <= c.val(), expr1.dx(i), value_type(0) ),
                   if_then_else( c.val() <= expr2.val(), value_type(0), expr2.fastAccessDx(i) ),
                   if_then_else( expr1.val() <= c.val(), expr1.fastAccessDx(i), value_type(0) ) )


#undef FAD_BINARYOP_MACRO

namespace Sacado {
  namespace Fad {

    template <typename ExprT1, typename ExprT2>
    class MultiplicationOp {};

    template <typename ExprT1, typename ExprT2>
    struct ExprSpec< MultiplicationOp< ExprT1, ExprT2 > > {
      typedef typename ExprSpec<ExprT1>::type type;
    };

    template <typename ExprT1, typename ExprT2>
    class Expr< MultiplicationOp< ExprT1, ExprT2 >,ExprSpecDefault > {

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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive() && expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        return expr1.val()*expr2.val();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        if (expr1.size() > 0 && expr2.size() > 0)
          return expr1.val()*expr2.dx(i) + expr1.dx(i)*expr2.val();
        else if (expr1.size() > 0)
          return expr1.dx(i)*expr2.val();
        else
          return expr1.val()*expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.val()*expr2.fastAccessDx(i) +
          expr1.fastAccessDx(i)*expr2.val();
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    struct ExprSpec< MultiplicationOp< ExprT1, ConstExpr<T2> > > {
      typedef typename ExprSpec<ExprT1>::type type;
    };

    template <typename ExprT1, typename T2>
    class Expr< MultiplicationOp< ExprT1, ConstExpr<T2> >,ExprSpecDefault > {

    public:

      typedef ConstExpr<T2> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr1.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        return expr1.val()*c.val();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return expr1.dx(i)*c.val();
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return expr1.fastAccessDx(i)*c.val();
      }

    protected:

      const ExprT1& expr1;
      ConstT c;
    };

    template <typename T1, typename ExprT2>
    struct ExprSpec< MultiplicationOp< ConstExpr<T1>, ExprT2 > > {
      typedef typename ExprSpec<ExprT2>::type type;
    };

    template <typename T1, typename ExprT2>
    class Expr< MultiplicationOp< ConstExpr<T1>, ExprT2 >,ExprSpecDefault > {

    public:

      typedef ConstExpr<T1> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const ConstT& c_, const ExprT2& expr2_) :
        c(c_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr2.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        return c.val()*expr2.val();
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        return c.val()*expr2.dx(i);
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        return c.val()*expr2.fastAccessDx(i);
      }

    protected:

      ConstT c;
      const ExprT2& expr2;
    };

    template <typename T1, typename T2>
    SACADO_INLINE_FUNCTION
    typename mpl::enable_if_c<
       ExprLevel< Expr<T1> >::value == ExprLevel< Expr<T2> >::value,
       Expr< MultiplicationOp< Expr<T1>, Expr<T2> > >
     >::type
    /*SACADO_FAD_OP_ENABLE_EXPR_EXPR(MultiplicationOp)*/
    operator* (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef MultiplicationOp< Expr<T1>, Expr<T2> > expr_t;

      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< MultiplicationOp< Expr<T>, Expr<T> > >
    operator* (const Expr<T>& expr1, const Expr<T>& expr2)
    {
      typedef MultiplicationOp< Expr<T>, Expr<T> > expr_t;

      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< MultiplicationOp< ConstExpr<typename Expr<T>::value_type>, Expr<T> > >
    operator* (const typename Expr<T>::value_type& c,
               const Expr<T>& expr)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef MultiplicationOp< ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(ConstT(c), expr);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< MultiplicationOp< Expr<T>, ConstExpr<typename Expr<T>::value_type> > >
    operator* (const Expr<T>& expr,
               const typename Expr<T>::value_type& c)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef MultiplicationOp< Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(expr, ConstT(c));
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    SACADO_FAD_OP_ENABLE_SCALAR_EXPR(MultiplicationOp)
    operator* (const typename Expr<T>::scalar_type& c,
               const Expr<T>& expr)
    {
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;
      typedef MultiplicationOp< ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(ConstT(c), expr);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    SACADO_FAD_OP_ENABLE_EXPR_SCALAR(MultiplicationOp)
    operator* (const Expr<T>& expr,
               const typename Expr<T>::scalar_type& c)
    {
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;
      typedef MultiplicationOp< Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(expr, ConstT(c));
    }
  }
}

// Special handling for std::pow() to provide specializations of PowerOp for
// "simd" value types that use if_then_else(). The only reason for not using
// if_then_else() always is to avoid evaluating the derivative if the value is
// zero to avoid throwing FPEs.
namespace Sacado {
  namespace Fad {

    template <typename ExprT1, typename ExprT2, typename Impl>
    class PowerOp {};

    template <typename ExprT1, typename ExprT2>
    struct ExprSpec< PowerOp< ExprT1, ExprT2 > > {
      typedef typename ExprSpec<ExprT1>::type type;
    };

    template <typename ExprT1, typename T2>
    struct ExprSpec<PowerOp< ExprT1, ConstExpr<T2> > > {
      typedef typename ExprSpec<ExprT1>::type type;
    };

    template <typename T1, typename ExprT2>
    struct ExprSpec< PowerOp< ConstExpr<T1>, ExprT2 > > {
      typedef typename ExprSpec<ExprT2>::type type;
    };

    //
    // Implementation for simd type using if_then_else()
    //
    template <typename ExprT1, typename ExprT2>
    class Expr< PowerOp< ExprT1, ExprT2, PowerImpl::Simd >, ExprSpecDefault > {

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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive() && expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        const int sz1 = expr1.size(), sz2 = expr2.size();
        if (sz1 > 0 && sz2 > 0)
          return if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type((expr2.dx(i)*log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*pow(expr1.val(),expr2.val())) );
        else if (sz1 > 0)
          // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
          // It seems less accurate and caused convergence problems in some codes
          return if_then_else( expr2.val() == scalar_type(1.0), expr1.dx(i), if_then_else(expr1.val() == scalar_type(0.0), value_type(0.0), value_type(expr2.val()*expr1.dx(i)/expr1.val()*pow(expr1.val(),expr2.val())) ));
        else
          return if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type(expr2.dx(i)*log(expr1.val())*pow(expr1.val(),expr2.val())) );
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        return if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type((expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val())) );
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< PowerOp< ExprT1, ConstExpr<T2>, PowerImpl::Simd >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T2> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr1.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), c.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using Sacado::if_then_else;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return if_then_else( c.val() == scalar_type(1.0), expr1.dx(i), if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type(c.val()*expr1.dx(i)/expr1.val()*pow(expr1.val(),c.val())) ));
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using Sacado::if_then_else;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return if_then_else( c.val() == scalar_type(1.0), expr1.fastAccessDx(i), if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type(c.val()*expr1.fastAccessDx(i)/expr1.val()*pow(expr1.val(),c.val()))));
      }

    protected:

      const ExprT1& expr1;
      ConstT c;
    };

    template <typename T1, typename ExprT2>
    class Expr< PowerOp< ConstExpr<T1>, ExprT2, PowerImpl::Simd >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T1> ConstT;
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


      SACADO_INLINE_FUNCTION
      Expr(const ConstT& c_, const ExprT2& expr2_) :
        c(c_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr2.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(c.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        return if_then_else( c.val() == scalar_type(0.0), value_type(0.0), value_type(expr2.dx(i)*log(c.val())*pow(c.val(),expr2.val())) );
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        return if_then_else( c.val() == scalar_type(0.0), value_type(0.0), value_type(expr2.fastAccessDx(i)*log(c.val())*pow(c.val(),expr2.val())) );
      }

    protected:

      ConstT c;
      const ExprT2& expr2;
    };

    //
    // Specialization for scalar types using ternary operator
    //

    template <typename ExprT1, typename ExprT2>
    class Expr< PowerOp< ExprT1, ExprT2, PowerImpl::Scalar >,
                ExprSpecDefault > {

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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive() && expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log;
        const int sz1 = expr1.size(), sz2 = expr2.size();
        if (sz1 > 0 && sz2 > 0)
          return expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type((expr2.dx(i)*log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*pow(expr1.val(),expr2.val()));
        else if (sz1 > 0)
          // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
          // It seems less accurate and caused convergence problems in some codes
          return expr2.val() == scalar_type(1.0) ? expr1.dx(i) : expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type(expr2.val()*expr1.dx(i)/expr1.val()*pow(expr1.val(),expr2.val()));
        else
          return expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type(expr2.dx(i)*log(expr1.val())*pow(expr1.val(),expr2.val()));
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type((expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val()));
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< PowerOp< ExprT1, ConstExpr<T2>, PowerImpl::Scalar >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T2> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr1.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), c.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return c.val() == scalar_type(1.0) ? expr1.dx(i) : expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type(c.val()*expr1.dx(i)/expr1.val()*pow(expr1.val(),c.val()));
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return c.val() == scalar_type(1.0) ? expr1.fastAccessDx(i) : expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type(c.val()*expr1.fastAccessDx(i)/expr1.val()*pow(expr1.val(),c.val()));
      }

    protected:

      const ExprT1& expr1;
      ConstT c;
    };

    template <typename T1, typename ExprT2>
    class Expr< PowerOp< ConstExpr<T1>, ExprT2, PowerImpl::Scalar >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T1> ConstT;
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


      SACADO_INLINE_FUNCTION
      Expr(const ConstT& c_, const ExprT2& expr2_) :
        c(c_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr2.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(c.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log;
        return c.val() == scalar_type(0.0) ? value_type(0.0) : value_type(expr2.dx(i)*log(c.val())*pow(c.val(),expr2.val()));
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return c.val() == scalar_type(0.0) ? value_type(0.0) : value_type(expr2.fastAccessDx(i)*log(c.val())*pow(c.val(),expr2.val()));
      }

    protected:

      ConstT c;
      const ExprT2& expr2;
    };

    //
    // Specialization for nested derivatives.  This version does not use
    // if_then_else/ternary-operator on the base so that nested derivatives
    // are correct.
    //
    template <typename ExprT1, typename ExprT2>
    class Expr< PowerOp< ExprT1, ExprT2, PowerImpl::Nested >,
                ExprSpecDefault > {

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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive() && expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log;
        const int sz1 = expr1.size(), sz2 = expr2.size();
        if (sz1 > 0 && sz2 > 0)
          return (expr2.dx(i)*log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*pow(expr1.val(),expr2.val());
        else if (sz1 > 0)
          return expr2.val() == scalar_type(0.0) ? value_type(0.0) : value_type((expr2.val()*expr1.dx(i))*pow(expr1.val(),expr2.val()-scalar_type(1.0)));
        else
          return expr2.dx(i)*log(expr1.val())*pow(expr1.val(),expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return (expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val());
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< PowerOp< ExprT1, ConstExpr<T2>, PowerImpl::Nested >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T2> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr1.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), c.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow;
        return c.val() == scalar_type(0.0) ? value_type(0.0) : value_type(c.val()*expr1.dx(i)*pow(expr1.val(),c.val()-scalar_type(1.0)));
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow;
        return c.val() == scalar_type(0.0) ? value_type(0.0) : value_type(c.val()*expr1.fastAccessDx(i)*pow(expr1.val(),c.val()-scalar_type(1.0)));
      }

    protected:

      const ExprT1& expr1;
      ConstT c;
    };

    template <typename T1, typename ExprT2>
    class Expr< PowerOp< ConstExpr<T1>, ExprT2, PowerImpl::Nested >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T1> ConstT;
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


      SACADO_INLINE_FUNCTION
      Expr(const ConstT& c_, const ExprT2& expr2_) :
        c(c_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr2.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(c.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log;
        return expr2.dx(i)*log(c.val())*pow(c.val(),expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return expr2.fastAccessDx(i)*log(c.val())*pow(c.val(),expr2.val());
      }

    protected:

      ConstT c;
      const ExprT2& expr2;
    };

    //
    // Specialization for nested derivatives.  This version does not use
    // if_then_else/ternary-operator on the base so that nested derivatives
    // are correct.
    //
    template <typename ExprT1, typename ExprT2>
    class Expr< PowerOp< ExprT1, ExprT2, PowerImpl::NestedSimd >,
                ExprSpecDefault > {

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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive() && expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        const int sz1 = expr1.size(), sz2 = expr2.size();
        if (sz1 > 0 && sz2 > 0)
          return (expr2.dx(i)*log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*pow(expr1.val(),expr2.val());
        else if (sz1 > 0)
          return if_then_else( expr2.val() == scalar_type(0.0), value_type(0.0), value_type((expr2.val()*expr1.dx(i))*pow(expr1.val(),expr2.val()-scalar_type(1.0))));
        else
          return expr2.dx(i)*log(expr1.val())*pow(expr1.val(),expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return (expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val());
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< PowerOp< ExprT1, ConstExpr<T2>, PowerImpl::NestedSimd >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T2> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr1.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(expr1.val(), c.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using Sacado::if_then_else;
        return if_then_else( c.val() == scalar_type(0.0), value_type(0.0), value_type(c.val()*expr1.dx(i)*pow(expr1.val(),c.val()-scalar_type(1.0))));
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using Sacado::if_then_else;
        return if_then_else( c.val() == scalar_type(0.0), value_type(0.0), value_type(c.val()*expr1.fastAccessDx(i)*pow(expr1.val(),c.val()-scalar_type(1.0))));
      }

    protected:

      const ExprT1& expr1;
      ConstT c;
    };

    template <typename T1, typename ExprT2>
    class Expr< PowerOp< ConstExpr<T1>, ExprT2, PowerImpl::NestedSimd >,
                ExprSpecDefault > {

    public:

      typedef ConstExpr<T1> ConstT;
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


      SACADO_INLINE_FUNCTION
      Expr(const ConstT& c_, const ExprT2& expr2_) :
        c(c_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr2.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using std::pow;
        return pow(c.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using std::pow; using std::log;
        return expr2.dx(i)*log(c.val())*pow(c.val(),expr2.val());
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return expr2.fastAccessDx(i)*log(c.val())*pow(c.val(),expr2.val());
      }

    protected:

      ConstT c;
      const ExprT2& expr2;
    };

    template <typename T1, typename T2>
    SACADO_INLINE_FUNCTION
    typename mpl::enable_if_c<
       ExprLevel< Expr<T1> >::value == ExprLevel< Expr<T2> >::value,
       Expr< PowerOp< Expr<T1>, Expr<T2> > >
     >::type
    /*SACADO_FAD_OP_ENABLE_EXPR_EXPR(PowerOp)*/
    pow (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef PowerOp< Expr<T1>, Expr<T2> > expr_t;

      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< PowerOp< Expr<T>, Expr<T> > >
    pow (const Expr<T>& expr1, const Expr<T>& expr2)
    {
      typedef PowerOp< Expr<T>, Expr<T> > expr_t;

      return Expr<expr_t>(expr1, expr2);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< PowerOp< ConstExpr<typename Expr<T>::value_type>, Expr<T> > >
    pow (const typename Expr<T>::value_type& c,
         const Expr<T>& expr)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef PowerOp< ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(ConstT(c), expr);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    Expr< PowerOp< Expr<T>, ConstExpr<typename Expr<T>::value_type> > >
    pow (const Expr<T>& expr,
         const typename Expr<T>::value_type& c)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef PowerOp< Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(expr, ConstT(c));
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    SACADO_FAD_OP_ENABLE_SCALAR_EXPR(PowerOp)
    pow (const typename Expr<T>::scalar_type& c,
         const Expr<T>& expr)
    {
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;
      typedef PowerOp< ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(ConstT(c), expr);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    SACADO_FAD_OP_ENABLE_EXPR_SCALAR(PowerOp)
    pow (const Expr<T>& expr,
         const typename Expr<T>::scalar_type& c)
    {
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;
      typedef PowerOp< Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(expr, ConstT(c));
    }
  }
}

//--------------------------if_then_else operator -----------------------
// Can't use the above macros because it is a ternary operator (sort of).
// Also, relies on C++11


namespace Sacado {
  namespace Fad {

    template <typename CondT, typename ExprT1, typename ExprT2>
    class IfThenElseOp {};

    template <typename CondT, typename ExprT1, typename ExprT2>
    struct ExprSpec< IfThenElseOp< CondT, ExprT1, ExprT2 > > {
      typedef typename ExprSpec<ExprT1>::type type;
    };

    template <typename CondT, typename ExprT1, typename ExprT2>
    class Expr< IfThenElseOp< CondT, ExprT1, ExprT2 >,ExprSpecDefault > {

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

      SACADO_INLINE_FUNCTION
      Expr(const CondT& cond_, const ExprT1& expr1_, const ExprT2& expr2_) :
        cond(cond_), expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive() && expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.val(), expr2.val() );
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.dx(i), expr2.dx(i) );
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.fastAccessDx(i), expr2.fastAccessDx(i) );
      }

    protected:

      const CondT&  cond;
      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename CondT, typename ExprT1, typename T2>
    struct ExprSpec< IfThenElseOp< CondT, ExprT1, ConstExpr<T2> > > {
      typedef typename ExprSpec<ExprT1>::type type;
    };

    template <typename CondT, typename ExprT1, typename T2>
    class Expr< IfThenElseOp< CondT, ExprT1, ConstExpr<T2> >,ExprSpecDefault > {

    public:

      typedef ConstExpr<T2> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const CondT& cond_, const ExprT1& expr1_, const ConstT& c_) :
        cond(cond_), expr1(expr1_), c(c_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr1.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.val(), c.val() );
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.dx(i), value_type(0.0) );
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.fastAccessDx(i), value_type(0.0) );
      }

    protected:

      const CondT&  cond;
      const ExprT1& expr1;
      ConstT c;
    };

    template <typename CondT, typename T1, typename ExprT2>
    struct ExprSpec< IfThenElseOp< CondT, ConstExpr<T1>, ExprT2 > > {
      typedef typename ExprSpec<ExprT2>::type type;
    };

    template <typename CondT, typename T1, typename ExprT2>
    class Expr< IfThenElseOp< CondT, ConstExpr<T1>, ExprT2 >,ExprSpecDefault > {

    public:

      typedef ConstExpr<T1> ConstT;
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

      SACADO_INLINE_FUNCTION
      Expr(const CondT& cond_, const ConstT& c_, const ExprT2& expr2_) :
        cond(cond_), c(c_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      bool isPassive() const {
        return expr2.isPassive();
      }

      SACADO_INLINE_FUNCTION
      bool updateValue() const { return expr2.updateValue(); }

      SACADO_INLINE_FUNCTION
      void cache() const {}

      SACADO_INLINE_FUNCTION
      const value_type val() const {
        using Sacado::if_then_else;
        return if_then_else( cond, c.val(), expr2.val() );
      }

      SACADO_INLINE_FUNCTION
      const value_type dx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, value_type(0.0), expr2.dx(i) );
      }

      SACADO_INLINE_FUNCTION
      const value_type fastAccessDx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, value_type(0.0), expr2.fastAccessDx(i) );
      }

    protected:

      const CondT&  cond;
      ConstT c;
      const ExprT2& expr2;
    };

    template <typename CondT, typename T1, typename T2>
    SACADO_INLINE_FUNCTION
    typename mpl::enable_if_c< IsFadExpr<T1>::value && IsFadExpr<T2>::value &&
                               ExprLevel<T1>::value == ExprLevel<T2>::value,
                               Expr< IfThenElseOp< CondT, T1, T2 > >
                             >::type
    if_then_else (const CondT& cond, const T1& expr1, const T2& expr2)
    {
      typedef IfThenElseOp< CondT, T1, T2 > expr_t;

      return Expr<expr_t>(cond, expr1, expr2);
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    Expr< IfThenElseOp< CondT, Expr<T>, Expr<T> > >
    if_then_else (const CondT& cond, const Expr<T>& expr1, const Expr<T>& expr2)
    {
      typedef IfThenElseOp< CondT, Expr<T>, Expr<T> > expr_t;

      return Expr<expr_t>(cond, expr1, expr2);
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    Expr< IfThenElseOp< CondT, ConstExpr<typename Expr<T>::value_type>,
                        Expr<T> > >
    if_then_else (const CondT& cond, const typename Expr<T>::value_type& c,
                  const Expr<T>& expr)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef IfThenElseOp< CondT, ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(cond, ConstT(c), expr);
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    Expr< IfThenElseOp< CondT, Expr<T>,
                        ConstExpr<typename Expr<T>::value_type> > >
    if_then_else (const CondT& cond, const Expr<T>& expr,
                  const typename Expr<T>::value_type& c)
    {
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef IfThenElseOp< CondT, Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(cond, expr, ConstT(c));
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    typename mpl::disable_if<
      std::is_same< typename Expr<T>::value_type,
                    typename Expr<T>::scalar_type>,
      Expr< IfThenElseOp< CondT, ConstExpr<typename Expr<T>::scalar_type>,
                          Expr<T> > >
      >::type
    if_then_else (const CondT& cond, const typename Expr<T>::scalar_type& c,
                  const Expr<T>& expr)
    {
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;
      typedef IfThenElseOp< CondT, ConstT, Expr<T> > expr_t;

      return Expr<expr_t>(cond, ConstT(c), expr);
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    typename mpl::disable_if<
      std::is_same< typename Expr<T>::value_type,
                    typename Expr<T>::scalar_type>,
      Expr< IfThenElseOp< CondT, Expr<T>,
                          ConstExpr<typename Expr<T>::scalar_type> > >
      >::type
    if_then_else (const CondT& cond, const Expr<T>& expr,
                  const typename Expr<T>::scalar_type& c)
    {
      typedef ConstExpr<typename Expr<T>::scalar_type> ConstT;
      typedef IfThenElseOp< CondT, Expr<T>, ConstT > expr_t;

      return Expr<expr_t>(cond, expr, ConstT(c));
    }
  }
}

//-------------------------- Relational Operators -----------------------

namespace Sacado {
  namespace Fad {
    template <typename T1, typename T2 = T1>
    struct ConditionalReturnType {
      typedef decltype( std::declval<T1>() == std::declval<T2>() ) type;
    };
  }
}

#define FAD_RELOP_MACRO(OP)                                             \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
    template <typename ExprT1, typename ExprT2>                         \
    SACADO_INLINE_FUNCTION                                              \
    typename ConditionalReturnType<typename Expr<ExprT1>::value_type,   \
                                   typename Expr<ExprT2>::value_type>::type \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return expr1.val() OP expr2.val();                                \
    }                                                                   \
                                                                        \
    template <typename ExprT2>                                          \
    SACADO_INLINE_FUNCTION                                              \
    typename ConditionalReturnType<typename Expr<ExprT2>::value_type>::type \
    operator OP (const typename Expr<ExprT2>::value_type& a,            \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return a OP expr2.val();                                          \
    }                                                                   \
                                                                        \
    template <typename ExprT1>                                          \
    SACADO_INLINE_FUNCTION                                              \
    typename ConditionalReturnType<typename Expr<ExprT1>::value_type>::type \
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

  namespace Fad {

    template <typename ExprT>
    SACADO_INLINE_FUNCTION
    bool operator ! (const Expr<ExprT>& expr)
    {
      return ! expr.val();
    }

  } // namespace Fad

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace Fad {

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
  namespace Fad {                                                       \
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

  namespace Fad {

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
