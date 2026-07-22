// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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

#ifndef SACADO_FAD_OPS_MP_VECTOR_HPP
#define SACADO_FAD_OPS_MP_VECTOR_HPP

#include "Sacado_Fad_Ops.hpp"
#include "Sacado_mpl_enable_if.hpp"
#include "Sacado_Fad_Expr_MP_Vector.hpp"

#define FAD_UNARYOP_MACRO(OPNAME,OP,MPVALUE,VALUE,DX,FASTACCESSDX)      \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
                                                                        \
    template <typename ExprT>                                           \
    class Expr< OP<ExprT>,ExprSpecMPVector > {                          \
    public:                                                             \
                                                                        \
      typedef typename ExprT::value_type value_type;                    \
      typedef typename ExprT::scalar_type scalar_type;                  \
      typedef typename ExprT::base_expr_type base_expr_type;            \
                                                                        \
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      Expr(const ExprT& expr_) : expr(expr_)  {}                        \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const { return expr.size(); }                          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const { return expr.hasFastAccess(); }       \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool isPassive() const { return expr.isPassive();}                \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool updateValue() const { return expr.updateValue(); }           \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return MPVALUE;                                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type val(int j) const {                                       \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type dx(int i, int j) const {                                 \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type fastAccessDx(int i, int j) const {                       \
        return FASTACCESSDX;                                            \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ExprT& expr;                                                \
    };                                                                  \
                                                                        \
  }                                                                     \
}

FAD_UNARYOP_MACRO(operator+,
                  UnaryPlusOp,
                  expr.val(),
                  expr.val(j),
                  expr.dx(i,j),
                  expr.fastAccessDx(i,j))
FAD_UNARYOP_MACRO(operator-,
                  UnaryMinusOp,
                  -expr.val(),
                  -expr.val(j),
                  -expr.dx(i,j),
                  -expr.fastAccessDx(i,j))
FAD_UNARYOP_MACRO(exp,
                  ExpOp,
                  std::exp(expr.val()),
                  std::exp(expr.val(j)),
                  std::exp(expr.val(j))*expr.dx(i,j),
                  std::exp(expr.val(j))*expr.fastAccessDx(i,j))
FAD_UNARYOP_MACRO(log,
                  LogOp,
                  std::log(expr.val()),
                  std::log(expr.val(j)),
                  expr.dx(i,j)/expr.val(j),
                  expr.fastAccessDx(i,j)/expr.val(j))
FAD_UNARYOP_MACRO(log10,
                  Log10Op,
                  std::log10(expr.val()),
                  std::log10(expr.val(j)),
                  expr.dx(i,j)/( std::log(val_type(10))*expr.val(j)),
                  expr.fastAccessDx(i,j) / ( std::log(val_type(10))*expr.val(j)))
FAD_UNARYOP_MACRO(sqrt,
                  SqrtOp,
                  std::sqrt(expr.val()),
                  std::sqrt(expr.val(j)),
                  expr.dx(i,j)/(val_type(2)* std::sqrt(expr.val(j))),
                  expr.fastAccessDx(i,j)/(val_type(2)* std::sqrt(expr.val(j))))
FAD_UNARYOP_MACRO(cos,
                  CosOp,
                  std::cos(expr.val()),
                  std::cos(expr.val(j)),
                  -expr.dx(i,j)* std::sin(expr.val(j)),
                  -expr.fastAccessDx(i,j)* std::sin(expr.val(j)))
FAD_UNARYOP_MACRO(sin,
                  SinOp,
                  std::sin(expr.val()),
                  std::sin(expr.val(j)),
                  expr.dx(i,j)* std::cos(expr.val(j)),
                  expr.fastAccessDx(i,j)* std::cos(expr.val(j)))
FAD_UNARYOP_MACRO(tan,
                  TanOp,
                  std::tan(expr.val()),
                  std::tan(expr.val(j)),
                  expr.dx(i,j)*
                    (val_type(1)+ std::tan(expr.val(j))* std::tan(expr.val(j))),
                  expr.fastAccessDx(i,j)*
                    (val_type(1)+ std::tan(expr.val(j))* std::tan(expr.val(j))))
FAD_UNARYOP_MACRO(acos,
                  ACosOp,
                  std::acos(expr.val()),
                  std::acos(expr.val(j)),
                  -expr.dx(i,j)/ std::sqrt(val_type(1)-expr.val(j)*expr.val(j)),
                  -expr.fastAccessDx(i,j) /
                    std::sqrt(val_type(1)-expr.val(j)*expr.val(j)))
FAD_UNARYOP_MACRO(asin,
                  ASinOp,
                  std::asin(expr.val()),
                  std::asin(expr.val(j)),
                  expr.dx(i,j)/ std::sqrt(val_type(1)-expr.val(j)*expr.val(j)),
                  expr.fastAccessDx(i,j) /
                    std::sqrt(val_type(1)-expr.val(j)*expr.val(j)))
FAD_UNARYOP_MACRO(atan,
                  ATanOp,
                  std::atan(expr.val()),
                  std::atan(expr.val(j)),
                  expr.dx(i,j)/(val_type(1)+expr.val(j)*expr.val(j)),
                  expr.fastAccessDx(i,j)/(val_type(1)+expr.val(j)*expr.val(j)))
FAD_UNARYOP_MACRO(cosh,
                  CoshOp,
                  std::cosh(expr.val()),
                  std::cosh(expr.val(j)),
                  expr.dx(i,j)* std::sinh(expr.val(j)),
                  expr.fastAccessDx(i,j)* std::sinh(expr.val(j)))
FAD_UNARYOP_MACRO(sinh,
                  SinhOp,
                  std::sinh(expr.val()),
                  std::sinh(expr.val(j)),
                  expr.dx(i,j)* std::cosh(expr.val(j)),
                  expr.fastAccessDx(i,j)* std::cosh(expr.val(j)))
FAD_UNARYOP_MACRO(tanh,
                  TanhOp,
                  std::tanh(expr.val()),
                  std::tanh(expr.val(j)),
                  expr.dx(i,j)/( std::cosh(expr.val(j))* std::cosh(expr.val(j))),
                  expr.fastAccessDx(i,j) /
                    ( std::cosh(expr.val(j))* std::cosh(expr.val(j))))
FAD_UNARYOP_MACRO(acosh,
                  ACoshOp,
                  std::acosh(expr.val()),
                  std::acosh(expr.val(j)),
                  expr.dx(i,j)/ std::sqrt((expr.val(j)-val_type(1)) *
                                       (expr.val(j)+val_type(1))),
                  expr.fastAccessDx(i,j)/ std::sqrt((expr.val(j)-val_type(1)) *
                                                 (expr.val(j)+val_type(1))))
FAD_UNARYOP_MACRO(asinh,
                  ASinhOp,
                  std::asinh(expr.val()),
                  std::asinh(expr.val(j)),
                  expr.dx(i,j)/ std::sqrt(val_type(1)+expr.val(j)*expr.val(j)),
                  expr.fastAccessDx(i,j)/ std::sqrt(val_type(1)+
                                                 expr.val(j)*expr.val(j)))
FAD_UNARYOP_MACRO(atanh,
                  ATanhOp,
                  std::atanh(expr.val()),
                  std::atanh(expr.val(j)),
                  expr.dx(i,j)/(val_type(1)-expr.val(j)*expr.val(j)),
                  expr.fastAccessDx(i,j)/(val_type(1)-
                                                 expr.val(j)*expr.val(j)))
FAD_UNARYOP_MACRO(abs,
                  AbsOp,
                  std::abs(expr.val()),
                  std::abs(expr.val(j)),
                  expr.val(j) >= 0 ? val_type(+expr.dx(i,j)) :
                    val_type(-expr.dx(i,j)),
                  expr.val(j) >= 0 ? val_type(+expr.fastAccessDx(i,j)) :
                    val_type(-expr.fastAccessDx(i,j)))
FAD_UNARYOP_MACRO(fabs,
                  FAbsOp,
                  std::fabs(expr.val()),
                  std::fabs(expr.val(j)),
                  expr.val(j) >= 0 ? val_type(+expr.dx(i,j)) :
                    val_type(-expr.dx(i,j)),
                  expr.val(j) >= 0 ? val_type(+expr.fastAccessDx(i,j)) :
                    val_type(-expr.fastAccessDx(i,j)))
FAD_UNARYOP_MACRO(cbrt,
                  CbrtOp,
                  std::cbrt(expr.val()),
                  std::cbrt(expr.val(j)),
                  expr.dx(i,j)/(val_type(3)*std::cbrt(expr.val(j)*expr.val(j))),
                  expr.fastAccessDx(i,j)/(val_type(3)*std::cbrt(expr.val(j)*expr.val(j))))

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,MPVALUE,VALUE,DX,FASTACCESSDX,MPVAL_CONST_DX_1,MPVAL_CONST_DX_2,VAL_CONST_DX_1,VAL_CONST_DX_2,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
                                                                        \
    template <typename ExprT1, typename ExprT2>                         \
    class Expr< OP< ExprT1, ExprT2 >,ExprSpecMPVector > {               \
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
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :                \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        int sz1 = expr1.size(), sz2 = expr2.size();                     \
        return sz1 > sz2 ? sz1 : sz2;                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess() && expr2.hasFastAccess();          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool isPassive() const {                                          \
        return expr1.isPassive() && expr2.isPassive();                  \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool updateValue() const {                                        \
        return expr1.updateValue() && expr2.updateValue();              \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const value_type val() const {                                    \
        return MPVALUE;                                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type val(int j) const {                                 \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type dx(int i, int j) const {                           \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type fastAccessDx(int i, int j) const {                 \
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
    class Expr< OP< ExprT1, ConstExpr<T2> >,ExprSpecMPVector > {        \
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
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      Expr(const ExprT1& expr1_, const ConstT& c_) :                    \
        expr1(expr1_), c(c_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr1.size();                                            \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool isPassive() const {                                          \
        return expr1.isPassive();                                       \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool updateValue() const { return expr1.updateValue(); }          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const value_type val() const {                                    \
        return MPVAL_CONST_DX_2;                                        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type val(int j) const {                                 \
        return VAL_CONST_DX_2;                                          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type dx(int i, int j) const {                           \
        return CONST_DX_2;                                              \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type fastAccessDx(int i, int j) const {                 \
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
    class Expr< OP< ConstExpr<T1>, ExprT2 >,ExprSpecMPVector > {        \
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
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      Expr(const ConstT& c_, const ExprT2& expr2_) :                    \
        c(c_), expr2(expr2_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr2.size();                                            \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr2.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool isPassive() const {                                          \
        return expr2.isPassive();                                       \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool updateValue() const { return expr2.updateValue(); }          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const value_type val() const {                                    \
        return MPVAL_CONST_DX_1;                                        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type val(int j) const {                                 \
        return VAL_CONST_DX_1;                                          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type dx(int i, int j) const {                           \
        return CONST_DX_1;                                              \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      const val_type fastAccessDx(int i, int j) const {                 \
        return CONST_FASTACCESSDX_1;                                    \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      ConstT c;                                                         \
      const ExprT2& expr2;                                              \
    };                                                                  \
                                                                        \
  }                                                                     \
}


FAD_BINARYOP_MACRO(operator+,
                   AdditionOp,
                   expr1.val() + expr2.val(),
                   expr1.val(j) + expr2.val(j),
                   expr1.dx(i,j) + expr2.dx(i,j),
                   expr1.fastAccessDx(i,j) + expr2.fastAccessDx(i,j),
                   c.val() + expr2.val(),
                   expr1.val() + c.val(),
                   c.val(j) + expr2.val(j),
                   expr1.val(j) + c.val(j),
                   expr2.dx(i,j),
                   expr1.dx(i,j),
                   expr2.fastAccessDx(i,j),
                   expr1.fastAccessDx(i,j))
FAD_BINARYOP_MACRO(operator-,
                   SubtractionOp,
                   expr1.val() - expr2.val(),
                   expr1.val(j) - expr2.val(j),
                   expr1.dx(i,j) - expr2.dx(i,j),
                   expr1.fastAccessDx(i,j) - expr2.fastAccessDx(i,j),
                   c.val() - expr2.val(),
                   expr1.val() - c.val(),
                   c.val(j) - expr2.val(j),
                   expr1.val(j) - c.val(j),
                   -expr2.dx(i,j),
                   expr1.dx(i,j),
                   -expr2.fastAccessDx(i,j),
                   expr1.fastAccessDx(i,j))
// FAD_BINARYOP_MACRO(operator*,
//                 MultiplicationOp,
//                 expr1.val() * expr2.val(),
//                 expr1.val(j) * expr2.val(j),
//                 expr1.val(j)*expr2.dx(i,j) + expr1.dx(i,j)*expr2.val(j),
//                 expr1.val(j)*expr2.fastAccessDx(i,j) +
//                   expr1.fastAccessDx(i,j)*expr2.val(j),
//                 c.val() * expr2.val(),
//                 expr1.val() * c.val(),
//                 c.val(j) * expr2.val(j),
//                 expr1.val(j) * c.val(j),
//                 c.val(j)*expr2.dx(i,j),
//                 expr1.dx(i,j)*c.val(j),
//                 c.val(j)*expr2.fastAccessDx(i,j),
//                 expr1.fastAccessDx(i,j)*c.val(j))
FAD_BINARYOP_MACRO(operator/,
                   DivisionOp,
                   expr1.val() / expr2.val(),
                   expr1.val(j) / expr2.val(j),
                   (expr1.dx(i,j)*expr2.val(j) - expr2.dx(i,j)*expr1.val(j)) /
                     (expr2.val(j)*expr2.val(j)),
                   (expr1.fastAccessDx(i,j)*expr2.val(j) -
                      expr2.fastAccessDx(i,j)*expr1.val(j)) /
                      (expr2.val(j)*expr2.val(j)),
                   c.val() / expr2.val(),
                   expr1.val() / c.val(),
                   c.val(j) / expr2.val(j),
                   expr1.val(j) / c.val(j),
                   -expr2.dx(i,j)*c.val(j) / (expr2.val(j)*expr2.val(j)),
                   expr1.dx(i,j)/c.val(j),
                   -expr2.fastAccessDx(i,j)*c.val(j) / (expr2.val(j)*expr2.val(j)),
                   expr1.fastAccessDx(i,j)/c.val(j))
FAD_BINARYOP_MACRO(atan2,
                   Atan2Op,
                   std::atan2(expr1.val(), expr2.val()),
                   std::atan2(expr1.val(j), expr2.val(j)),
                   (expr2.val(j)*expr1.dx(i,j) - expr1.val(j)*expr2.dx(i,j))/
                        (expr1.val(j)*expr1.val(j) + expr2.val(j)*expr2.val(j)),
                   (expr2.val(j)*expr1.fastAccessDx(i,j) - expr1.val(j)*expr2.fastAccessDx(i,j))/
                        (expr1.val(j)*expr1.val(j) + expr2.val(j)*expr2.val(j)),
                   std::atan2(c.val(), expr2.val()),
                   std::atan2(expr1.val(), c.val()),
                   std::atan2(c.val(j), expr2.val(j)),
                   std::atan2(expr1.val(j), c.val(j)),
                   (-c.val(j)*expr2.dx(i,j)) / (c.val(j)*c.val(j) + expr2.val(j)*expr2.val(j)),
                   (c.val(j)*expr1.dx(i,j))/ (expr1.val(j)*expr1.val(j) + c.val(j)*c.val(j)),
                   (-c.val(j)*expr2.fastAccessDx(i,j))/ (c.val(j)*c.val(j) + expr2.val(j)*expr2.val(j)),
                   (c.val(j)*expr1.fastAccessDx(i,j))/ (expr1.val(j)*expr1.val(j) + c.val(j)*c.val(j)))
FAD_BINARYOP_MACRO(pow,
                   PowerOp,
                   std::pow(expr1.val(), expr2.val()),
                   std::pow(expr1.val(j), expr2.val(j)),
                   expr1.val(j) == val_type(0) ? val_type(0) : val_type((expr2.dx(i,j)*std::log(expr1.val(j))+expr2.val(j)*expr1.dx(i,j)/expr1.val(j))*std::pow(expr1.val(j),expr2.val(j))),
                   expr1.val(j) == val_type(0) ? val_type(0.0) : val_type((expr2.fastAccessDx(i,j)*std::log(expr1.val(j))+expr2.val(j)*expr1.fastAccessDx(i,j)/expr1.val(j))*std::pow(expr1.val(j),expr2.val(j))),
                   std::pow(c.val(), expr2.val()),
                   std::pow(expr1.val(), c.val()),
                   std::pow(c.val(j), expr2.val(j)),
                   std::pow(expr1.val(j), c.val(j)),
                   c.val(j) == val_type(0) ? val_type(0) : val_type(expr2.dx(i,j)*std::log(c.val(j))*std::pow(c.val(j),expr2.val(j))),
                   expr1.val(j) == val_type(0) ? val_type(0.0) : val_type(c.val(j)*expr1.dx(i,j)/expr1.val(j)*std::pow(expr1.val(j),c.val(j))),
                   c.val(j) == val_type(0) ? val_type(0) : val_type(expr2.fastAccessDx(i,j)*std::log(c.val(j))*std::pow(c.val(j),expr2.val(j))),
                   expr1.val(j) == val_type(0) ? val_type(0.0) : val_type(c.val(j)*expr1.fastAccessDx(i,j)/expr1.val(j)*std::pow(expr1.val(j),c.val(j))))
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   std::max(expr1.val(), expr2.val()),
                   std::max(expr1.val(j), expr2.val(j)),
                   expr1.val(j) >= expr2.val(j) ? expr1.dx(i,j) : expr2.dx(i,j),
                   expr1.val(j) >= expr2.val(j) ? expr1.fastAccessDx(i,j) :
                                                expr2.fastAccessDx(i,j),
                   std::max(c.val(), expr2.val()),
                   std::max(expr1.val(), c.val()),
                   std::max(c.val(j), expr2.val(j)),
                   std::max(expr1.val(j), c.val(j)),
                   c.val(j) >= expr2.val(j) ? val_type(0) : expr2.dx(i,j),
                   expr1.val(j) >= c.val(j) ? expr1.dx(i,j) : val_type(0),
                   c.val(j) >= expr2.val(j) ? val_type(0) : expr2.fastAccessDx(i,j),
                   expr1.val(j) >= c.val(j) ? expr1.fastAccessDx(i,j) : val_type(0))
FAD_BINARYOP_MACRO(min,
                   MinOp,
                   std::min(expr1.val(), expr2.val()),
                   std::min(expr1.val(j), expr2.val(j)),
                   expr1.val(j) <= expr2.val(j) ? expr1.dx(i,j) : expr2.dx(i,j),
                   expr1.val(j) <= expr2.val(j) ? expr1.fastAccessDx(i,j) :
                                                expr2.fastAccessDx(i,j),
                   std::min(c.val(), expr2.val()),
                   std::min(expr1.val(), c.val()),
                   std::min(c.val(j), expr2.val(j)),
                   std::min(expr1.val(j), c.val(j)),
                   c.val(j) <= expr2.val(j) ? val_type(0) : expr2.dx(i,j),
                   expr1.val(j) <= c.val(j) ? expr1.dx(i,j) : val_type(0),
                   c.val(j) <= expr2.val(j) ? val_type(0) : expr2.fastAccessDx(i,j),
                   expr1.val(j) <= c.val(j) ? expr1.fastAccessDx(i,j) : val_type(0))


#undef FAD_BINARYOP_MACRO

namespace Sacado {
  namespace Fad {

    template <typename ExprT1, typename ExprT2>
    class Expr< MultiplicationOp< ExprT1, ExprT2 >,ExprSpecMPVector > {

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

      typedef typename value_type::value_type val_type;

      KOKKOS_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ExprT2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive() && expr2.isPassive();
      }

      KOKKOS_INLINE_FUNCTION
      bool updateValue() const {
        return expr1.updateValue() && expr2.updateValue();
      }

      KOKKOS_INLINE_FUNCTION
      const value_type val() const {
        return expr1.val()*expr2.val();
      }

      KOKKOS_INLINE_FUNCTION
      const val_type val(int j) const {
        return expr1.val(j)*expr2.val(j);
      }

      KOKKOS_INLINE_FUNCTION
      const val_type dx(int i, int j) const {
        if (expr1.size() > 0 && expr2.size() > 0)
          return expr1.val(j)*expr2.dx(i,j) + expr1.dx(i,j)*expr2.val(j);
        else if (expr1.size() > 0)
          return expr1.dx(i,j)*expr2.val(j);
        else
          return expr1.val(j)*expr2.dx(i,j);
      }

      KOKKOS_INLINE_FUNCTION
      const val_type fastAccessDx(int i, int j) const {
        return expr1.val(j)*expr2.fastAccessDx(i,j) +
          expr1.fastAccessDx(i,j)*expr2.val(j);
      }

    protected:

      const ExprT1& expr1;
      const ExprT2& expr2;

    };

    template <typename ExprT1, typename T2>
    class Expr< MultiplicationOp< ExprT1, ConstExpr<T2> >,ExprSpecMPVector > {

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

      typedef typename value_type::value_type val_type;

      KOKKOS_INLINE_FUNCTION
      Expr(const ExprT1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      bool isPassive() const {
        return expr1.isPassive();
      }

      KOKKOS_INLINE_FUNCTION
      bool updateValue() const { return expr1.updateValue(); }

       KOKKOS_INLINE_FUNCTION
      const value_type val() const {
        return expr1.val()*c.val();
      }

      KOKKOS_INLINE_FUNCTION
      const val_type val(int j) const {
        return expr1.val(j)*c.val(j);
      }

      KOKKOS_INLINE_FUNCTION
      const val_type dx(int i, int j) const {
        return expr1.dx(i,j)*c.val(j);
      }

      KOKKOS_INLINE_FUNCTION
      const val_type fastAccessDx(int i, int j) const {
        return expr1.fastAccessDx(i,j)*c.val(j);
      }

    protected:

      const ExprT1& expr1;
      ConstT c;
    };

    template <typename T1, typename ExprT2>
    class Expr< MultiplicationOp< ConstExpr<T1>, ExprT2 >,ExprSpecMPVector > {

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

      typedef typename value_type::value_type val_type;

      KOKKOS_INLINE_FUNCTION
      Expr(const ConstT& c_, const ExprT2& expr2_) :
        c(c_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      bool isPassive() const {
        return expr2.isPassive();
      }

      KOKKOS_INLINE_FUNCTION
      bool updateValue() const { return expr2.updateValue(); }

      KOKKOS_INLINE_FUNCTION
      const value_type val() const {
        return c.val()*expr2.val();
      }

      KOKKOS_INLINE_FUNCTION
      const val_type val(int j) const {
        return c.val(j)*expr2.val(j);
      }

      KOKKOS_INLINE_FUNCTION
      const val_type dx(int i, int j) const {
        return c.val(j)*expr2.dx(i,j);
      }

      KOKKOS_INLINE_FUNCTION
      const val_type fastAccessDx(int i, int j) const {
        return c.val(j)*expr2.fastAccessDx(i,j);
      }

    protected:

      ConstT c;
      const ExprT2& expr2;
    };

    template <typename ExprT>
    KOKKOS_INLINE_FUNCTION
    bool toBool(const Expr<ExprT,ExprSpecMPVector>& x) {
      bool is_zero = (x.val() == 0.0);
      for (int i=0; i<x.size(); i++)
        for (int j=0; j<x.val().size(); ++j)
          is_zero = is_zero && (x.dx(i,j) == 0.0);
      return !is_zero;
    }

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT,ExprSpecMPVector>& x) {
      os << x.val() << " [";

      for (int i=0; i< x.size(); i++) {
        os << " [";
        for (int j=0; j<x.val().size(); ++j) {
          os << " " << x.dx(i,j);
        }
        os << " ]";
      }

      os << " ]";
      return os;
    }
  }
}

#endif // SACADO_FAD_OPS_MP_VECTOR_HPP
