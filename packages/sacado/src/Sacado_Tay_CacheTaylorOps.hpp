// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TAY_CACHETAYLOROPS_HPP
#define SACADO_TAY_CACHETAYLOROPS_HPP

#include "Sacado_Tay_CacheTaylorExpr.hpp"

#include <cmath>
#include <valarray>
#include <algorithm>    // for std::min and std::max
#include <ostream>      // for std::ostream

namespace Sacado {

  namespace Tay {

    // ---------------------- Unary Addition operator ------------------------

    template <typename ExprT>
    class UnaryPlusOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      UnaryPlusOp(const ExprT& expr) {}

      void allocateCache(int d) const {}

      value_type computeCoeff(int i, const ExprT& expr) const {
        return expr.coeff(i);
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const {
        return expr.fastAccessCoeff(i);
      }

    }; // class UnaryPlusOp

    // ---------------------- Unary Subtraction operator ---------------------

    template <typename ExprT>
    class UnaryMinusOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      UnaryMinusOp(const ExprT& expr) {}

      void allocateCache(int d) const {}

      value_type computeCoeff(int i, const ExprT& expr) const {
        return -expr.coeff(i);
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const {
        return -expr.fastAccessCoeff(i);
      }

    }; // class UnaryPlusOp

    // -------------------------- exp() function -----------------------------

    template <typename ExprT>
    class ExpOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      ExpOp(const ExprT& expr) :
        c(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::exp(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*c[k-j]*expr.coeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::exp(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*c[k-j]*expr.fastAccessCoeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ExpOp

    // -------------------------- log() function -----------------------------

    template <typename ExprT>
    class LogOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      LogOp(const ExprT& expr) :
        c(),
        dc(-1)
      {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::log(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            c[k] = value_type(k)*expr.coeff(k);
            for (int j=1; j<=k-1; j++)
              c[k] -= value_type(j)*expr.coeff(k-j)*c[j];
            c[k] /= (value_type(k)*expr.fastAccessCoeff(0));
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::log(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            c[k] = value_type(k)*expr.fastAccessCoeff(k);
            for (int j=1; j<=k-1; j++)
              c[k] -= value_type(j)*expr.fastAccessCoeff(k-j)*c[j];
            c[k] /= (value_type(k)*expr.fastAccessCoeff(0));
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class LogOp

    // -------------------------- sqrt() function -----------------------------

    template <typename ExprT>
    class SqrtOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      SqrtOp(const ExprT& expr) :
        c(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::sqrt(expr.fastAccessCoeff(0));
            dc = 0;
          }
          value_type tmp = value_type(2)*c[0];
          for (int k=dc+1; k<=i; k++) {
            c[k] = expr.coeff(k);
            for (int j=1; j<=k-1; j++)
              c[k] -= c[j]*c[k-j];
            c[k] /= tmp;
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::sqrt(expr.fastAccessCoeff(0));
            dc = 0;
          }
          value_type tmp = value_type(2)*c[0];
          for (int k=dc+1; k<=i; k++) {
            c[k] = expr.fastAccessCoeff(k);
            for (int j=1; j<=k-1; j++)
              c[k] -= c[j]*c[k-j];
            c[k] /= tmp;
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class SqrtOp

    // -------------------------- cos() function -----------------------------

    template <typename ExprT>
    class CosOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      CosOp(const ExprT& expr) :
        c(),
        s(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
        s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cos(expr.fastAccessCoeff(0));
            s[0] = std::sin(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] -= value_type(j)*expr.coeff(j)*s[k-j];
              s[k] += value_type(j)*expr.coeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cos(expr.fastAccessCoeff(0));
            s[0] = std::sin(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] -= value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
              s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class CosOp

    // -------------------------- sin() function -----------------------------

    template <typename ExprT>
    class SinOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      SinOp(const ExprT& expr) :
        c(),
        s(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
        s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cos(expr.fastAccessCoeff(0));
            s[0] = std::sin(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] -= value_type(j)*expr.coeff(j)*s[k-j];
              s[k] += value_type(j)*expr.coeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return s[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cos(expr.fastAccessCoeff(0));
            s[0] = std::sin(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] -= value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
              s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return s[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class SinOp

    // -------------------------- cosh() function -----------------------------

    template <typename ExprT>
    class CoshOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      CoshOp(const ExprT& expr) :
        c(),
        s(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
        s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cosh(expr.fastAccessCoeff(0));
            s[0] = std::sinh(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] += value_type(j)*expr.coeff(j)*s[k-j];
              s[k] += value_type(j)*expr.coeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cosh(expr.fastAccessCoeff(0));
            s[0] = std::sinh(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] += value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
              s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class CoshOp

    // -------------------------- sinh() function -----------------------------

    template <typename ExprT>
    class SinhOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      SinhOp(const ExprT& expr) :
        c(),
        s(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
        s.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cosh(expr.fastAccessCoeff(0));
            s[0] = std::sinh(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] += value_type(j)*expr.coeff(j)*s[k-j];
              s[k] += value_type(j)*expr.coeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return s[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::cosh(expr.fastAccessCoeff(0));
            s[0] = std::sinh(expr.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++) {
              c[k] += value_type(j)*expr.fastAccessCoeff(j)*s[k-j];
              s[k] += value_type(j)*expr.fastAccessCoeff(j)*c[k-j];
            }
            c[k] /= value_type(k);
            s[k] /= value_type(k);
          }
          dc = i;
        }
        return s[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable std::valarray<value_type> s;
      mutable int dc;

    }; // class SinhOp

    // -------------------------- fabs() function -----------------------------

    template <typename ExprT>
    class FAbsOp {
    public:

      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;
      typedef typename ExprT::base_expr_type base_expr_type;

      FAbsOp(const ExprT& expr) {}

      void allocateCache(int d) const {}

      value_type computeCoeff(int i, const ExprT& expr) const {
        if (expr.fastAccessCoeff(0) > 0)
          return expr.coeff(i);
        else
          return -expr.coeff(i);
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT& expr) const
      {
        if (expr.fastAccessCoeff(0) > 0)
          return expr.fastAccessCoeff(i);
        else
          return -expr.fastAccessCoeff(i);
      }

    }; // class FAbsOp

  } // namespace Tay

} // namespace Sacado

#define TAYLOR_UNARYOP_MACRO(OPNAME,OP)                                 \
namespace Sacado {                                                      \
  namespace Tay {                                                       \
    template <typename T>                                               \
    inline Expr< UnaryExpr< Expr<T>, OP > >                             \
    OPNAME (const Expr<T>& expr)                                        \
    {                                                                   \
      typedef UnaryExpr< Expr<T>, OP > expr_t;                          \
                                                                        \
      return Expr<expr_t>(expr_t(expr));                                \
    }                                                                   \
  }                                                                     \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Tay::OPNAME;                                            \
}

TAYLOR_UNARYOP_MACRO(operator+, UnaryPlusOp)
TAYLOR_UNARYOP_MACRO(operator-, UnaryMinusOp)
TAYLOR_UNARYOP_MACRO(exp, ExpOp)
TAYLOR_UNARYOP_MACRO(log, LogOp)
TAYLOR_UNARYOP_MACRO(sqrt, SqrtOp)
TAYLOR_UNARYOP_MACRO(cos, CosOp)
TAYLOR_UNARYOP_MACRO(sin, SinOp)
TAYLOR_UNARYOP_MACRO(cosh, CoshOp)
TAYLOR_UNARYOP_MACRO(sinh, SinhOp)
TAYLOR_UNARYOP_MACRO(abs, FAbsOp)
TAYLOR_UNARYOP_MACRO(fabs, FAbsOp)

#undef TAYLOR_UNARYOP_MACRO

namespace Sacado {

  namespace Tay {

    // ---------------------- Addition operator -----------------------------

    template <typename ExprT1, typename ExprT2>
    class AdditionOp {
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

      AdditionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        return expr1.coeff(i) + expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        return expr1.fastAccessCoeff(i) + expr2.fastAccessCoeff(i);
      }

    }; // class AdditionOp

    template <typename ExprT1>
    class AdditionOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;
      typedef typename ExprT1::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      AdditionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return expr1.coeff(i) + expr2.coeff(i);
        else
          return expr1.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return expr1.fastAccessCoeff(i) + expr2.fastAccessCoeff(i);
        else
          return expr1.fastAccessCoeff(i);
      }

    }; // class AdditionOp

    template <typename ExprT2>
    class AdditionOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;
      typedef typename ExprT2::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      AdditionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return expr1.coeff(i) + expr2.coeff(i);
        else
          return expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return expr1.fastAccessCoeff(i) + expr2.fastAccessCoeff(i);
        else
          return expr2.fastAccessCoeff(i);
      }

    }; // class AdditionOp

    // ---------------------- Subtraction operator ---------------------------

    template <typename ExprT1, typename ExprT2>
    class SubtractionOp {
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

      SubtractionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        return expr1.coeff(i) - expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        return expr1.fastAccessCoeff(i) - expr2.fastAccessCoeff(i);
      }

    }; // class SubtractionOp

    template <typename ExprT1>
    class SubtractionOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;
      typedef typename ExprT1::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      SubtractionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return expr1.coeff(i) - expr2.coeff(i);
        else
          return expr1.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return expr1.fastAccessCoeff(i) - expr2.fastAccessCoeff(i);
        else
          return expr1.fastAccessCoeff(i);
      }

    }; // class SubtractionOp

    template <typename ExprT2>
    class SubtractionOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;
      typedef typename ExprT2::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      SubtractionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return expr1.coeff(i) - expr2.coeff(i);
        else
          return -expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return expr1.fastAccessCoeff(i) - expr2.fastAccessCoeff(i);
        else
          return -expr2.fastAccessCoeff(i);
      }

    }; // class SubtractionOp

    // ---------------------- Multiplication operator -------------------------

    template <typename ExprT1, typename ExprT2>
    class MultiplicationOp {
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

      MultiplicationOp(const ExprT1& expr1, const ExprT2 expr2) :
        c(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          for (int k=dc+1; k<=i; k++) {
            for (int j=0; j<=k; j++)
              c[k] += expr1.coeff(j)*expr2.coeff(k-j);
          }
          dc = i;
        }
        return c[i];
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          for (int k=dc+1; k<=i; k++) {
            for (int j=0; j<=k; j++)
              c[k] += expr1.fastAccessCoeff(j)*expr2.fastAccessCoeff(k-j);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class MultiplicationOp

    template <typename ExprT1>
    class MultiplicationOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;
      typedef typename ExprT1::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      MultiplicationOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        return expr1.coeff(i)*expr2.value();
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        return expr1.fastAccessCoeff(i)*expr2.value();
      }

    }; // class MultiplicationOp

    template <typename ExprT2>
    class MultiplicationOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;
      typedef typename ExprT2::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      MultiplicationOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        return expr1.value()*expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        return expr1.value()*expr2.fastAccessCoeff(i);
      }

    }; // class MultiplicationOp

    // ---------------------- Division operator -------------------------

    template <typename ExprT1, typename ExprT2>
    class DivisionOp {
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

      DivisionOp(const ExprT1& expr1, const ExprT2 expr2) :
        c(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          for (int k=dc+1; k<=i; k++) {
            c[k] = expr1.coeff(k);
            for (int j=1; j<=k; j++)
              c[k] -= expr2.coeff(j)*c[k-j];
            c[k] /= expr2.fastAccessCoeff(0);
          }
          dc = i;
        }
        return c[i];
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          for (int k=dc+1; k<=i; k++) {
            c[k] = expr1.coeff(k);
            for (int j=1; j<=k; j++)
              c[k] -= expr2.fastAccessCoeff(j)*c[k-j];
            c[k] /= expr2.fastAccessCoeff(0);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class DivisionOp

    template <typename ExprT1>
    class DivisionOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;
      typedef typename ExprT1::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      DivisionOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        return expr1.coeff(i)/expr2.value();
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        return expr1.fastAccessCoeff(i)/expr2.value();
      }

    }; // class DivisionOp

    template <typename ExprT2>
    class DivisionOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;
      typedef typename ExprT2::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      DivisionOp(const ExprT1& expr1, const ExprT2 expr2) :
        c(),
        dc(-1) {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = expr1.fastAccessCoeff(0) / expr2.fastAccessCoeff(0);
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] -= expr2.coeff(j)*c[k-j];
            c[k] /= expr2.fastAccessCoeff(0);
          }
          dc = i;
        }
        return c[i];
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = expr1.fastAccessCoeff(0) / expr2.fastAccessCoeff(0);
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            c[k] = expr1.coeff(k);
            for (int j=1; j<=k; j++)
              c[k] -= expr2.fastAccessCoeff(j)*c[k-j];
            c[k] /= expr2.fastAccessCoeff(0);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class DivisionOp

    // ---------------------- Max operator -----------------------------

    template <typename ExprT1, typename ExprT2>
    class MaxOp {
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

      MaxOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return std::max(expr1.coeff(0), expr2.coeff(0));
        else
          return expr1.coeff(0) >= expr2.coeff(0) ? expr1.coeff(i) :
            expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return std::max(expr1.fastAccessCoeff(0), expr2.fastAccessCoeff(0));
        else
          return expr1.fastAccessCoeff(0) >= expr2.fastAccessCoeff(0) ?
            expr1.fastAccessoeff(i) : expr2.fastAccessCoeff(i);
      }

    }; // class MaxOp

    template <typename ExprT1>
    class MaxOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;
      typedef typename ExprT1::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      MaxOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return std::max(expr1.coeff(0), expr2.value());
        else
          return expr1.coeff(0) >= expr2.value() ? expr1.coeff(i) :
            value_type(0);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return std::max(expr1.fastAccessCoeff(0), expr2.value());
        else
          return expr1.fastAccessCoeff(0) >= expr2.value() ?
            expr1.fastAccessCoeff(i) : value_type(0);
      }

    }; // class MaxOp

    template <typename ExprT2>
    class MaxOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;
      typedef typename ExprT2::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      MaxOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return std::max(expr1.value(), expr2.coeff(0));
        else
          return expr1.value() >= expr2.coeff(0) ? value_type(0) :
            expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return std::max(expr1.value(), expr2.fastAccessCoeff(0));
        else
          return expr1.value() >= expr2.fastAccessCoeff(0) ? value_type(0) :
            expr2.fastAccessCoeff(i);
      }

    }; // class MaxOp

    // ---------------------- Min operator -----------------------------

    template <typename ExprT1, typename ExprT2>
    class MinOp {
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

      MinOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return min(expr1.coeff(0), expr2.coeff(0));
        else
          return expr1.coeff(0) <= expr2.coeff(0) ? expr1.coeff(i) :
            expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return min(expr1.fastAccessCoeff(0), expr2.fastAccessCoeff(0));
        else
          return expr1.fastAccessCoeff(0) <= expr2.fastAccessCoeff(0) ?
            expr1.fastAccessCoeff(i) : expr2.fastAccessCoeff(i);
      }

    }; // class MinOp

    template <typename ExprT1>
    class MinOp<ExprT1, ConstExpr<typename ExprT1::value_type> > {
    public:

      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;
      typedef typename ExprT1::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT1::value_type> ExprT2;

      MinOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return std::min(expr1.coeff(0), expr2.value());
        else
          return expr1.coeff(0) <= expr2.value() ? expr1.coeff(i) :
            value_type(0);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return std::min(expr1.fastAccessCoeff(0), expr2.value());
        else
          return expr1.fastAccessCoeff(0) <= expr2.value() ?
            expr1.fastAccessCoeff(i) : value_type(0);
      }

    }; // class MinOp

    template <typename ExprT2>
    class MinOp<ConstExpr<typename ExprT2::value_type>, ExprT2 > {
    public:

      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;
      typedef typename ExprT2::base_expr_type base_expr_type;
      typedef ConstExpr<typename ExprT2::value_type> ExprT1;

      MinOp(const ExprT1& expr1, const ExprT2 expr2) {}

      void allocateCache(int d) const {}

      value_type
      computeCoeff(int i, const ExprT1& expr1,
                   const ExprT2& expr2) const {
        if (i == 0)
          return std::min(expr1.value(), expr2.coeff(0));
        else
          return expr1.value() <= expr2.coeff(0) ? value_type(0) :
            expr2.coeff(i);
      }

      value_type
      computeFastAccessCoeff(int i, const ExprT1& expr1,
                             const ExprT2& expr2) const {
        if (i == 0)
          return std::min(expr1.value(), expr2.fastAccessCoeff(0));
        else
          return expr1.value() <= expr2.fastAccessCoeff(0) ? value_type(0) :
            expr2.fastAccessCoeff(i);
      }

    }; // class MinOp

    // ------------------------ quadrature function ---------------------------

    template <typename ExprT1, typename ExprT2>
    class ASinQuadOp {
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

      ASinQuadOp(const ExprT1& expr1, const ExprT2& expr2) :
        c(),
        dc(-1)
      {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT1& expr1,
                              const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::asin(expr1.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*expr2.coeff(k-j)*expr1.coeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT1& expr1,
                                        const ExprT2& expr2) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::asin(expr1.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*expr2.fastAccessCoeff(k-j)*
                expr1.fastAccessCoeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ASinQuadOp

    template <typename ExprT1, typename ExprT2>
    class ACosQuadOp {
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

      ACosQuadOp(const ExprT1& expr1, const ExprT2& expr2) :
        c(),
        dc(-1)
      {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT1& expr1,
                              const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::acos(expr1.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*expr2.coeff(k-j)*expr1.coeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT1& expr1,
                                        const ExprT2& expr2) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::acos(expr1.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*expr2.fastAccessCoeff(k-j)*
                expr1.fastAccessCoeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ACosQuadOp

    template <typename ExprT1, typename ExprT2>
    class ATanQuadOp {
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

      ATanQuadOp(const ExprT1& expr1, const ExprT2& expr2) :
        c(),
        dc(-1)
      {}

      void allocateCache(int d) const {
        c.resize(d+1,value_type(0));
      }

      value_type computeCoeff(int i, const ExprT1& expr1,
                              const ExprT2& expr2) const {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::atan(expr1.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*expr2.coeff(k-j)*expr1.coeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

      value_type computeFastAccessCoeff(int i,
                                        const ExprT1& expr1,
                                        const ExprT2& expr2) const
      {
        if (static_cast<int>(i) > dc) {
          if (dc < 0) {
            c[0] = std::atan(expr1.fastAccessCoeff(0));
            dc = 0;
          }
          for (int k=dc+1; k<=i; k++) {
            for (int j=1; j<=k; j++)
              c[k] += value_type(j)*expr2.fastAccessCoeff(k-j)*
                expr1.fastAccessCoeff(j);
            c[k] /= value_type(k);
          }
          dc = i;
        }
        return c[i];
      }

    protected:

      mutable std::valarray<value_type> c;
      mutable int dc;

    }; // class ATanQuadOp

  } // namespace Tay

} // namespace Sacado

#define TAYLOR_BINARYOP_MACRO(OPNAME,OP)                                \
namespace Sacado {                                                      \
  namespace Tay {                                                       \
    template <typename T1, typename T2>                                 \
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, OP > >                 \
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)               \
    {                                                                   \
      typedef BinaryExpr< Expr<T1>, Expr<T2>, OP > expr_t;              \
                                                                        \
      return Expr<expr_t>(expr_t(expr1, expr2));                        \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    inline Expr< BinaryExpr< ConstExpr<typename Expr<T>::value_type>,   \
                             Expr<T>, OP > >                            \
    OPNAME (const typename Expr<T>::value_type& c,                      \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;           \
      typedef BinaryExpr< ConstT, Expr<T>, OP > expr_t;                 \
                                                                        \
      return Expr<expr_t>(expr_t(ConstT(c), expr));                     \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    inline Expr< BinaryExpr< Expr<T>,                                   \
                             ConstExpr<typename Expr<T>::value_type>,   \
                             OP > >                                     \
    OPNAME (const Expr<T>& expr,                                        \
            const typename Expr<T>::value_type& c)                      \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;           \
      typedef BinaryExpr< Expr<T>, ConstT, OP > expr_t;                 \
                                                                        \
      return Expr<expr_t>(expr_t(expr, ConstT(c)));                     \
    }                                                                   \
  }                                                                     \
}

TAYLOR_BINARYOP_MACRO(operator+, AdditionOp)
TAYLOR_BINARYOP_MACRO(operator-, SubtractionOp)
TAYLOR_BINARYOP_MACRO(operator*, MultiplicationOp)
TAYLOR_BINARYOP_MACRO(operator/, DivisionOp)

#undef TAYLOR_BINARYOP_MACRO

  // The general definition of max/min works for Taylor variables too, except
  // we need to add a case when the argument types are different.  This
  // can't conflict with the general definition, so we need to use
  // Substitution Failure Is Not An Error
#include "Sacado_mpl_disable_if.hpp"

#define TAYLOR_SFINAE_BINARYOP_MACRO(OPNAME,OP)                         \
namespace Sacado {                                                      \
  namespace Tay {                                                       \
    template <typename T1, typename T2>                                 \
    inline                                                              \
    typename                                                            \
    mpl::disable_if< std::is_same<T1,T2>,                               \
                     Expr<BinaryExpr<Expr<T1>, Expr<T2>, OP> > >::type  \
    OPNAME (const Expr<T1>& expr1, const Expr<T2>& expr2)               \
    {                                                                   \
      typedef BinaryExpr< Expr<T1>, Expr<T2>, OP > expr_t;              \
                                                                        \
      return Expr<expr_t>(expr_t(expr1, expr2));                        \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    inline Expr< BinaryExpr< ConstExpr<typename Expr<T>::value_type>,   \
                             Expr<T>, OP > >                            \
    OPNAME (const typename Expr<T>::value_type& c,                      \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;           \
      typedef BinaryExpr< ConstT, Expr<T>, OP > expr_t;                 \
                                                                        \
      return Expr<expr_t>(expr_t(ConstT(c), expr));                     \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    inline Expr< BinaryExpr< Expr<T>,                                   \
                             ConstExpr<typename Expr<T>::value_type>,   \
                             OP > >                                     \
    OPNAME (const Expr<T>& expr,                                        \
            const typename Expr<T>::value_type& c)                      \
    {                                                                   \
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;           \
      typedef BinaryExpr< Expr<T>, ConstT, OP > expr_t;                 \
                                                                        \
      return Expr<expr_t>(expr_t(expr, ConstT(c)));                     \
    }                                                                   \
  }                                                                     \
}

TAYLOR_SFINAE_BINARYOP_MACRO(max, MaxOp)
TAYLOR_SFINAE_BINARYOP_MACRO(min, MinOp)

#undef TAYLOR_SFINAE_BINARYOP_MACRO

namespace std {
  using Sacado::Tay::min;
  using Sacado::Tay::max;
}

namespace Sacado {

  namespace Tay {

    template <typename T1, typename T2>
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, ASinQuadOp > >
    asin_quad (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef BinaryExpr< Expr<T1>, Expr<T2>, ASinQuadOp > expr_t;

      return Expr<expr_t>(expr_t(expr1, expr2));
    }

    template <typename T1, typename T2>
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, ACosQuadOp > >
    acos_quad (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef BinaryExpr< Expr<T1>, Expr<T2>, ACosQuadOp > expr_t;

      return Expr<expr_t>(expr_t(expr1, expr2));
    }

    template <typename T1, typename T2>
    inline Expr< BinaryExpr< Expr<T1>, Expr<T2>, ATanQuadOp > >
    atan_quad (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      typedef BinaryExpr< Expr<T1>, Expr<T2>, ATanQuadOp > expr_t;

      return Expr<expr_t>(expr_t(expr1, expr2));
    }

    template <typename ExprT1, typename ExprT2>
    struct PowExprType {
      typedef UnaryExpr< ExprT1, LogOp > T3;
      typedef BinaryExpr< ExprT2, Expr<T3>, MultiplicationOp > T4;
      typedef UnaryExpr< Expr<T4>, ExpOp > T5;

      typedef Expr<T5> expr_type;
    };

     template <typename ExprT2>
     struct PowExprType< typename ExprT2::value_type, ExprT2 > {
       typedef typename ExprT2::value_type T1;
      typedef BinaryExpr< ExprT2, ConstExpr<T1>, MultiplicationOp > T4;
      typedef UnaryExpr< Expr<T4>, ExpOp > T5;

      typedef Expr<T5> expr_type;
    };

    template <typename ExprT1>
    struct PowExprType<ExprT1,typename ExprT1::value_type> {
      typedef typename ExprT1::value_type T2;
      typedef UnaryExpr< ExprT1, LogOp > T3;
      typedef BinaryExpr< ConstExpr<T2>, Expr<T3>, MultiplicationOp > T4;
      typedef UnaryExpr< Expr<T4>, ExpOp > T5;

      typedef Expr<T5> expr_type;
    };

    template <typename T1, typename T2>
    inline typename PowExprType< Expr<T1>, Expr<T2> >::expr_type
    pow (const Expr<T1>& expr1, const Expr<T2>& expr2)
    {
      // pow(x,y) = exp(y*log(x))
      return exp(expr2*log(expr1));
    }

    template <typename T>
    inline typename PowExprType< typename Expr<T>::value_type, Expr<T> >::expr_type
    pow (const typename Expr<T>::value_type& c, const Expr<T>& expr)
    {
      // pow(x,y) = exp(y*log(x))
      return exp(expr*std::log(c));
    }

    template <typename T>
    inline typename PowExprType< Expr<T>, typename Expr<T>::value_type >::expr_type
    pow (const Expr<T>& expr, const typename Expr<T>::value_type& c)
    {
      // pow(x,y) = exp(y*log(x))
      return exp(c*log(expr));
    }

    template <typename T>
    struct Log10ExprType {
      typedef UnaryExpr< Expr<T>, LogOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< Expr<T1>, ConstT, DivisionOp > T2;
      typedef Expr<T2> expr_type;
    };

    template <typename T>
    inline typename Log10ExprType<T>::expr_type
    log10 (const Expr<T>& expr)
    {
      // log10(x) = log(x)/log(10)
      return log(expr)/std::log(typename Expr<T>::value_type(10));
    }

    template <typename T>
    struct TanExprType {
      typedef UnaryExpr< Expr<T>, SinOp > T1;
      typedef UnaryExpr< Expr<T>, CosOp > T2;
      typedef BinaryExpr< Expr<T1>, Expr<T2>, DivisionOp > T3;
      typedef Expr<T3> expr_type;
    };

    template <typename T>
    inline typename TanExprType<T>::expr_type
    tan (const Expr<T>& expr)
    {
      // tan(x) = sin(x)/cos(x)
      return sin(expr)/cos(expr);
    }

    template <typename T>
    struct ASinExprType {
      typedef BinaryExpr< Expr<T>, Expr<T>, MultiplicationOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< ConstT, Expr<T1>, SubtractionOp > T2;
      typedef UnaryExpr< Expr<T2>, SqrtOp > T3;
      typedef BinaryExpr< ConstT, Expr<T3>, DivisionOp > T4;
      typedef BinaryExpr< Expr<T>, Expr<T4>, ASinQuadOp > T5;
      typedef Expr<T5> expr_type;
    };

    template <typename T>
    inline typename ASinExprType<T>::expr_type
    asin (const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type value_type;

      // asin(x) = integral of 1/sqrt(1-x*x)
      return asin_quad(expr, value_type(1)/sqrt(value_type(1)-expr*expr));
    }

    template <typename T>
    struct ACosExprType {
      typedef BinaryExpr< Expr<T>, Expr<T>, MultiplicationOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< ConstT, Expr<T1>, SubtractionOp > T2;
      typedef UnaryExpr< Expr<T2>, SqrtOp > T3;
      typedef BinaryExpr< ConstT, Expr<T3>, DivisionOp > T4;
      typedef BinaryExpr< Expr<T>, Expr<T4>, ACosQuadOp > T5;
      typedef Expr<T5> expr_type;
    };

    template <typename T>
    inline typename ACosExprType<T>::expr_type
    acos (const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type value_type;

      // acos(x) = integral of -1/sqrt(1-x*x)
      return acos_quad(expr, value_type(-1)/sqrt(value_type(1)-expr*expr));
    }

    template <typename T>
    struct ATanExprType {
      typedef BinaryExpr< Expr<T>, Expr<T>, MultiplicationOp > T1;
      typedef ConstExpr<typename Expr<T>::value_type> ConstT;
      typedef BinaryExpr< ConstT, Expr<T1>, AdditionOp > T2;
      typedef BinaryExpr< ConstT, Expr<T2>, DivisionOp > T3;
      typedef BinaryExpr< Expr<T>, Expr<T3>, ATanQuadOp > T4;
      typedef Expr<T4> expr_type;
    };

    template <typename T>
    inline typename ATanExprType<T>::expr_type
    atan (const Expr<T>& expr)
    {
      typedef typename Expr<T>::value_type value_type;

      // atan(x) = integral of 1/(1+x*x)
      return atan_quad(expr, value_type(1)/(value_type(1)+expr*expr));
    }

    template <typename T>
    struct TanhExprType {
      typedef UnaryExpr< Expr<T>, SinhOp > T1;
      typedef UnaryExpr< Expr<T>, CoshOp > T2;
      typedef BinaryExpr< Expr<T1>, Expr<T2>, DivisionOp > T3;
      typedef Expr<T3> expr_type;
    };

    template <typename T>
    inline typename TanhExprType<T>::expr_type
    tanh (const Expr<T>& expr)
    {
      // tanh(x) = sinh(x)/cosh(x)
      return sinh(expr)/cosh(expr);
    }

  } // namespace Tay

} // namespace Sacado

namespace std {
  using Sacado::Tay::pow;
  using Sacado::Tay::log10;
  using Sacado::Tay::tan;
  using Sacado::Tay::asin;
  using Sacado::Tay::acos;
  using Sacado::Tay::atan;
  using Sacado::Tay::tanh;
}

//-------------------------- Relational Operators -----------------------

#define TAYLOR_RELOP_MACRO(OP)                                          \
namespace Sacado {                                                      \
  namespace Tay {                                                       \
    template <typename ExprT1, typename ExprT2>                         \
    inline bool                                                         \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return expr1.fastAccessCoeff(0) OP expr2.fastAccessCoeff(0);      \
    }                                                                   \
                                                                        \
    template <typename ExprT2>                                          \
    inline bool                                                         \
    operator OP (const typename Expr<ExprT2>::value_type& a,            \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return a OP expr2.fastAccessCoeff(0);                             \
    }                                                                   \
                                                                        \
    template <typename ExprT1>                                          \
    inline bool                                                         \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const typename Expr<ExprT1>::value_type& b)            \
    {                                                                   \
      return expr1.fastAccessCoeff(0) OP b;                             \
    }                                                                   \
  }                                                                     \
}

TAYLOR_RELOP_MACRO(==)
TAYLOR_RELOP_MACRO(!=)
TAYLOR_RELOP_MACRO(<)
TAYLOR_RELOP_MACRO(>)
TAYLOR_RELOP_MACRO(<=)
TAYLOR_RELOP_MACRO(>=)
TAYLOR_RELOP_MACRO(<<=)
TAYLOR_RELOP_MACRO(>>=)
TAYLOR_RELOP_MACRO(&)
TAYLOR_RELOP_MACRO(|)

#undef TAYLOR_RELOP_MACRO

namespace Sacado {

  namespace Tay {

    template <typename ExprT>
    inline bool operator ! (const Expr<ExprT>& expr)
    {
      return ! expr.fastAccessCoeff(0);
    }

  } // namespace Tay

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace Tay {

    template <typename ExprT>
    bool toBool2(const Expr<ExprT>& x) {
      bool is_zero = true;
      for (int i=0; i<=x.degree(); i++)
        is_zero = is_zero && (x.coeff(i) == 0.0);
      return !is_zero;
    }

  } // namespace Tay

} // namespace Sacado

#define TAY_BOOL_MACRO(OP)                                              \
namespace Sacado {                                                      \
  namespace Tay {                                                       \
    template <typename ExprT1, typename ExprT2>                         \
    inline bool                                                         \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return toBool2(expr1) OP toBool2(expr2);                          \
    }                                                                   \
                                                                        \
    template <typename ExprT2>                                          \
    inline bool                                                         \
    operator OP (const typename Expr<ExprT2>::value_type& a,            \
                 const Expr<ExprT2>& expr2)                             \
    {                                                                   \
      return a OP toBool2(expr2);                                       \
    }                                                                   \
                                                                        \
    template <typename ExprT1>                                          \
    inline bool                                                         \
    operator OP (const Expr<ExprT1>& expr1,                             \
                 const typename Expr<ExprT1>::value_type& b)            \
    {                                                                   \
      return toBool2(expr1) OP b;                                       \
    }                                                                   \
  }                                                                     \
}

TAY_BOOL_MACRO(&&)
TAY_BOOL_MACRO(||)

#undef TAY_BOOL_MACRO

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace Tay {

    template <typename ExprT>
    std::ostream& operator << (std::ostream& os, const Expr<ExprT>& x) {
      os.setf(std::ios::fixed, std::ios::floatfield);
      os.width(12);
      os << "[";

      for (int i=0; i<=x.degree(); i++) {
        os.width(12);
        os << x.coeff(i);
      }

      os << "]";
      return os;
    }

    //! Compute Taylor series of n-th derivative of x
    template <typename T>
    CacheTaylor<T> diff(const CacheTaylor<T>& x, int n = 1) {
      const int d = x.degree();
      if (n <= 0)
        return x;
      else if (n > d) {
        CacheTaylor<T> y(0);
        return y;
      }
      CacheTaylor<T> y(d-n, 0);
      int c = 1;
      for (int i=1; i<=n; ++i)
        c *= i;
      for (int i=n; i<=d; ++i) {
        y.fastAccessCoeff(i-n) = x.fastAccessCoeff(i) * T(c);
        c = (c / (i-n+1)) * (i+1);
      }
      return y;
    }

  } // namespace Tay

} // namespace Sacado


#endif // SACADO_TAY_CACHETAYLOROPS_HPP
