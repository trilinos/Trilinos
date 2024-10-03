// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_OPS_HPP
#define SACADO_FAD_EXP_OPS_HPP

#include <type_traits>
#include <ostream>

#include "Sacado_Fad_Exp_Expression.hpp"
#include "Sacado_Fad_Exp_ExpressionTraits.hpp"
#include "Sacado_Fad_Exp_Ops_Fwd.hpp"
#include "Sacado_cmath.hpp"

#include "Sacado_mpl_has_equal_to.hpp"

#if defined(HAVE_SACADO_KOKKOS)
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Error.hpp"
#endif

#define FAD_UNARYOP_MACRO(OPNAME,OP,USING,VALUE,DX,FASTACCESSDX)        \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
                                                                        \
    template <typename T, typename ExprSpec>                            \
    class OP {};                                                        \
                                                                        \
    template <typename T>                                               \
    class OP< T,ExprSpecDefault > :                                     \
      public Expr< OP<T,ExprSpecDefault> > {                            \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T>::type ExprT;                   \
      typedef typename ExprT::value_type value_type;                    \
      typedef typename ExprT::scalar_type scalar_type;                  \
                                                                        \
      typedef ExprSpecDefault expr_spec_type;                           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      explicit OP(const T& expr_) : expr(expr_)  {}                     \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const { return expr.size(); }                          \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr.hasFastAccess();                                    \
      }                                                                 \
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
      const T& expr;                                                    \
    };                                                                  \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    OP< typename Expr<T>::derived_type,                                 \
        typename T::expr_spec_type >                                    \
    OPNAME (const Expr<T>& expr)                                        \
    {                                                                   \
      typedef OP< typename Expr<T>::derived_type,                       \
                  typename T::expr_spec_type > expr_t;                  \
                                                                        \
      return expr_t(expr.derived());                                    \
    }                                                                   \
                                                                        \
    template <typename T, typename E>                                   \
    struct ExprLevel< OP< T,E > > {                                     \
      static const unsigned value = ExprLevel<T>::value;                \
    };                                                                  \
                                                                        \
    template <typename T, typename E>                                   \
    struct IsFadExpr< OP< T,E > > {                                     \
      static const unsigned value = true;                               \
    };                                                                  \
                                                                        \
  }                                                                     \
  }                                                                     \
                                                                        \
  template <typename T, typename E>                                     \
  struct IsExpr< Fad::Exp::OP< T,E > > {                                \
    static const bool value = true;                                     \
  };                                                                    \
                                                                        \
  template <typename T, typename E>                                     \
  struct BaseExprType< Fad::Exp::OP< T,E > > {                          \
    typedef typename BaseExprType<T>::type type;                        \
  };                                                                    \
                                                                        \
  template <typename T, typename E>                                     \
  struct IsSimdType< Fad::Exp::OP< T,E > > {                            \
    static const bool value =                                           \
      IsSimdType< typename Fad::Exp::OP< T,E >::scalar_type >::value;   \
  };                                                                    \
                                                                        \
  template <typename T, typename E>                                     \
  struct ValueType< Fad::Exp::OP< T,E > > {                             \
    typedef typename Fad::Exp::OP< T,E >::value_type type;              \
  };                                                                    \
                                                                        \
  template <typename T, typename E>                                     \
  struct ScalarType< Fad::Exp::OP< T,E > > {                            \
    typedef typename Fad::Exp::OP< T,E >::scalar_type type;             \
  };                                                                    \
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
                  tan(expr.val()),
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

// Special handling for safe_sqrt() to provide specializations of SafeSqrtOp for
// "simd" value types that use if_then_else(). The only reason for not using
// if_then_else() always is to avoid evaluating the derivative if the value is
// zero to avoid throwing FPEs.
namespace Sacado {
  namespace Fad {
  namespace Exp {

    template <typename T, typename ExprSpec, bool is_simd>
    class SafeSqrtOp {};

    //
    // Implementation for simd type using if_then_else()
    //
    template <typename T>
    class SafeSqrtOp< T,ExprSpecDefault,true > :
      public Expr< SafeSqrtOp<T,ExprSpecDefault> > {
    public:

      typedef typename std::remove_cv<T>::type ExprT;
      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      explicit SafeSqrtOp(const T& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

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

      const T& expr;
    };

    //
    // Specialization for scalar types using ternary operator
    //
    template <typename T>
    class SafeSqrtOp< T,ExprSpecDefault,false > :
      public Expr< SafeSqrtOp<T,ExprSpecDefault> > {
    public:

      typedef typename std::remove_cv<T>::type ExprT;
      typedef typename ExprT::value_type value_type;
      typedef typename ExprT::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      explicit SafeSqrtOp(const T& expr_) : expr(expr_)  {}

      SACADO_INLINE_FUNCTION
      int size() const { return expr.size(); }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr.hasFastAccess();
      }

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

      const T& expr;
    };

    template <typename T>
    SACADO_INLINE_FUNCTION
    SafeSqrtOp< typename Expr<T>::derived_type,
                typename T::expr_spec_type >
    safe_sqrt (const Expr<T>& expr)
    {
      typedef SafeSqrtOp< typename Expr<T>::derived_type,
                          typename T::expr_spec_type > expr_t;

      return expr_t(expr.derived());
    }

    template <typename T, typename E>
    struct ExprLevel< SafeSqrtOp< T,E > > {
      static const unsigned value = ExprLevel<T>::value;
    };

    template <typename T, typename E>
    struct IsFadExpr< SafeSqrtOp< T,E > > {
      static const unsigned value = true;
    };

  }
  }

  template <typename T, typename E>
  struct IsExpr< Fad::Exp::SafeSqrtOp< T,E > > {
    static const bool value = true;
  };

  template <typename T, typename E>
  struct BaseExprType< Fad::Exp::SafeSqrtOp< T,E > > {
    typedef typename BaseExprType<T>::type type;
  };

  template <typename T, typename E>
  struct IsSimdType< Fad::Exp::SafeSqrtOp< T,E > > {
    static const bool value =
      IsSimdType< typename Fad::Exp::SafeSqrtOp< T,E >::scalar_type >::value;
  };

  template <typename T, typename E>
  struct ValueType< Fad::Exp::SafeSqrtOp< T,E > > {
    typedef typename Fad::Exp::SafeSqrtOp< T,E >::value_type type;
  };

  template <typename T, typename E>
  struct ScalarType< Fad::Exp::SafeSqrtOp< T,E > > {
    typedef typename Fad::Exp::SafeSqrtOp< T,E >::scalar_type type;
  };

}

#undef FAD_UNARYOP_MACRO

#define FAD_BINARYOP_MACRO(OPNAME,OP,USING,VALUE,DX,CDX1,CDX2,FASTACCESSDX,VAL_CONST_DX_1,VAL_CONST_DX_2,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
                                                                        \
    template <typename T1, typename T2,                                 \
              bool is_const_T1, bool is_const_T2,                       \
              typename ExprSpec >                                       \
    class OP {};                                                        \
                                                                        \
    template <typename T1, typename T2>                                 \
    class OP< T1, T2, false, false, ExprSpecDefault > :                 \
      public Expr< OP< T1, T2, false, false, ExprSpecDefault > > {      \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T1>::type ExprT1;                 \
      typedef typename std::remove_cv<T2>::type ExprT2;                 \
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
      typedef ExprSpecDefault expr_spec_type;                           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const T2& expr2_) :                          \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      int size() const {                                                \
        const int sz1 = expr1.size(), sz2 = expr2.size();               \
        return sz1 > sz2 ? sz1 : sz2;                                   \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess() && expr2.hasFastAccess();          \
      }                                                                 \
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
        const int sz1 = expr1.size(), sz2 = expr2.size();               \
        if (sz1 > 0 && sz2 > 0)                                         \
          return DX;                                                    \
        else if (sz1 > 0)                                               \
          return CDX2;                                                  \
        else                                                            \
          return CDX1;                                                  \
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
      const T1& expr1;                                                  \
      const T2& expr2;                                                  \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    class OP< T1, T2, false, true, ExprSpecDefault >                    \
      : public Expr< OP< T1, T2, false, true, ExprSpecDefault > > {     \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T1>::type ExprT1;                 \
      typedef T2 ConstT;                                                \
      typedef typename ExprT1::value_type value_type;                   \
      typedef typename ExprT1::scalar_type scalar_type;                 \
                                                                        \
      typedef ExprSpecDefault expr_spec_type;                           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const ConstT& c_) :                          \
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
      value_type val() const {                                          \
        USING                                                           \
        return VAL_CONST_DX_2;                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type dx(int i) const {                                      \
        USING                                                           \
        return CONST_DX_2;                                              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type fastAccessDx(int i) const {                            \
        USING                                                           \
        return CONST_FASTACCESSDX_2;                                    \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const T1& expr1;                                                  \
      const ConstT& c;                                                  \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    class OP< T1, T2, true, false, ExprSpecDefault >                    \
      : public Expr< OP< T1, T2, true, false, ExprSpecDefault > > {     \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T2>::type ExprT2;                 \
      typedef T1 ConstT;                                                \
      typedef typename ExprT2::value_type value_type;                   \
      typedef typename ExprT2::scalar_type scalar_type;                 \
                                                                        \
      typedef ExprSpecDefault expr_spec_type;                           \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      OP(const ConstT& c_, const T2& expr2_) :                          \
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
      value_type val() const {                                          \
        USING                                                           \
        return VAL_CONST_DX_1;                                          \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type dx(int i) const {                                      \
        USING                                                           \
        return CONST_DX_1;                                              \
      }                                                                 \
                                                                        \
      SACADO_INLINE_FUNCTION                                            \
      value_type fastAccessDx(int i) const {                            \
        USING                                                           \
        return CONST_FASTACCESSDX_1;                                    \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ConstT& c;                                                  \
      const T2& expr2;                                                  \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_EXPR_EXPR(OP)                              \
    OPNAME (const T1& expr1, const T2& expr2)                           \
    {                                                                   \
      typedef OP< typename Expr<T1>::derived_type,                      \
                  typename Expr<T2>::derived_type,                      \
                  false, false, typename T1::expr_spec_type > expr_t;   \
                                                                        \
      return expr_t(expr1.derived(), expr2.derived());                  \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    OP< typename T::value_type, typename Expr<T>::derived_type,         \
        true, false, typename T::expr_spec_type >                       \
    OPNAME (const typename T::value_type& c,                            \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef OP< ConstT, typename Expr<T>::derived_type,               \
                  true, false, typename T::expr_spec_type > expr_t;     \
                                                                        \
      return expr_t(c, expr.derived());                                 \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    OP< typename Expr<T>::derived_type, typename T::value_type,         \
        false, true, typename T::expr_spec_type >                       \
    OPNAME (const Expr<T>& expr,                                        \
            const typename T::value_type& c)                            \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef OP< typename Expr<T>::derived_type, ConstT,               \
                  false, true, typename T::expr_spec_type > expr_t;     \
                                                                        \
      return expr_t(expr.derived(), c);                                 \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_SCALAR_EXPR(OP)                            \
    OPNAME (const typename T::scalar_type& c,                           \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef typename T::scalar_type ConstT;                           \
      typedef OP< ConstT, typename Expr<T>::derived_type,               \
                  true, false, typename T::expr_spec_type > expr_t;     \
                                                                        \
      return expr_t(c, expr.derived());                                 \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_EXPR_SCALAR(OP)                            \
    OPNAME (const Expr<T>& expr,                                        \
            const typename T::scalar_type& c)                           \
    {                                                                   \
      typedef typename T::scalar_type ConstT;                           \
      typedef OP< typename Expr<T>::derived_type, ConstT,               \
                  false, true, typename T::expr_spec_type > expr_t;     \
                                                                        \
      return expr_t(expr.derived(), c);                                 \
    }                                                                   \
                                                                        \
    template <typename T1, typename T2, bool c1, bool c2, typename E>   \
    struct ExprLevel< OP< T1, T2, c1, c2, E > > {                       \
      static constexpr unsigned value_1 = ExprLevel<T1>::value;         \
      static constexpr unsigned value_2 = ExprLevel<T2>::value;         \
      static constexpr unsigned value =                                 \
        value_1 >= value_2 ? value_1 : value_2;                         \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2, bool c1, bool c2, typename E>   \
    struct IsFadExpr< OP< T1, T2, c1, c2, E > > {                       \
      static constexpr unsigned value = true;                           \
    };                                                                  \
                                                                        \
  }                                                                     \
  }                                                                     \
                                                                        \
  template <typename T1, typename T2, bool c1, bool c2, typename E>     \
  struct IsExpr< Fad::Exp::OP< T1, T2, c1, c2, E > > {                  \
    static constexpr bool value = true;                                 \
  };                                                                    \
                                                                        \
  template <typename T1, typename T2, bool c1, bool c2, typename E>     \
  struct BaseExprType< Fad::Exp::OP< T1, T2, c1, c2, E > > {            \
    typedef typename BaseExprType<T1>::type base_expr_1;                \
    typedef typename BaseExprType<T2>::type base_expr_2;                \
    typedef typename Sacado::Promote<base_expr_1,                       \
                                     base_expr_2>::type type;           \
  };                                                                    \
                                                                        \
  template <typename T1, typename T2, bool c1, bool c2, typename E>     \
  struct IsSimdType< Fad::Exp::OP< T1, T2, c1, c2, E > > {              \
    static const bool value =                                           \
      IsSimdType< typename Fad::Exp::OP< T1, T2, c1, c2, E >::value_type >::value; \
  };                                                                    \
                                                                        \
  template <typename T1, typename T2, bool c1, bool c2, typename E>     \
  struct ValueType< Fad::Exp::OP< T1, T2, c1, c2, E > > {               \
    typedef typename Fad::Exp::OP< T1, T2, c1, c2, E >::value_type type;\
  };                                                                    \
                                                                        \
  template <typename T1, typename T2, bool c1, bool c2, typename E>     \
  struct ScalarType< Fad::Exp::OP< T1, T2, c1, c2, E > > {              \
    typedef typename Fad::Exp::OP< T1, T2, c1, c2, E >::scalar_type type;\
  };                                                                    \
                                                                        \
}


FAD_BINARYOP_MACRO(operator+,
                   AdditionOp,
                   ;,
                   expr1.val() + expr2.val(),
                   expr1.dx(i) + expr2.dx(i),
                   expr2.dx(i),
                   expr1.dx(i),
                   expr1.fastAccessDx(i) + expr2.fastAccessDx(i),
                   c + expr2.val(),
                   expr1.val() + c,
                   expr2.dx(i),
                   expr1.dx(i),
                   expr2.fastAccessDx(i),
                   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator-,
                   SubtractionOp,
                   ;,
                   expr1.val() - expr2.val(),
                   expr1.dx(i) - expr2.dx(i),
                   -expr2.dx(i),
                   expr1.dx(i),
                   expr1.fastAccessDx(i) - expr2.fastAccessDx(i),
                   c - expr2.val(),
                   expr1.val() - c,
                   -expr2.dx(i),
                   expr1.dx(i),
                   -expr2.fastAccessDx(i),
                   expr1.fastAccessDx(i))
FAD_BINARYOP_MACRO(operator*,
                   MultiplicationOp,
                   ;,
                   expr1.val() * expr2.val(),
                   expr1.val()*expr2.dx(i) + expr1.dx(i)*expr2.val(),
                   expr1.val()*expr2.dx(i),
                   expr1.dx(i)*expr2.val(),
                   expr1.val()*expr2.fastAccessDx(i) +
                     expr1.fastAccessDx(i)*expr2.val(),
                   c * expr2.val(),
                   expr1.val() * c,
                   c*expr2.dx(i),
                   expr1.dx(i)*c,
                   c*expr2.fastAccessDx(i),
                   expr1.fastAccessDx(i)*c)
FAD_BINARYOP_MACRO(operator/,
                   DivisionOp,
                   ;,
                   expr1.val() / expr2.val(),
                   (expr1.dx(i)*expr2.val() - expr2.dx(i)*expr1.val()) /
                     (expr2.val()*expr2.val()),
                   -expr2.dx(i)*expr1.val() / (expr2.val()*expr2.val()),
                   expr1.dx(i)/expr2.val(),
                   (expr1.fastAccessDx(i)*expr2.val() -
                      expr2.fastAccessDx(i)*expr1.val()) /
                      (expr2.val()*expr2.val()),
                   c / expr2.val(),
                   expr1.val() / c,
                   -expr2.dx(i)*c / (expr2.val()*expr2.val()),
                   expr1.dx(i)/c,
                   -expr2.fastAccessDx(i)*c / (expr2.val()*expr2.val()),
                   expr1.fastAccessDx(i)/c)
FAD_BINARYOP_MACRO(atan2,
                   Atan2Op,
                   using std::atan2;,
                   atan2(expr1.val(), expr2.val()),
                   (expr2.val()*expr1.dx(i) - expr1.val()*expr2.dx(i))/
                        (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   -expr1.val()*expr2.dx(i)/
                        (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   expr2.val()*expr1.dx(i)/
                        (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   (expr2.val()*expr1.fastAccessDx(i) - expr1.val()*expr2.fastAccessDx(i))/
                        (expr1.val()*expr1.val() + expr2.val()*expr2.val()),
                   atan2(c, expr2.val()),
                   atan2(expr1.val(), c),
                   (-c*expr2.dx(i)) / (c*c + expr2.val()*expr2.val()),
                   (c*expr1.dx(i))/ (expr1.val()*expr1.val() + c*c),
                   (-c*expr2.fastAccessDx(i))/ (c*c + expr2.val()*expr2.val()),
                   (c*expr1.fastAccessDx(i))/ (expr1.val()*expr1.val() + c*c))
// FAD_BINARYOP_MACRO(pow,
//                    PowerOp,
//                    using std::pow; using std::log; using Sacado::if_then_else;,
//                    pow(expr1.val(), expr2.val()),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type((expr2.dx(i)*log(expr1.val())+expr2.val()*expr1.dx(i)/expr1.val())*pow(expr1.val(),expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type(expr2.dx(i)*log(expr1.val())*pow(expr1.val(),expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type(expr2.val()*expr1.dx(i)/expr1.val()*pow(expr1.val(),expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type((expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val())) ),
//                    pow(c, expr2.val()),
//                    pow(expr1.val(), c),
//                    if_then_else( c == value_type(0.0), value_type(0.0), value_type(expr2.dx(i)*log(c)*pow(c,expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type(c*expr1.dx(i)/expr1.val()*pow(expr1.val(),c)) ),
//                    if_then_else( c == value_type(0.0), value_type(0.0), value_type(expr2.fastAccessDx(i)*log(c)*pow(c,expr2.val())) ),
//                    if_then_else( expr1.val() == value_type(0.0), value_type(0.0), value_type(c*expr1.fastAccessDx(i)/expr1.val()*pow(expr1.val(),c))) )
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   using Sacado::if_then_else;,
                   if_then_else( expr1.val() >= expr2.val(),  expr1.val(), expr2.val() ),
                   if_then_else( expr1.val() >= expr2.val(), expr1.dx(i), expr2.dx(i) ),
                   if_then_else( expr1.val() >= expr2.val(), value_type(0.0), expr2.dx(i) ),
                   if_then_else( expr1.val() >= expr2.val(), expr1.dx(i), value_type(0.0) ),
                   if_then_else( expr1.val() >= expr2.val(), expr1.fastAccessDx(i), expr2.fastAccessDx(i) ),
                   if_then_else( c >= expr2.val(), value_type(c),  expr2.val() ),
                   if_then_else( expr1.val() >= c, expr1.val(), value_type(c) ),
                   if_then_else( c >= expr2.val(), value_type(0.0),  expr2.dx(i) ),
                   if_then_else( expr1.val() >= c, expr1.dx(i), value_type(0.0) ),
                   if_then_else( c >= expr2.val(), value_type(0.0), expr2.fastAccessDx(i) ),
                   if_then_else( expr1.val() >= c, expr1.fastAccessDx(i), value_type(0.0) ) )
FAD_BINARYOP_MACRO(min,
                   MinOp,
                   using Sacado::if_then_else;,
                   if_then_else( expr1.val() <= expr2.val(), expr1.val(), expr2.val() ),
                   if_then_else( expr1.val() <= expr2.val(), expr1.dx(i), expr2.dx(i) ),
                   if_then_else( expr1.val() <= expr2.val(), value_type(0.0), expr2.dx(i) ),
                   if_then_else( expr1.val() <= expr2.val(), expr1.dx(i), value_type(0.0) ),
                   if_then_else( expr1.val() <= expr2.val(), expr1.fastAccessDx(i), expr2.fastAccessDx(i) ),
                   if_then_else( c <= expr2.val(), value_type(c), expr2.val() ),
                   if_then_else( expr1.val() <= c, expr1.val(), value_type(c) ),
                   if_then_else( c <= expr2.val(), value_type(0), expr2.dx(i) ),
                   if_then_else( expr1.val() <= c, expr1.dx(i), value_type(0) ),
                   if_then_else( c <= expr2.val(), value_type(0), expr2.fastAccessDx(i) ),
                   if_then_else( expr1.val() <= c, expr1.fastAccessDx(i), value_type(0) ) )


// Special handling for std::pow() to provide specializations of PowerOp for
// "simd" value types that use if_then_else(). The only reason for not using
// if_then_else() always is to avoid evaluating the derivative if the value is
// zero to avoid throwing FPEs.
namespace Sacado {
  namespace Fad {
  namespace Exp {

    template <typename T1, typename T2,
              bool is_const_T1, bool is_const_T2,
              typename ExprSpec, typename Impl >
    class PowerOp {};

    //
    // Implementation for simd type using if_then_else()
    //
    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, false, ExprSpecDefault,
                   PowerImpl::Simd > :
      public Expr< PowerOp< T1, T2, false, false, ExprSpecDefault,
                   PowerImpl::Simd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const T2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        const int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
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
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        return if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type((expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val())));
      }

    protected:

      const T1& expr1;
      const T2& expr2;

    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, true, ExprSpecDefault,
                   PowerImpl::Simd > :
      public Expr< PowerOp< T1, T2, false, true, ExprSpecDefault,
                   PowerImpl::Simd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const ConstT& c_) :
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
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), c);
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow; using Sacado::if_then_else;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return if_then_else( c == scalar_type(1.0), expr1.dx(i), if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type(c*expr1.dx(i)/expr1.val()*pow(expr1.val(),c)) ));
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow; using Sacado::if_then_else;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return if_then_else( c == scalar_type(1.0), expr1.fastAccessDx(i), if_then_else( expr1.val() == scalar_type(0.0), value_type(0.0), value_type(c*expr1.fastAccessDx(i)/expr1.val()*pow(expr1.val(),c)) ));
      }

    protected:

      const T1& expr1;
      const ConstT& c;
    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, true, false, ExprSpecDefault,
                   PowerImpl::Simd > :
      public Expr< PowerOp< T1, T2, true, false, ExprSpecDefault,
                   PowerImpl::Simd > > {
    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const ConstT& c_, const T2& expr2_) :
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
      value_type val() const {
        using std::pow;
        return pow(c, expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        return if_then_else( c == scalar_type(0.0), value_type(0.0), value_type(expr2.dx(i)*log(c)*pow(c,expr2.val())) );
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        return if_then_else( c == scalar_type(0.0), value_type(0.0), value_type(expr2.fastAccessDx(i)*log(c)*pow(c,expr2.val())) );
      }

    protected:

      const ConstT& c;
      const T2& expr2;
    };

    //
    // Specialization for scalar types using ternary operator
    //

    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, false, ExprSpecDefault,
                   PowerImpl::Scalar > :
      public Expr< PowerOp< T1, T2, false, false, ExprSpecDefault,
                            PowerImpl::Scalar > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const T2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        const int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
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
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type((expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val()));
      }

    protected:

      const T1& expr1;
      const T2& expr2;

    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, true, ExprSpecDefault,
                   PowerImpl::Scalar > :
      public Expr< PowerOp< T1, T2, false, true, ExprSpecDefault,
                            PowerImpl::Scalar > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const ConstT& c_) :
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
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), c);
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return c == scalar_type(1.0) ? expr1.dx(i) : expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type(c*expr1.dx(i)/expr1.val()*pow(expr1.val(),c));
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return c == scalar_type(1.0) ? expr1.fastAccessDx(i) : expr1.val() == scalar_type(0.0) ? value_type(0.0) : value_type(c*expr1.fastAccessDx(i)/expr1.val()*pow(expr1.val(),c));
      }

    protected:

      const T1& expr1;
      const ConstT& c;
    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, true, false, ExprSpecDefault,
                   PowerImpl::Scalar > :
      public Expr< PowerOp< T1, T2, true, false, ExprSpecDefault,
                            PowerImpl::Scalar > > {
    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const ConstT& c_, const T2& expr2_) :
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
      value_type val() const {
        using std::pow;
        return pow(c, expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow; using std::log;
        return c == scalar_type(0.0) ? value_type(0.0) : value_type(expr2.dx(i)*log(c)*pow(c,expr2.val()));
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return c == scalar_type(0.0) ? value_type(0.0) : value_type(expr2.fastAccessDx(i)*log(c)*pow(c,expr2.val()));
      }

    protected:

      const ConstT& c;
      const T2& expr2;
    };

    //
    // Specialization for nested derivatives.  This version does not use
    // if_then_else/ternary-operator on the base so that nested derivatives
    // are correct.
    //
    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, false, ExprSpecDefault,
                   PowerImpl::Nested > :
      public Expr< PowerOp< T1, T2, false, false, ExprSpecDefault,
                   PowerImpl::Nested > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const T2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        const int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
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
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return (expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val());
      }

    protected:

      const T1& expr1;
      const T2& expr2;

    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, true, ExprSpecDefault,
                   PowerImpl::Nested > :
      public Expr< PowerOp< T1, T2, false, true, ExprSpecDefault,
                   PowerImpl::Nested > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const ConstT& c_) :
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
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), c);
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow;
        return c == scalar_type(0.0) ? value_type(0.0) : value_type(c*expr1.dx(i)*pow(expr1.val(),c-scalar_type(1.0)));
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow;
        return c == scalar_type(0.0) ? value_type(0.0) : value_type(c*expr1.fastAccessDx(i)*pow(expr1.val(),c-scalar_type(1.0)));
      }

    protected:

      const T1& expr1;
      const ConstT& c;
    };

    template <typename T1, typename T2>
    class PowerOp<T1, T2, true, false, ExprSpecDefault,
                   PowerImpl::Nested > :
      public Expr< PowerOp< T1, T2, true, false, ExprSpecDefault,
                   PowerImpl::Nested > > {
    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const ConstT& c_, const T2& expr2_) :
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
      value_type val() const {
        using std::pow;
        return pow(c, expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow; using std::log;
        return expr2.dx(i)*log(c)*pow(c,expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return expr2.fastAccessDx(i)*log(c)*pow(c,expr2.val());
      }

    protected:

      const ConstT& c;
      const T2& expr2;
    };

    //
    // Specialization for nested derivatives.  This version does not use
    // if_then_else/ternary-operator on the base so that nested derivatives
    // are correct.
    //
    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, false, ExprSpecDefault,
                   PowerImpl::NestedSimd > :
      public Expr< PowerOp< T1, T2, false, false, ExprSpecDefault,
                   PowerImpl::NestedSimd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const T2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      SACADO_INLINE_FUNCTION
      int size() const {
        const int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      SACADO_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
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
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log; using Sacado::if_then_else;
        return (expr2.fastAccessDx(i)*log(expr1.val())+expr2.val()*expr1.fastAccessDx(i)/expr1.val())*pow(expr1.val(),expr2.val());
      }

    protected:

      const T1& expr1;
      const T2& expr2;

    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, true, ExprSpecDefault,
                   PowerImpl::NestedSimd > :
      public Expr< PowerOp< T1, T2, false, true, ExprSpecDefault,
                   PowerImpl::NestedSimd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const ConstT& c_) :
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
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), c);
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow; using Sacado::if_then_else;
        return if_then_else( c == scalar_type(0.0), value_type(0.0), value_type(c*expr1.dx(i)*pow(expr1.val(),c-scalar_type(1.0))));
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow; using Sacado::if_then_else;
        return if_then_else( c == scalar_type(0.0), value_type(0.0), value_type(c*expr1.fastAccessDx(i)*pow(expr1.val(),c-scalar_type(1.0))));
      }

    protected:

      const T1& expr1;
      const ConstT& c;
    };

    template <typename T1, typename T2>
    class PowerOp<T1, T2, true, false, ExprSpecDefault,
                   PowerImpl::NestedSimd > :
      public Expr< PowerOp< T1, T2, true, false, ExprSpecDefault,
                   PowerImpl::NestedSimd > > {
    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      PowerOp(const ConstT& c_, const T2& expr2_) :
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
      value_type val() const {
        using std::pow;
        return pow(c, expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using std::pow; using std::log;
        return expr2.dx(i)*log(c)*pow(c,expr2.val());
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using std::pow; using std::log;
        return expr2.fastAccessDx(i)*log(c)*pow(c,expr2.val());
      }

    protected:

      const ConstT& c;
      const T2& expr2;
    };

    template <typename T1, typename T2>
    SACADO_INLINE_FUNCTION
    SACADO_FAD_EXP_OP_ENABLE_EXPR_EXPR(PowerOp)
    pow (const T1& expr1, const T2& expr2)
    {
      typedef PowerOp< typename Expr<T1>::derived_type,
                  typename Expr<T2>::derived_type,
                  false, false, typename T1::expr_spec_type > expr_t;

      return expr_t(expr1.derived(), expr2.derived());
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    PowerOp< typename T::value_type, typename Expr<T>::derived_type,
        true, false, typename T::expr_spec_type >
    pow (const typename T::value_type& c,
            const Expr<T>& expr)
    {
      typedef typename T::value_type ConstT;
      typedef PowerOp< ConstT, typename Expr<T>::derived_type,
                  true, false, typename T::expr_spec_type > expr_t;

      return expr_t(c, expr.derived());
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    PowerOp< typename Expr<T>::derived_type, typename T::value_type,
        false, true, typename T::expr_spec_type >
    pow (const Expr<T>& expr,
         const typename T::value_type& c)
    {
      typedef typename T::value_type ConstT;
      typedef PowerOp< typename Expr<T>::derived_type, ConstT,
              false, true, typename T::expr_spec_type > expr_t;

      return expr_t(expr.derived(), c);
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    SACADO_FAD_EXP_OP_ENABLE_SCALAR_EXPR(PowerOp)
    pow (const typename T::scalar_type& c,
         const Expr<T>& expr)
    {
      typedef typename T::scalar_type ConstT;
      typedef PowerOp< ConstT, typename Expr<T>::derived_type,
                       true, false, typename T::expr_spec_type > expr_t;

      return expr_t(c, expr.derived());
    }

    template <typename T>
    SACADO_INLINE_FUNCTION
    SACADO_FAD_EXP_OP_ENABLE_EXPR_SCALAR(PowerOp)
    pow (const Expr<T>& expr,
         const typename T::scalar_type& c)
    {
      typedef typename T::scalar_type ConstT;
      typedef PowerOp< typename Expr<T>::derived_type, ConstT,
                       false, true, typename T::expr_spec_type > expr_t;

      return expr_t(expr.derived(), c);
    }

    template <typename T1, typename T2, bool c1, bool c2, typename E>
    struct ExprLevel< PowerOp< T1, T2, c1, c2, E > > {
      static constexpr unsigned value_1 = ExprLevel<T1>::value;
      static constexpr unsigned value_2 = ExprLevel<T2>::value;
      static constexpr unsigned value =
        value_1 >= value_2 ? value_1 : value_2;
    };

    template <typename T1, typename T2, bool c1, bool c2, typename E>
    struct IsFadExpr< PowerOp< T1, T2, c1, c2, E > > {
      static constexpr unsigned value = true;
    };

  }
  }

  template <typename T1, typename T2, bool c1, bool c2, typename E>
  struct IsExpr< Fad::Exp::PowerOp< T1, T2, c1, c2, E > > {
    static constexpr bool value = true;
  };

  template <typename T1, typename T2, bool c1, bool c2, typename E>
  struct BaseExprType< Fad::Exp::PowerOp< T1, T2, c1, c2, E > > {
    typedef typename BaseExprType<T1>::type base_expr_1;
    typedef typename BaseExprType<T2>::type base_expr_2;
    typedef typename Sacado::Promote<base_expr_1,
                                     base_expr_2>::type type;
  };

  template <typename T1, typename T2, bool c1, bool c2, typename E>
  struct IsSimdType< Fad::Exp::PowerOp< T1, T2, c1, c2, E > > {
    static const bool value =
      IsSimdType< typename Fad::Exp::PowerOp< T1, T2, c1, c2, E >::value_type >::value;
  };

  template <typename T1, typename T2, bool c1, bool c2, typename E>
  struct ValueType< Fad::Exp::PowerOp< T1, T2, c1, c2, E > > {
    typedef typename Fad::Exp::PowerOp< T1, T2, c1, c2, E >::value_type type;
  };

  template <typename T1, typename T2, bool c1, bool c2, typename E>
  struct ScalarType< Fad::Exp::PowerOp< T1, T2, c1, c2, E > > {
    typedef typename Fad::Exp::PowerOp< T1, T2, c1, c2, E >::scalar_type type;
  };

}

//--------------------------if_then_else operator -----------------------
// Can't use the above macros because it is a ternary operator (sort of).
// Also, relies on C++11

namespace Sacado {
  namespace Fad {
  namespace Exp {

    template <typename CondT, typename T1, typename T2,
              bool is_const_T1, bool is_const_T2,
              typename ExprSpec>
    class IfThenElseOp {};

    template <typename CondT, typename T1, typename T2>
    class IfThenElseOp< CondT, T1, T2, false, false, ExprSpecDefault > :
      public Expr< IfThenElseOp< CondT, T1, T2, false, false, ExprSpecDefault > > {

    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      IfThenElseOp(const CondT& cond_, const T1& expr1_, const T2& expr2_) :
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
      value_type val() const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.val(), expr2.val() );
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.dx(i), expr2.dx(i) );
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.fastAccessDx(i), expr2.fastAccessDx(i) );
      }

    protected:

      const CondT&  cond;
      const T1& expr1;
      const T2& expr2;

    };

    template <typename CondT, typename T1, typename T2>
    class IfThenElseOp< CondT, T1, T2, false, true, ExprSpecDefault > :
      public Expr< IfThenElseOp< CondT, T1, T2, false, true, ExprSpecDefault > > {

    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      IfThenElseOp(const CondT& cond_, const T1& expr1_, const ConstT& c_) :
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
      value_type val() const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.val(), c );
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.dx(i), value_type(0.0) );
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, expr1.fastAccessDx(i), value_type(0.0) );
      }

    protected:

      const CondT&  cond;
      const T1& expr1;
      const ConstT& c;
    };

    template <typename CondT, typename T1, typename T2>
    class IfThenElseOp< CondT, T1, T2, true, false, ExprSpecDefault > :
      public Expr< IfThenElseOp< CondT, T1, T2, true, false, ExprSpecDefault > > {

    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef ExprSpecDefault expr_spec_type;

      SACADO_INLINE_FUNCTION
      IfThenElseOp(const CondT& cond_, const ConstT& c_, const T2& expr2_) :
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
      value_type val() const {
        using Sacado::if_then_else;
        return if_then_else( cond, c, expr2.val() );
      }

      SACADO_INLINE_FUNCTION
      value_type dx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, value_type(0.0), expr2.dx(i) );
      }

      SACADO_INLINE_FUNCTION
      value_type fastAccessDx(int i) const {
        using Sacado::if_then_else;
        return if_then_else( cond, value_type(0.0), expr2.fastAccessDx(i) );
      }

    protected:

      const CondT&  cond;
      const ConstT& c;
      const T2& expr2;
    };

    template <typename CondT, typename T1, typename T2>
    SACADO_INLINE_FUNCTION
    typename mpl::enable_if_c< IsFadExpr<T1>::value && IsFadExpr<T2>::value &&
                               ExprLevel<T1>::value == ExprLevel<T2>::value,
                               IfThenElseOp< CondT,
                                             typename Expr<T1>::derived_type,
                                             typename Expr<T2>::derived_type,
                                             false, false,
                                             typename T1::expr_spec_type >
                             >::type
    if_then_else (const CondT& cond, const T1& expr1, const T2& expr2)
    {
      typedef IfThenElseOp< CondT, typename Expr<T1>::derived_type,
                            typename Expr<T2>::derived_type,
                            false, false, typename T1::expr_spec_type > expr_t;

      return expr_t(cond, expr1.derived(), expr2.derived());
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    IfThenElseOp< CondT, typename T::value_type, typename Expr<T>::derived_type,
                  true, false, typename T::expr_spec_type >
    if_then_else (const CondT& cond, const typename T::value_type& c,
                  const Expr<T>& expr)
    {
      typedef typename T::value_type ConstT;
      typedef IfThenElseOp< CondT, ConstT, typename Expr<T>::derived_type,
                            true, false, typename T::expr_spec_type > expr_t;

      return expr_t(cond, c, expr.derived());
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    IfThenElseOp< CondT, typename Expr<T>::derived_type, typename T::value_type,
                  false, true, typename T::expr_spec_type >
    if_then_else (const CondT& cond, const Expr<T>& expr,
                  const typename T::value_type& c)
    {
      typedef typename T::value_type ConstT;
      typedef IfThenElseOp< CondT, typename Expr<T>::derived_type, ConstT,
                            false, true, typename T::expr_spec_type > expr_t;

      return expr_t(cond, expr.derived(), c);
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    typename mpl::disable_if<
      std::is_same< typename T::value_type,
                    typename T::scalar_type >,
      IfThenElseOp< CondT, typename T::scalar_type,
                    typename Expr<T>::derived_type,
                    true, false, typename T::expr_spec_type >
      >::type
    if_then_else (const CondT& cond, const typename Expr<T>::scalar_type& c,
                  const Expr<T>& expr)
    {
      typedef typename T::scalar_type ConstT;
      typedef IfThenElseOp< CondT, ConstT, typename Expr<T>::derived_type,
                            true, false, typename T::expr_spec_type > expr_t;

      return expr_t(cond, c, expr.derived());
    }

    template <typename CondT, typename T>
    SACADO_INLINE_FUNCTION
    typename mpl::disable_if<
      std::is_same< typename T::value_type,
                    typename T::scalar_type >,
      IfThenElseOp< CondT, typename Expr<T>::derived_type,
                    typename T::scalar_type,
                    false, true, typename T::expr_spec_type >
      >::type
    if_then_else (const CondT& cond, const Expr<T>& expr,
                  const typename Expr<T>::scalar_type& c)
    {
      typedef typename T::scalar_type ConstT;
      typedef IfThenElseOp< CondT, typename Expr<T>::derived_type, ConstT,
                            false, true, typename T::expr_spec_type > expr_t;

      return expr_t(cond, expr.derived(), c);
    }

    template <typename CondT, typename T1, typename T2, bool c1, bool c2,
              typename E>
    struct ExprLevel< IfThenElseOp< CondT, T1, T2, c1, c2, E > > {
      static constexpr unsigned value_1 = ExprLevel<T1>::value;
      static constexpr unsigned value_2 = ExprLevel<T2>::value;
      static constexpr unsigned value =
        value_1 >= value_2 ? value_1 : value_2;
    };

    template <typename CondT, typename T1, typename T2, bool c1, bool c2,
              typename E>
    struct IsFadExpr< IfThenElseOp< CondT, T1, T2, c1, c2, E > > {
      static constexpr unsigned value = true;
    };

  }
  }

  template <typename CondT, typename T1, typename T2, bool c1, bool c2,
            typename E>
  struct IsExpr< Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E > > {
    static constexpr bool value = true;
  };

  template <typename CondT, typename T1, typename T2, bool c1, bool c2,
            typename E>
  struct BaseExprType< Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E > > {
    typedef typename BaseExprType<T1>::type base_expr_1;
    typedef typename BaseExprType<T2>::type base_expr_2;
    typedef typename Sacado::Promote<base_expr_1,
                                     base_expr_2>::type type;
  };

  template <typename CondT, typename T1, typename T2, bool c1, bool c2,
            typename E>
  struct IsSimdType< Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E > > {
    static const bool value =
      IsSimdType< typename Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E >::value_type >::value;
  };

  template <typename CondT, typename T1, typename T2, bool c1, bool c2,
            typename E>
  struct ValueType< Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E > > {
    typedef typename Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E >::value_type type;
  };

  template <typename CondT, typename T1, typename T2, bool c1, bool c2,
            typename E>
  struct ScalarType< Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E > > {
    typedef typename Fad::Exp::IfThenElseOp< CondT, T1, T2, c1, c2, E >::scalar_type type;
  };
}

#undef FAD_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

namespace Sacado {
  namespace Fad {
  namespace Exp {
  namespace Impl {
    // Helper trait to determine return type of logical comparison operations
    // (==, !=, ...), usually bool but maybe something else for SIMD types.
    // Need to use SFINAE to restrict to types that define == operator in the
    // conditional overloads below, otherwise instantiating ConditionaReturnType
    // may fail during overload resolution.
    template <typename T1, typename T2 = T1,
              bool = mpl::has_equal_to<T1,T2>::value>
    struct ConditionalReturnType {};

    template <typename T1, typename T2>
    struct ConditionalReturnType<T1,T2,true> {
      typedef decltype( std::declval<T1>() == std::declval<T2>() ) type;
    };
  }
  }
  }
}

#define FAD_RELOP_MACRO(OP)                                             \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    typename mpl::enable_if_c<                                          \
       IsFadExpr<T1>::value && IsFadExpr<T2>::value &&                  \
       ExprLevel<T1>::value == ExprLevel<T2>::value,                    \
       typename Impl::ConditionalReturnType<typename T1::value_type,    \
                                            typename T2::value_type>::type \
       >::type                                                          \
    operator OP (const T1& expr1, const T2& expr2)                      \
    {                                                                   \
      return expr1.derived().val() OP expr2.derived().val();            \
    }                                                                   \
                                                                        \
    template <typename T2>                                              \
    SACADO_INLINE_FUNCTION                                              \
    typename Impl::ConditionalReturnType<typename T2::value_type>::type \
    operator OP (const typename T2::value_type& a,                      \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return a OP expr2.derived().val();                                \
    }                                                                   \
                                                                        \
    template <typename T1>                                              \
    SACADO_INLINE_FUNCTION                                              \
    typename Impl::ConditionalReturnType<typename T1::value_type>::type \
    operator OP (const Expr<T1>& expr1,                                 \
                 const typename T1::value_type& b)                      \
    {                                                                   \
      return expr1.derived().val() OP b;                                \
    }                                                                   \
  }                                                                     \
  }                                                                     \
}                                                                       \

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
  namespace Exp {

    template <typename ExprT>
    SACADO_INLINE_FUNCTION
    bool operator ! (const Expr<ExprT>& expr)
    {
      return ! expr.derived().val();
    }

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace Fad {
  namespace Exp {

    template <typename T>
    SACADO_INLINE_FUNCTION
    bool toBool(const Expr<T>& xx) {
      const typename Expr<T>::derived_type& x = xx.derived();
      bool is_zero = (x.val() == 0.0);
      for (int i=0; i<x.size(); i++)
        is_zero = is_zero && (x.dx(i) == 0.0);
      return !is_zero;
    }

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#define FAD_BOOL_MACRO(OP)                                              \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<T1>& expr1,                                 \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return toBool(expr1) OP toBool(expr2);                            \
    }                                                                   \
                                                                        \
    template <typename T2>                                              \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename Expr<T2>::value_type& a,                \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return a OP toBool(expr2);                                        \
    }                                                                   \
                                                                        \
    template <typename T1>                                              \
    SACADO_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<T1>& expr1,                                 \
                 const typename Expr<T1>::value_type& b)                \
    {                                                                   \
      return toBool(expr1) OP b;                                        \
    }                                                                   \
  }                                                                     \
  }                                                                     \
}

FAD_BOOL_MACRO(&&)
FAD_BOOL_MACRO(||)

#undef FAD_BOOL_MACRO

//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace Fad {
  namespace Exp {

    template <typename T>
    std::ostream& operator << (std::ostream& os, const Expr<T>& xx) {
      const typename Expr<T>::derived_type& x = xx.derived();
      os << x.val() << " [";

      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#if defined(HAVE_SACADO_KOKKOS)

//-------------------------- Atomic Operators -----------------------

namespace Sacado {

  namespace Fad {
  namespace Exp {

    // Overload of Kokkos::atomic_add for Fad types.
    template <typename S>
    SACADO_INLINE_FUNCTION
    void atomic_add(GeneralFad<S>* dst, const GeneralFad<S>& x) {
      using Kokkos::atomic_add;

      const int xsz = x.size();
      const int sz = dst->size();

      // We currently cannot handle resizing since that would need to be
      // done atomically.
      if (xsz > sz)
        Kokkos::abort(
          "Sacado error: Fad resize within atomic_add() not supported!");

      if (xsz != sz && sz > 0 && xsz > 0)
        Kokkos::abort(
          "Sacado error: Fad assignment of incompatiable sizes!");


      if (sz > 0 && xsz > 0) {
        SACADO_FAD_DERIV_LOOP(i,sz)
          atomic_add(&(dst->fastAccessDx(i)), x.fastAccessDx(i));
      }
      SACADO_FAD_THREAD_SINGLE
        atomic_add(&(dst->val()), x.val());
    }

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // HAVE_SACADO_KOKKOS

#endif // SACADO_FAD_OPS_HPP
