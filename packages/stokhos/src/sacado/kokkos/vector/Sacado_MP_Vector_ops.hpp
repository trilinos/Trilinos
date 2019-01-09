// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Sacado_cmath.hpp"
#include <ostream>      // for std::ostream

#ifdef __CUDACC__
#include <math_functions.h>
#endif

/*
namespace Sacado {
  namespace MP {

    template <typename T>
    class LogOp :
      public Expr< LogOp< T > > {
    public:

      typedef typename T::value_type value_type;
      typedef typename T::storage_type storage_type;

      KOKKOS_INLINE_FUNCTION
      LogOp(const T& expr_) : expr(expr_)  {}

      KOKKOS_INLINE_FUNCTION
      std::string name() const {
        return std::string("log") + expr.name();
      }

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess(int sz) const {
        return expr.hasFastAccess(sz);
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        return std::log(expr.val());
      }

      KOKKOS_INLINE_FUNCTION
      value_type coeff(int i) const {
        return std::log(expr.coeff(i));
      }

      KOKKOS_INLINE_FUNCTION
      value_type fastAccessCoeff(int i) const {
        return std::log(expr.fastAccessCoeff(i));
      }

    protected:

      const T& expr;

    };

    template <typename T, typename N>
    KOKKOS_INLINE_FUNCTION
    LogOp< T >
    log (const Expr<T>& expr)
    {
      typedef LogOp< typename Expr<T>::derived_type > expr_t;

      return expr_t(expr.derived());
    }
  }
}
*/

#define MP_UNARYOP_MACRO(OPNAME,OP,OPER)                                \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
                                                                        \
    template <typename T>                                               \
    class OP :                                                          \
      public Expr< OP< T > > {                                          \
    public:                                                             \
                                                                        \
      typedef typename remove_volatile<T>::type Tnv;                    \
      typedef typename Tnv::value_type value_type;                      \
      typedef typename Tnv::storage_type storage_type;                  \
      typedef typename Tnv::base_expr_type base_expr_type;              \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T& expr_) : expr(expr_)  {}                              \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      std::string name() const {                                        \
        return std::string(#OPER) + expr.name();                        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr.size();                                             \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess(int sz) const {                                \
        return expr.hasFastAccess(sz);                                  \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return OPER(expr.val());                                        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type coeff(int i) const {                                   \
        return OPER(expr.coeff(i));                                     \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type fastAccessCoeff(int i) const {                         \
        return OPER(expr.fastAccessCoeff(i));                           \
      }                                                                 \
                                                                        \
      template <int i>                                                  \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type getCoeff() const {                                     \
        return OPER(expr.template getCoeff<i>());                       \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      typename const_expr_ref<T>::type expr;                            \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< T >                                                             \
    OPNAME (const Expr<T>& expr)                                        \
    {                                                                   \
      typedef OP< typename Expr<T>::derived_type > expr_t;              \
                                                                        \
      return expr_t(expr.derived());                                    \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< volatile T >                                                    \
    OPNAME (const volatile Expr<T>& expr)                               \
    {                                                                   \
      typedef typename Expr<T>::derived_type derived;                   \
      typedef OP< typename add_volatile<derived>::type > expr_t;        \
                                                                        \
      return expr_t(expr.derived());                                    \
    }                                                                   \
  }                                                                     \
                                                                        \
  template <typename T>                                                 \
  struct IsExpr< MP::OP<T> > {                                          \
    static const bool value = true;                                     \
  };                                                                    \
                                                                        \
  template <typename T>                                                 \
  struct BaseExprType< MP::OP<T> > {                                    \
    typedef typename MP::OP<T>::base_expr_type type;                    \
  };                                                                    \
}

MP_UNARYOP_MACRO(operator+, UnaryPlusOp, +)
MP_UNARYOP_MACRO(operator-, UnaryMinusOp, -)
MP_UNARYOP_MACRO(exp, ExpOp, std::exp)
MP_UNARYOP_MACRO(log, LogOp, std::log)
MP_UNARYOP_MACRO(log10, Log10Op, std::log10)
MP_UNARYOP_MACRO(sqrt, SqrtOp, std::sqrt)
MP_UNARYOP_MACRO(cbrt, CbrtOp, std::cbrt)
MP_UNARYOP_MACRO(cos, CosOp, std::cos)
MP_UNARYOP_MACRO(sin, SinOp, std::sin)
MP_UNARYOP_MACRO(tan, TanOp, std::tan)
MP_UNARYOP_MACRO(acos, ACosOp, std::acos)
MP_UNARYOP_MACRO(asin, ASinOp, std::asin)
MP_UNARYOP_MACRO(atan, ATanOp, std::atan)
MP_UNARYOP_MACRO(cosh, CoshOp, std::cosh)
MP_UNARYOP_MACRO(sinh, SinhOp, std::sinh)
MP_UNARYOP_MACRO(tanh, TanhOp, std::tanh)
MP_UNARYOP_MACRO(acosh, ACoshOp, std::acosh)
MP_UNARYOP_MACRO(asinh, ASinhOp, std::asinh)
MP_UNARYOP_MACRO(atanh, ATanhOp, std::atanh)
MP_UNARYOP_MACRO(abs, AbsOp, std::abs)
MP_UNARYOP_MACRO(fabs, FAbsOp, std::fabs)

#undef MP_UNARYOP_MACRO

#define MP_BINARYOP_MACRO(OPNAME,OP,OPER)                               \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename T1, typename T2>                                 \
    class OP :                                                          \
      public Expr< OP< T1, T2> > {                                      \
                                                                        \
    public:                                                             \
                                                                        \
      typedef typename remove_volatile<T1>::type Tnv1;                  \
      typedef typename remove_volatile<T2>::type Tnv2;                  \
      typedef typename Tnv1::value_type value_type_1;                   \
      typedef typename Tnv2::value_type value_type_2;                   \
      typedef typename Sacado::Promote<value_type_1,                    \
                                       value_type_2>::type value_type;  \
                                                                        \
      typedef typename Tnv1::storage_type storage_type;                 \
      typedef typename Tnv1::base_expr_type base_expr_type;             \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const T2& expr2_) :                          \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      std::string name() const {                                        \
        return expr1.name() + std::string(#OPER) + expr2.name();        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        int sz1 = expr1.size(), sz2 = expr2.size();                     \
        return sz1 > sz2 ? sz1 : sz2;                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess(int sz) const {                                \
        return expr1.hasFastAccess(sz) && expr2.hasFastAccess(sz);      \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return (expr1.val() OPER expr2.val());                          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type coeff(int i) const {                                   \
        return (expr1.coeff(i) OPER expr2.coeff(i));                    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type fastAccessCoeff(int i) const {                         \
        return (expr1.fastAccessCoeff(i) OPER expr2.fastAccessCoeff(i)); \
      }                                                                 \
                                                                        \
      template <int i>                                                  \
        KOKKOS_INLINE_FUNCTION                                          \
      value_type getCoeff() const {                                     \
        return expr1.template getCoeff<i>() OPER expr2.template getCoeff<i>(); \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      typename const_expr_ref<T1>::type expr1;                          \
      typename const_expr_ref<T2>::type expr2;                          \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename T1>                                              \
    class OP< T1, typename T1::value_type > :                           \
      public Expr< OP< T1, typename T1::value_type > > {                \
                                                                        \
    public:                                                             \
                                                                        \
      typedef typename remove_volatile<T1>::type Tnv1;                  \
      typedef typename Tnv1::value_type value_type;                     \
      typedef typename Tnv1::value_type ConstT;                         \
                                                                        \
      typedef typename Tnv1::storage_type storage_type;                 \
      typedef typename Tnv1::base_expr_type base_expr_type;             \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const ConstT& c_) :                          \
        expr1(expr1_), c(c_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      std::string name() const {                                        \
        return expr1.name() + std::string(#OPER) + std::string("c");    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr1.size();                                            \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess(int sz) const {                                \
        return expr1.hasFastAccess(sz);                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return (expr1.val() OPER c);                                    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type coeff(int i) const {                                   \
        return (expr1.coeff(i) OPER c);                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type fastAccessCoeff(int i) const {                         \
        return (expr1.fastAccessCoeff(i) OPER c);                       \
      }                                                                 \
                                                                        \
      template <int i>                                                  \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type getCoeff() const {                                     \
        return expr1.template getCoeff<i>() OPER c;                     \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      typename const_expr_ref<T1>::type expr1;                          \
      const ConstT& c;                                                  \
    };                                                                  \
                                                                        \
    template <typename T2>                                              \
    class OP< typename T2::value_type, T2 > :                           \
      public Expr< OP< typename T2::value_type, T2 > > {                \
                                                                        \
    public:                                                             \
                                                                        \
      typedef typename remove_volatile<T2>::type Tnv2;                  \
      typedef typename Tnv2::value_type value_type;                     \
      typedef typename Tnv2::value_type ConstT;                         \
                                                                        \
      typedef typename Tnv2::storage_type storage_type;                 \
      typedef typename Tnv2::base_expr_type base_expr_type;             \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const ConstT& c_, const T2& expr2_) :                          \
        c(c_), expr2(expr2_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      std::string name() const {                                        \
        return std::string("c") + std::string(#OPER) + expr2.name();    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const { return expr2.size(); }                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess(int sz) const {                                \
        return expr2.hasFastAccess(sz);                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return (c OPER expr2.val());                                    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type coeff(int i) const {                                   \
        return (c OPER expr2.coeff(i));                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type fastAccessCoeff(int i) const {                         \
        return (c OPER expr2.fastAccessCoeff(i));                       \
      }                                                                 \
                                                                        \
      template <int i>                                                  \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type getCoeff() const {                                     \
        return c OPER expr2.template getCoeff<i>();                     \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ConstT& c;                                                  \
      typename const_expr_ref<T2>::type expr2;                          \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< T1, T2 >                                                        \
    OPNAME (const Expr<T1>& expr1,                                      \
            const Expr<T2>& expr2)                                      \
    {                                                                   \
      typedef OP< typename Expr<T1>::derived_type,                      \
                  typename Expr<T2>::derived_type > expr_t;             \
                                                                        \
      return expr_t(expr1.derived(), expr2.derived());                  \
    }                                                                   \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< volatile T1, volatile T2 >                                      \
    OPNAME (const volatile Expr<T1>& expr1,                             \
            const volatile Expr<T2>& expr2)                             \
    {                                                                   \
      typedef typename Expr<T1>::derived_type derived1;                 \
      typedef typename Expr<T2>::derived_type derived2;                 \
      typedef OP< typename add_volatile<derived1>::type,                \
                  typename add_volatile<derived2>::type > expr_t;       \
                                                                        \
      return expr_t(expr1.derived(), expr2.derived());                  \
    }                                                                   \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< T1, volatile T2 >                                               \
    OPNAME (const Expr<T1>& expr1,                                      \
            const volatile Expr<T2>& expr2)                             \
    {                                                                   \
      typedef typename Expr<T1>::derived_type derived1;                 \
      typedef typename Expr<T2>::derived_type derived2;                 \
      typedef OP< derived1,                                             \
                  typename add_volatile<derived2>::type > expr_t;       \
                                                                        \
      return expr_t(expr1.derived(), expr2.derived());                  \
    }                                                                   \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< volatile T1, T2 >                                               \
    OPNAME (const volatile Expr<T1>& expr1,                             \
            const Expr<T2>& expr2)                                      \
    {                                                                   \
      typedef typename Expr<T1>::derived_type derived1;                 \
      typedef typename Expr<T2>::derived_type derived2;                 \
      typedef OP< typename add_volatile<derived1>::type,                \
                  derived2 > expr_t;                                    \
                                                                        \
      return expr_t(expr1.derived(), expr2.derived());                  \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< typename T::value_type, T >                                     \
    OPNAME (const typename T::value_type& c,                            \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef OP< ConstT, typename Expr<T>::derived_type > expr_t;      \
                                                                        \
      return expr_t(c, expr.derived());                                 \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< typename T::value_type, volatile T >                            \
    OPNAME (const typename T::value_type& c,                            \
            const volatile Expr<T>& expr)                               \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef typename Expr<T>::derived_type derived;                   \
      typedef OP< ConstT,                                               \
                  typename add_volatile<derived>::type > expr_t;        \
                                                                        \
      return expr_t(c, expr.derived());                                 \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< T, typename T::value_type >                                     \
    OPNAME (const Expr<T>& expr,                                        \
            const typename T::value_type& c)                            \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef OP< typename Expr<T>::derived_type, ConstT > expr_t;      \
                                                                        \
      return expr_t(expr.derived(), c);                                 \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< volatile T, typename T::value_type >                            \
    OPNAME (const volatile Expr<T>& expr,                               \
            const typename T::value_type& c)                            \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef typename Expr<T>::derived_type derived;                   \
      typedef OP< typename add_volatile<derived>::type,                 \
                  ConstT > expr_t;                                      \
                                                                        \
      return expr_t(expr.derived(), c);                                 \
    }                                                                   \
  }                                                                     \
                                                                        \
  template <typename T1, typename T2>                                   \
  struct IsExpr< MP::OP<T1,T2> > {                                      \
    static const bool value = true;                                     \
  };                                                                    \
                                                                        \
  template <typename T1, typename T2>                                   \
  struct BaseExprType< MP::OP<T1,T2> > {                                \
    typedef typename MP::OP<T1,T2>::base_expr_type type;                \
  };                                                                    \
}

MP_BINARYOP_MACRO(operator+, AdditionOp, +)
MP_BINARYOP_MACRO(operator-, SubtractionOp, -)
MP_BINARYOP_MACRO(operator*, MultiplicationOp, *)
MP_BINARYOP_MACRO(operator/, DivisionOp, /)

#undef MP_BINARYOP_MACRO

#define MP_BINARYOP_MACRO(OPNAME,OP,OPER)                               \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename T1, typename T2>                                 \
    class OP :                                                          \
      public Expr< OP< T1, T2 > > {                                     \
                                                                        \
    public:                                                             \
                                                                        \
      typedef typename T1::value_type value_type_1;                     \
      typedef typename T2::value_type value_type_2;                     \
      typedef typename Sacado::Promote<value_type_1,                    \
                                       value_type_2>::type value_type;  \
                                                                        \
      typedef typename T1::storage_type storage_type;                   \
      typedef typename T1::base_expr_type base_expr_type;               \
                                                                        \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const T2& expr2_) :                          \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      std::string name() const {                                        \
        return expr1.name() + std::string(#OPER) + expr2.name();        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        int sz1 = expr1.size(), sz2 = expr2.size();                     \
        return sz1 > sz2 ? sz1 : sz2;                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess(int sz) const {                                \
        return expr1.hasFastAccess(sz) && expr2.hasFastAccess(sz);      \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return OPER(expr1.val(), expr2.val());                          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type coeff(int i) const {                                   \
        return OPER(expr1.coeff(i), expr2.coeff(i));                    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type fastAccessCoeff(int i) const {                         \
        return OPER(expr1.fastAccessCoeff(i), expr2.fastAccessCoeff(i)); \
      }                                                                 \
                                                                        \
      template <int i>                                                  \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type getCoeff() const {                                     \
        return OPER(expr1.template getCoeff<i>(),                       \
                    expr2.template getCoeff<i>());                      \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      typename const_expr_ref<T1>::type expr1;                          \
      typename const_expr_ref<T2>::type expr2;                          \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename T1>                                              \
    class OP< T1, typename T1::value_type > :                           \
      public Expr< OP< T1, typename T1::value_type > > {                \
                                                                        \
    public:                                                             \
                                                                        \
      typedef typename T1::value_type value_type;                       \
      typedef typename T1::value_type ConstT;                           \
                                                                        \
      typedef typename T1::storage_type storage_type;                   \
      typedef typename T1::base_expr_type base_expr_type;               \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const ConstT& c_) :                          \
        expr1(expr1_), c(c_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      std::string name() const {                                        \
        return expr1.name() + std::string(#OPER) + std::string("c");    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const { return expr1.size(); }                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess(int sz) const {                                \
        return expr1.hasFastAccess(sz);                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return OPER(expr1.val(), c);                                    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type coeff(int i) const {                                   \
        return OPER(expr1.coeff(i), c);                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type fastAccessCoeff(int i) const {                         \
        return OPER(expr1.fastAccessCoeff(i), c);                       \
      }                                                                 \
                                                                        \
      template <int i>                                                  \
        KOKKOS_INLINE_FUNCTION                                          \
      value_type getCoeff() const {                                     \
        return OPER(expr1.template getCoeff<i>(), c);                   \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      typename const_expr_ref<T1>::type expr1;                          \
      const ConstT& c;                                                  \
    };                                                                  \
                                                                        \
    template <typename T2>                                              \
    class OP< typename T2::value_type, T2 > :                           \
      public Expr< OP< typename T2::value_type, T2 > > {                \
                                                                        \
    public:                                                             \
                                                                        \
      typedef typename T2::value_type value_type;                       \
      typedef typename T2::value_type ConstT;                           \
                                                                        \
      typedef typename T2::storage_type storage_type;                   \
      typedef typename T2::base_expr_type base_expr_type;               \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const ConstT& c_, const T2& expr2_) :                          \
        c(c_), expr2(expr2_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      std::string name() const {                                        \
        return std::string("c") + std::string(#OPER) + expr2.name();    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const { return expr2.size(); }                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess(int sz) const {                                \
        return expr2.hasFastAccess(sz);                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        return OPER(c, expr2.val());                                    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type coeff(int i) const {                                   \
        return OPER(c, expr2.coeff(i));                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type fastAccessCoeff(int i) const {                         \
        return OPER(c, expr2.fastAccessCoeff(i));                       \
      }                                                                 \
                                                                        \
      template <int i>                                                  \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type getCoeff() const {                                     \
        return OPER(c, expr2.template getCoeff<i>());                   \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const ConstT& c;                                                  \
      typename const_expr_ref<T2>::type expr2;                          \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< T1, T2 >                                                        \
    OPNAME (const Expr<T1>& expr1,                                      \
            const Expr<T2>& expr2)                                      \
    {                                                                   \
      typedef OP< typename Expr<T1>::derived_type,                      \
                  typename Expr<T2>::derived_type > expr_t;             \
                                                                        \
      return expr_t(expr1.derived(), expr2.derived());                  \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< typename T::value_type, T >                                     \
    OPNAME (const typename T::value_type& c,                            \
            const Expr<T>& expr)                                        \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef OP< ConstT, typename Expr<T>::derived_type > expr_t;      \
                                                                        \
      return expr_t(c, expr.derived());                                 \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    OP< T, typename T::value_type >                                     \
    OPNAME (const Expr<T>& expr,                                        \
            const typename T::value_type& c)                            \
    {                                                                   \
      typedef typename T::value_type ConstT;                            \
      typedef OP< typename Expr<T>::derived_type, ConstT > expr_t;      \
                                                                        \
      return expr_t(expr.derived(), c);                                 \
    }                                                                   \
  }                                                                     \
                                                                        \
  template <typename T1, typename T2>                                   \
  struct IsExpr< MP::OP<T1,T2> > {                                      \
    static const bool value = true;                                     \
  };                                                                    \
                                                                        \
  template <typename T1, typename T2>                                   \
  struct BaseExprType< MP::OP<T1,T2> > {                                \
    typedef typename MP::OP<T1,T2>::base_expr_type type;                \
  };                                                                    \
}

MP_BINARYOP_MACRO(atan2, Atan2Op, std::atan2)
MP_BINARYOP_MACRO(pow, PowerOp, std::pow)
#ifdef __CUDACC__
MP_BINARYOP_MACRO(max, MaxOp, ::max)
MP_BINARYOP_MACRO(min, MinOp, ::min)
#else
MP_BINARYOP_MACRO(max, MaxOp, std::max)
MP_BINARYOP_MACRO(min, MinOp, std::min)
#endif

#undef MP_BINARYOP_MACRO

//-------------------------- Relational Operators -----------------------

#define MP_RELOP_MACRO(OP)                                              \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<T1>& expr1,                                 \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return expr1.derived().val() OP expr2.derived().val();            \
    }                                                                   \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const volatile Expr<T1>& expr1,                        \
                 const volatile Expr<T2>& expr2)                        \
    {                                                                   \
      return expr1.derived().val() OP expr2.derived().val();            \
    }                                                                   \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const volatile Expr<T1>& expr1,                        \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return expr1.derived().val() OP expr2.derived().val();            \
    }                                                                   \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<T1>& expr1,                                 \
                 const volatile Expr<T2>& expr2)                        \
    {                                                                   \
      return expr1.derived().val() OP expr2.derived().val();            \
    }                                                                   \
                                                                        \
    template <typename T2>                                              \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename T2::value_type& a,                      \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return a OP expr2.derived().val();                                \
    }                                                                   \
                                                                        \
    template <typename T2>                                              \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename T2::value_type& a,                      \
                 const volatile Expr<T2>& expr2)                        \
    {                                                                   \
      return a OP expr2.derived().val();                                \
    }                                                                   \
                                                                        \
    template <typename T1>                                              \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<T1>& expr1,                                 \
                 const typename T1::value_type& b)                      \
    {                                                                   \
      return expr1.derived().val() OP b;                                \
    }                                                                   \
                                                                        \
    template <typename T1>                                              \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const volatile Expr<T1>& expr1,                        \
                 const typename T1::value_type& b)                      \
    {                                                                   \
      return expr1.derived().val() OP b;                                \
    }                                                                   \
  }                                                                     \
}

MP_RELOP_MACRO(==)
MP_RELOP_MACRO(!=)
MP_RELOP_MACRO(<)
MP_RELOP_MACRO(>)
MP_RELOP_MACRO(<=)
MP_RELOP_MACRO(>=)
MP_RELOP_MACRO(<<=)
MP_RELOP_MACRO(>>=)
MP_RELOP_MACRO(&)
MP_RELOP_MACRO(|)

#undef MP_RELOP_MACRO

namespace Sacado {

  namespace MP {

    template <typename T>
    KOKKOS_INLINE_FUNCTION
    bool operator ! (const Expr<T>& expr)
    {
      return ! expr.derived().val();
    }

  } // namespace MP

} // namespace Sacado


//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace MP {

    template <typename T>
    KOKKOS_INLINE_FUNCTION
    bool toBool(const Expr<T>& xx) {
      const typename Expr<T>::derived_type& x =
        xx.derived();
      bool is_zero = true;
      for (int i=0; i<x.size(); i++)
        is_zero = is_zero && (x.coeff(i) == 0.0);
      return !is_zero;
    }

  } // namespace MP

} // namespace Sacado

#define PCE_BOOL_MACRO(OP)                                              \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<T1>& expr1,                                 \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return toBool(expr1) OP toBool(expr2);                            \
    }                                                                   \
                                                                        \
    template <typename T2>                                              \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename T2::value_type& a,                      \
                 const Expr<T2>& expr2)                                 \
    {                                                                   \
      return a OP toBool(expr2);                                        \
    }                                                                   \
                                                                        \
    template <typename T1>                                              \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Expr<T1>& expr1,                                 \
                 const typename T1::value_type& b)                      \
    {                                                                   \
      return toBool(expr1) OP b;                                        \
    }                                                                   \
  }                                                                     \
}

PCE_BOOL_MACRO(&&)
PCE_BOOL_MACRO(||)

#undef PCE_BOOL_MACRO


//-------------------------- I/O Operators -----------------------

namespace Sacado {

  namespace MP {

    template <typename T>
    std::ostream& operator << (std::ostream& os,
                               const Expr<T>& x) {
      typedef typename T::storage_type storage_type;
      Vector<storage_type> a(x);
      os << a;
      return os;
    }

  } // namespace MP

} // namespace Sacado

//-------------------------- Standard library -----------------------
namespace std {
  template <typename T>
  bool isfinite(const Sacado::MP::Expr<T>& xx) {
    using std::isfinite;
    const typename Sacado::MP::Expr<T>::derived_type& x = xx.derived();
    for (int i=0; i<x.size(); i++)
      if (!isfinite(x.coeff(i)))
        return false;
      return true;
  }
}
