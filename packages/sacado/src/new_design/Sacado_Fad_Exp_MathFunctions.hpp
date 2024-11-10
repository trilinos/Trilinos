// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_MATHFUNCTIONS_HPP
#define SACADO_FAD_EXP_MATHFUNCTIONS_HPP

#include "Sacado_cmath.hpp"
#include "Sacado_SFINAE_Macros.hpp"

// Note:  Sacado::Fad::Ops are forward-declared here, instead of in macros
// below.
#include "Sacado_Fad_Exp_Ops_Fwd.hpp"

#define UNARYFUNC_MACRO(OP,FADOP)                                       \
namespace Sacado {                                                      \
                                                                        \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
    template <typename T> class Expr;                                   \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    FADOP< typename Expr<T>::derived_type,                              \
           typename T::expr_spec_type >                                 \
    OP (const Expr<T>&);                                                \
  }                                                                     \
  }                                                                     \
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Fad::Exp::OP;                                           \
}

UNARYFUNC_MACRO(exp, ExpOp)
UNARYFUNC_MACRO(log, LogOp)
UNARYFUNC_MACRO(log10, Log10Op)
UNARYFUNC_MACRO(sqrt, SqrtOp)
UNARYFUNC_MACRO(safe_sqrt, SafeSqrtOp)
UNARYFUNC_MACRO(cos, CosOp)
UNARYFUNC_MACRO(sin, SinOp)
UNARYFUNC_MACRO(tan, TanOp)
UNARYFUNC_MACRO(acos, ACosOp)
UNARYFUNC_MACRO(asin, ASinOp)
UNARYFUNC_MACRO(atan, ATanOp)
UNARYFUNC_MACRO(cosh, CoshOp)
UNARYFUNC_MACRO(sinh, SinhOp)
UNARYFUNC_MACRO(tanh, TanhOp)
UNARYFUNC_MACRO(acosh, ACoshOp)
UNARYFUNC_MACRO(asinh, ASinhOp)
UNARYFUNC_MACRO(atanh, ATanhOp)
UNARYFUNC_MACRO(abs, AbsOp)
UNARYFUNC_MACRO(fabs, FAbsOp)
UNARYFUNC_MACRO(cbrt, CbrtOp)

#undef UNARYFUNC_MACRO

#define BINARYFUNC_MACRO(OP,FADOP)                                      \
namespace Sacado {                                                      \
                                                                        \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
    template <typename T> class Expr;                                   \
    template <typename T> struct IsFadExpr;                             \
    template <typename T> struct ExprLevel;                             \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_EXPR_EXPR(FADOP)                           \
    OP (const T1&, const T2&);                                          \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    FADOP< typename T::value_type, typename Expr<T>::derived_type,      \
           true, false, typename T::expr_spec_type >                    \
    OP (const typename T::value_type&, const Expr<T>&);                 \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    FADOP< typename Expr<T>::derived_type, typename T::value_type,      \
           false, true, typename T::expr_spec_type >                    \
    OP (const Expr<T>&, const typename T::value_type&);                 \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_SCALAR_EXPR(FADOP)                         \
    OP (const typename T::scalar_type&, const Expr<T>&);                \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_EXPR_SCALAR(FADOP)                         \
    OP (const Expr<T>&, const typename T::scalar_type&);                \
  }                                                                     \
  }                                                                     \
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Fad::Exp::OP;                                           \
}

BINARYFUNC_MACRO(atan2, Atan2Op)
BINARYFUNC_MACRO(pow, PowerOp)
BINARYFUNC_MACRO(max, MaxOp)
BINARYFUNC_MACRO(min, MinOp)

#undef BINARYFUNC_MACRO

#if defined(HAVE_SACADO_KOKKOS)

namespace Sacado {
  namespace Fad {
  namespace Exp {
    template <typename S> class GeneralFad;
    template <typename ValT, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr;

    template <typename S>
    SACADO_INLINE_FUNCTION
    void atomic_add(GeneralFad<S>* dst, const GeneralFad<S>& x);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x);

    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_max_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_max_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_min_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_min_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_add_fetch(GeneralFad<S>* dst, const GeneralFad<S>& x);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_add_fetch(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_sub_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_sub_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_mul_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_mul_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_div_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_div_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);

    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_max(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_max(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_min(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_min(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_add(GeneralFad<S>* dst, const GeneralFad<S>& x);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_sub(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_sub(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_mul(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_mul(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_div(GeneralFad<S>* dest, const GeneralFad<S>& val);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_div(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val);
  }
  }
}

namespace Kokkos {
  using Sacado::Fad::Exp::atomic_add;
  using Sacado::Fad::Exp::atomic_max_fetch;
  using Sacado::Fad::Exp::atomic_min_fetch;
  using Sacado::Fad::Exp::atomic_add_fetch;
  using Sacado::Fad::Exp::atomic_sub_fetch;
  using Sacado::Fad::Exp::atomic_mul_fetch;
  using Sacado::Fad::Exp::atomic_div_fetch;
  using Sacado::Fad::Exp::atomic_fetch_max;
  using Sacado::Fad::Exp::atomic_fetch_min;
  using Sacado::Fad::Exp::atomic_fetch_add;
  using Sacado::Fad::Exp::atomic_fetch_sub;
  using Sacado::Fad::Exp::atomic_fetch_mul;
  using Sacado::Fad::Exp::atomic_fetch_div;
}

#endif

#endif // SACADO_FAD_EXP_MATHFUNCTIONS_HPP
