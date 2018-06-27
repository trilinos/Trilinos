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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_MATHFUNCTIONS_HPP
#define SACADO_FAD_EXP_MATHFUNCTIONS_HPP

#include "Sacado_cmath.hpp"
#include "Sacado_SFINAE_Macros.hpp"

#define UNARYFUNC_MACRO(OP,FADOP)                                       \
namespace Sacado {                                                      \
                                                                        \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
    template <typename T> class Expr;                                   \
    template <typename T, typename E> class FADOP;                      \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
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
    template <typename T1, typename T2, bool, bool, typename E>         \
    class FADOP;                                                        \
    template <typename T> struct IsFadExpr;                             \
    template <typename T> struct ExprLevel;                             \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_EXPR_EXPR(FADOP)                           \
    OP (const T1&, const T2&);                                          \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    FADOP< typename T::value_type, typename Expr<T>::derived_type,      \
           true, false, typename T::expr_spec_type >                    \
    OP (const typename T::value_type&, const Expr<T>&);                 \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    FADOP< typename Expr<T>::derived_type, typename T::value_type,      \
           false, true, typename T::expr_spec_type >                    \
    OP (const Expr<T>&, const typename T::value_type&);                 \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    SACADO_FAD_EXP_OP_ENABLE_SCALAR_EXPR(FADOP)                         \
    OP (const typename T::scalar_type&, const Expr<T>&);                \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
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

#if defined(HAVE_SACADO_KOKKOSCORE)

namespace Sacado {
  namespace Fad {
  namespace Exp {
    template <typename S> class GeneralFad;
    template <typename ValT, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr;

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    void atomic_add(GeneralFad<S>* dst, const GeneralFad<S>& x);

    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    KOKKOS_INLINE_FUNCTION
    void atomic_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x);
  }
  }
}

namespace Kokkos {
  using Sacado::Fad::Exp::atomic_add;
}

#endif

#endif // SACADO_FAD_EXP_MATHFUNCTIONS_HPP
