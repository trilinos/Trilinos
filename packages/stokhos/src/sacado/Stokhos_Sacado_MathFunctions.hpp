// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SACADO_MATH_FUNCTIONS_HPP
#define STOKHOS_SACADO_MATH_FUNCTIONS_HPP

#define UNARYFUNC_MACRO(OP,FADOP)                                       \
namespace Sacado {                                                      \
                                                                        \
  namespace PCE {                                                       \
    template <typename T, typename S> class OrthogPoly;                 \
    template <typename T, typename S>                                   \
    OrthogPoly<T,S> OP (const OrthogPoly<T,S>&);                        \
  }                                                                     \
                                                                        \
  namespace ETPCE {                                                     \
    template <typename T> class FADOP;                                  \
    template <typename T> class Expr;                                   \
    template <typename T>                                               \
    Expr< FADOP< Expr<T> > > OP (const Expr<T>&);                       \
  }                                                                     \
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::PCE::OP;                                                \
  using Sacado::ETPCE::OP;                                              \
}

UNARYFUNC_MACRO(exp, ExpOp)
UNARYFUNC_MACRO(log, LogOp)
UNARYFUNC_MACRO(log10, Log10Op)
UNARYFUNC_MACRO(sqrt, SqrtOp)
UNARYFUNC_MACRO(cbrt, CbrtOp)
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

#undef UNARYFUNC_MACRO

#define BINARYFUNC_MACRO(OP,FADOP)                                      \
namespace Sacado {                                                      \
                                                                        \
  namespace PCE {                                                       \
    template <typename T, typename S> class OrthogPoly;                 \
    template <typename T, typename S>                                   \
    OrthogPoly<T,S> OP (const OrthogPoly<T,S>&,                         \
                        const OrthogPoly<T,S>&);                        \
    template <typename T, typename S>                                   \
    OrthogPoly<T,S> OP (const T&,                                       \
                        const OrthogPoly<T,S>&);                        \
    template <typename T, typename S>                                   \
    OrthogPoly<T,S> OP (const OrthogPoly<T,S>&,                         \
                        const T&);                                      \
  }                                                                     \
                                                                        \
  namespace ETPCE {                                                     \
    template <typename T1, typename T2> class FADOP;                    \
    template <typename T> class Expr;                                   \
    template <typename T> class ConstExpr;                              \
    template <typename T1, typename T2>                                 \
    Expr< FADOP< Expr<T1>, Expr<T2> > >                                 \
    OP (const Expr<T1>&, const Expr<T2>&);                              \
                                                                        \
    template <typename T>                                               \
    Expr< FADOP< Expr<T>, Expr<T> > >                                   \
    OP (const Expr<T>&, const Expr<T>&);                                \
                                                                        \
    template <typename T>                                               \
    Expr< FADOP< typename Expr<T>::value_type, Expr<T> > >              \
    OP (const typename Expr<T>::value_type&, const Expr<T>&);           \
                                                                        \
    template <typename T>                                               \
    Expr< FADOP< Expr<T>, typename Expr<T>::value_type > >              \
    OP (const Expr<T>&, const typename Expr<T>::value_type&);           \
  }                                                                     \
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::PCE::OP;                                                \
  using Sacado::ETPCE::OP;                                              \
}

BINARYFUNC_MACRO(atan2, Atan2Op)
BINARYFUNC_MACRO(pow, PowerOp)
BINARYFUNC_MACRO(max, MaxOp)
BINARYFUNC_MACRO(min, MinOp)

#undef BINARYFUNC_MACRO

#endif // STOKHOS_SACADO_MATH_FUNCTIONS_HPP
