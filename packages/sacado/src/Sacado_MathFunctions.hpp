// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MATHFUNCTIONS_HPP
#define SACADO_MATHFUNCTIONS_HPP

#include <cmath>

#include "Sacado_ConfigDefs.h"
#include "Sacado_Base.hpp"
#include "Sacado_cmath.hpp"
#include "Sacado_SFINAE_Macros.hpp"

// Note:  Sacado::Fad::Ops are forward-declared here, instead of in macros
// below.
#include "Sacado_Fad_Ops_Fwd.hpp"

#define SACADO_UNARYFUNC_MACRO(OP,FADOP)                                \
namespace Sacado {                                                      \
                                                                        \
  namespace Fad {                                                       \
    template <typename T> class Expr;                                   \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    FADOP< typename Expr<T>::derived_type,                              \
           typename T::expr_spec_type >                                 \
    OP (const Expr<T>&);                                                \
  }                                                                     \
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Fad::OP;                                                \
}

SACADO_UNARYFUNC_MACRO(exp, ExpOp)
SACADO_UNARYFUNC_MACRO(expm1, ExpM1Op)
SACADO_UNARYFUNC_MACRO(log, LogOp)
SACADO_UNARYFUNC_MACRO(log1p, Log1POp)
SACADO_UNARYFUNC_MACRO(log10, Log10Op)
SACADO_UNARYFUNC_MACRO(sqrt, SqrtOp)
SACADO_UNARYFUNC_MACRO(safe_sqrt, SafeSqrtOp)
SACADO_UNARYFUNC_MACRO(cos, CosOp)
SACADO_UNARYFUNC_MACRO(sin, SinOp)
SACADO_UNARYFUNC_MACRO(tan, TanOp)
SACADO_UNARYFUNC_MACRO(acos, ACosOp)
SACADO_UNARYFUNC_MACRO(asin, ASinOp)
SACADO_UNARYFUNC_MACRO(atan, ATanOp)
SACADO_UNARYFUNC_MACRO(cosh, CoshOp)
SACADO_UNARYFUNC_MACRO(sinh, SinhOp)
SACADO_UNARYFUNC_MACRO(tanh, TanhOp)
SACADO_UNARYFUNC_MACRO(acosh, ACoshOp)
SACADO_UNARYFUNC_MACRO(asinh, ASinhOp)
SACADO_UNARYFUNC_MACRO(atanh, ATanhOp)
SACADO_UNARYFUNC_MACRO(abs, AbsOp)
SACADO_UNARYFUNC_MACRO(fabs, FAbsOp)
SACADO_UNARYFUNC_MACRO(cbrt, CbrtOp)

#undef SACADO_UNARYFUNC_MACRO

#define SACADO_BINARYFUNC_MACRO(OP,FADOP)                               \
namespace Sacado {                                                      \
                                                                        \
  namespace Fad {                                                       \
                                                                        \
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
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Fad::OP;                                                \
}

SACADO_BINARYFUNC_MACRO(atan2, Atan2Op)
SACADO_BINARYFUNC_MACRO(pow, PowerOp)
SACADO_BINARYFUNC_MACRO(max, MaxOp)
SACADO_BINARYFUNC_MACRO(min, MinOp)

#undef SACADO_BINARYFUNC_MACRO

#if defined(HAVE_SACADO_KOKKOS)

namespace Sacado {
  namespace Fad {
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

namespace Kokkos {
  using Sacado::Fad::atomic_add;
  using Sacado::Fad::atomic_max_fetch;
  using Sacado::Fad::atomic_min_fetch;
  using Sacado::Fad::atomic_add_fetch;
  using Sacado::Fad::atomic_sub_fetch;
  using Sacado::Fad::atomic_mul_fetch;
  using Sacado::Fad::atomic_div_fetch;
  using Sacado::Fad::atomic_fetch_max;
  using Sacado::Fad::atomic_fetch_min;
  using Sacado::Fad::atomic_fetch_add;
  using Sacado::Fad::atomic_fetch_sub;
  using Sacado::Fad::atomic_fetch_mul;
  using Sacado::Fad::atomic_fetch_div;
}

#endif

#define SACADO_UNARYFUNC_MACRO(OP,FADOP)                                \
namespace Sacado {                                                      \
                                                                        \
  namespace ELRFad {                                                    \
    template <typename T> class FADOP;                                  \
    template <typename T> class Expr;                                   \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T> > > OP (const Expr<T>&);                       \
  }                                                                     \
                                                                        \
  namespace CacheFad {                                                  \
    template <typename T> class FADOP;                                  \
    template <typename T> class Expr;                                   \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T> > > OP (const Expr<T>&);                       \
  }                                                                     \
                                                                        \
  namespace ELRCacheFad {                                               \
    template <typename T> class FADOP;                                  \
    template <typename T> class Expr;                                   \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T> > > OP (const Expr<T>&);                       \
  }                                                                     \
                                                                        \
  namespace LFad {                                                      \
    template <typename T> class FADOP;                                  \
    template <typename T> class Expr;                                   \
    template <typename T>                                               \
    Expr< FADOP< Expr<T> > > OP (const Expr<T>&);                       \
  }                                                                     \
                                                                        \
  namespace Tay {                                                       \
    template <typename T> class Taylor;                                 \
    template <typename T> Taylor<T> OP (const Base< Taylor<T> >&);      \
  }                                                                     \
                                                                        \
  namespace FlopCounterPack {                                           \
    template <typename T> class ScalarFlopCounter;                      \
    template <typename T>                                               \
    ScalarFlopCounter<T> OP (const Base< ScalarFlopCounter<T> >&);      \
  }                                                                     \
                                                                        \
  namespace Rad {                                                       \
    template <typename T> class ADvari;                                 \
    template <typename T> class IndepADvar;                             \
    template <typename T> ADvari<T>& OP (const Base< ADvari<T> >&);     \
    template <typename T> ADvari<T>& OP (const Base< IndepADvar<T> >&); \
  }                                                                     \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Fad::OP;                                                \
  using Sacado::ELRFad::OP;                                             \
  using Sacado::CacheFad::OP;                                           \
  using Sacado::ELRCacheFad::OP;                                        \
  using Sacado::LFad::OP;                                               \
  using Sacado::Tay::OP;                                                \
  using Sacado::FlopCounterPack::OP;                                    \
  using Sacado::Rad::OP;                                                \
}

SACADO_UNARYFUNC_MACRO(exp, ExpOp)
SACADO_UNARYFUNC_MACRO(expm1, ExpM1Op)
SACADO_UNARYFUNC_MACRO(log, LogOp)
SACADO_UNARYFUNC_MACRO(log1p, Log1POp)
SACADO_UNARYFUNC_MACRO(log10, Log10Op)
SACADO_UNARYFUNC_MACRO(sqrt, SqrtOp)
SACADO_UNARYFUNC_MACRO(cos, CosOp)
SACADO_UNARYFUNC_MACRO(sin, SinOp)
SACADO_UNARYFUNC_MACRO(tan, TanOp)
SACADO_UNARYFUNC_MACRO(acos, ACosOp)
SACADO_UNARYFUNC_MACRO(asin, ASinOp)
SACADO_UNARYFUNC_MACRO(atan, ATanOp)
SACADO_UNARYFUNC_MACRO(cosh, CoshOp)
SACADO_UNARYFUNC_MACRO(sinh, SinhOp)
SACADO_UNARYFUNC_MACRO(tanh, TanhOp)
SACADO_UNARYFUNC_MACRO(acosh, ACoshOp)
SACADO_UNARYFUNC_MACRO(asinh, ASinhOp)
SACADO_UNARYFUNC_MACRO(atanh, ATanhOp)
SACADO_UNARYFUNC_MACRO(abs, AbsOp)
SACADO_UNARYFUNC_MACRO(fabs, FAbsOp)
SACADO_UNARYFUNC_MACRO(cbrt, CbrtOp)

#undef SACADO_UNARYFUNC_MACRO

#define SACADO_BINARYFUNC_MACRO(OP,FADOP)                               \
namespace Sacado {                                                      \
                                                                        \
  namespace Fad {                                                       \
    template <typename T> class SimpleFad;                              \
    template <typename T>                                               \
    SimpleFad<T>                                                        \
    OP (const SimpleFad<T>&, const SimpleFad<T>&);                      \
                                                                        \
    template <typename T>                                               \
    SimpleFad<T>                                                        \
    OP (const SimpleFad<T>&,                                            \
        const typename SimpleFad<T>::value_type&);                      \
                                                                        \
    template <typename T>                                               \
    SimpleFad<T>                                                        \
    OP (const typename SimpleFad<T>::value_type&,                       \
        const SimpleFad<T>&);                                           \
  }                                                                     \
                                                                        \
  namespace ELRFad {                                                    \
    template <typename T1, typename T2> class FADOP;                    \
    template <typename T> class Expr;                                   \
    template <typename T> class ConstExpr;                              \
    template <typename T> struct IsFadExpr;                             \
    template <typename T> struct ExprLevel;                             \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_EXPR(FADOP)                               \
    OP (const T1&, const T2&);                                          \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T>, Expr<T> > >                                   \
    OP (const Expr<T>&, const Expr<T>&);                                \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< ConstExpr<typename Expr<T>::value_type>, Expr<T> > >   \
    OP (const typename Expr<T>::value_type&, const Expr<T>&);           \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T>, ConstExpr<typename Expr<T>::value_type> > >   \
    OP (const Expr<T>&, const typename Expr<T>::value_type&);           \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_SCALAR_EXPR(FADOP)                             \
    OP (const typename Expr<T>::scalar_type&, const Expr<T>&);          \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_SCALAR(FADOP)                             \
    OP (const Expr<T>&, const typename Expr<T>::scalar_type&);          \
  }                                                                     \
                                                                        \
  namespace CacheFad {                                                  \
    template <typename T1, typename T2> class FADOP;                    \
    template <typename T> class Expr;                                   \
    template <typename T> class ConstExpr;                              \
    template <typename T> struct IsFadExpr;                             \
    template <typename T> struct ExprLevel;                             \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_EXPR(FADOP)                               \
    OP (const T1&, const T2&);                                          \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T>, Expr<T> > >                                   \
    OP (const Expr<T>&, const Expr<T>&);                                \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< ConstExpr<typename Expr<T>::value_type>, Expr<T> > >   \
    OP (const typename Expr<T>::value_type&, const Expr<T>&);           \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T>, ConstExpr<typename Expr<T>::value_type> > >   \
    OP (const Expr<T>&, const typename Expr<T>::value_type&);           \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_SCALAR_EXPR(FADOP)                             \
    OP (const typename Expr<T>::scalar_type&, const Expr<T>&);          \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_SCALAR(FADOP)                             \
    OP (const Expr<T>&, const typename Expr<T>::scalar_type&);          \
  }                                                                     \
                                                                        \
  namespace ELRCacheFad {                                               \
    template <typename T1, typename T2> class FADOP;                    \
    template <typename T> class Expr;                                   \
    template <typename T> class ConstExpr;                              \
    template <typename T> struct IsFadExpr;                             \
    template <typename T> struct ExprLevel;                             \
    template <typename T1, typename T2>                                 \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_EXPR(FADOP)                               \
    OP (const T1&, const T2&);                                          \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T>, Expr<T> > >                                   \
    OP (const Expr<T>&, const Expr<T>&);                                \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< ConstExpr<typename Expr<T>::value_type>, Expr<T> > >   \
    OP (const typename Expr<T>::value_type&, const Expr<T>&);           \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    Expr< FADOP< Expr<T>, ConstExpr<typename Expr<T>::value_type> > >   \
    OP (const Expr<T>&, const typename Expr<T>::value_type&);           \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_SCALAR_EXPR(FADOP)                             \
    OP (const typename Expr<T>::scalar_type&, const Expr<T>&);          \
                                                                        \
    template <typename T>                                               \
    SACADO_INLINE_FUNCTION                                              \
    SACADO_FAD_OP_ENABLE_EXPR_SCALAR(FADOP)                             \
    OP (const Expr<T>&, const typename Expr<T>::scalar_type&);          \
  }                                                                     \
                                                                        \
  namespace LFad {                                                      \
    template <typename T1, typename T2> class FADOP;                    \
    template <typename T> class Expr;                                   \
                                                                        \
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
  namespace Tay {                                                       \
    template <typename T> class Taylor;                                 \
    template <typename T> Taylor<T> OP (                                \
      const Base< Taylor<T> >&,                                         \
      const Base< Taylor<T> >&);                                        \
    template <typename T> Taylor<T> OP (                                \
      const typename Taylor<T>::value_type&,                            \
      const Base< Taylor<T> >&);                                        \
    template <typename T> Taylor<T> OP (                                \
      const Base< Taylor<T> >&,                                         \
      const typename Taylor<T>::value_type&);                           \
  }                                                                     \
                                                                        \
  namespace FlopCounterPack {                                           \
    template <typename T> class ScalarFlopCounter;                      \
    template <typename T>                                               \
    ScalarFlopCounter<T> OP (                                           \
      const Base< ScalarFlopCounter<T> >&,                              \
      const Base< ScalarFlopCounter<T> >&);                             \
    template <typename T>                                               \
    ScalarFlopCounter<T> OP (                                           \
      const typename ScalarFlopCounter<T>::value_type&,                 \
      const Base< ScalarFlopCounter<T> >&);                             \
    template <typename T>                                               \
    ScalarFlopCounter<T> OP (                                           \
      const Base< ScalarFlopCounter<T> >&,                              \
      const typename ScalarFlopCounter<T>::value_type&);                \
    template <typename T>                                               \
    ScalarFlopCounter<T> OP (                                           \
      const int&,                                                       \
      const Base< ScalarFlopCounter<T> >&);                             \
    template <typename T>                                               \
    ScalarFlopCounter<T> OP (                                           \
      const Base< ScalarFlopCounter<T> >&,                              \
      const int&);                                                      \
  }                                                                     \
                                                                        \
  namespace Rad {                                                       \
    template <typename T> class ADvari;                                 \
    template <typename T> class IndepADvar;                             \
    template <typename T> class DoubleAvoid;                            \
    template <typename T> ADvari<T>& OP (const Base< ADvari<T> >&,      \
                                         const Base< ADvari<T> >&);     \
    template <typename T> ADvari<T>& OP (const Base< IndepADvar<T> >&,  \
                                         const Base< ADvari<T> >&);     \
    template <typename T> ADvari<T>& OP (T,                             \
                                         const Base< ADvari<T> >&);     \
    template <typename T> ADvari<T>& OP (typename DoubleAvoid<T>::dtype,\
                                         const Base< ADvari<T> >&);     \
    template <typename T> ADvari<T>& OP (typename DoubleAvoid<T>::itype,\
                                         const Base< ADvari<T> >&);     \
    template <typename T> ADvari<T>& OP (typename DoubleAvoid<T>::ltype,\
                                         const Base< ADvari<T> >&);     \
    template <typename T> ADvari<T>& OP (const Base< ADvari<T> >&,      \
                                         const Base< IndepADvar<T> >&); \
    template <typename T> ADvari<T>& OP (const Base< ADvari<T> >&,      \
                                         T);                            \
    template <typename T> ADvari<T>& OP (const Base< ADvari<T> >&,      \
                                         typename DoubleAvoid<T>::dtype);\
    template <typename T> ADvari<T>& OP (const Base< ADvari<T> >&,      \
                                         typename DoubleAvoid<T>::itype);\
    template <typename T> ADvari<T>& OP (const Base< ADvari<T> >&,      \
                                         typename DoubleAvoid<T>::ltype);\
    template <typename T> ADvari<T>& OP (const Base< IndepADvar<T> >&,  \
                                         const Base< IndepADvar<T> >&); \
    template <typename T> ADvari<T>& OP (T,                             \
                                         const Base< IndepADvar<T> >&); \
    template <typename T> ADvari<T>& OP (typename DoubleAvoid<T>::dtype,\
                                         const Base< IndepADvar<T> >&); \
    template <typename T> ADvari<T>& OP (typename DoubleAvoid<T>::itype,\
                                         const Base< IndepADvar<T> >&); \
    template <typename T> ADvari<T>& OP (typename DoubleAvoid<T>::ltype,\
                                         const Base< IndepADvar<T> >&); \
    template <typename T> ADvari<T>& OP (const Base< IndepADvar<T> >&,  \
                                         T);                            \
    template <typename T> ADvari<T>& OP (const Base< IndepADvar<T> >&,  \
                                         typename DoubleAvoid<T>::dtype);\
    template <typename T> ADvari<T>& OP (const Base< IndepADvar<T> >&,  \
                                         typename DoubleAvoid<T>::itype);\
    template <typename T> ADvari<T>& OP (const Base< IndepADvar<T> >&,  \
                                         typename DoubleAvoid<T>::ltype);\
  }                                                                     \
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::Fad::OP;                                                \
  using Sacado::ELRFad::OP;                                             \
  using Sacado::CacheFad::OP;                                           \
  using Sacado::ELRCacheFad::OP;                                        \
  using Sacado::LFad::OP;                                               \
  using Sacado::Tay::OP;                                                \
  using Sacado::FlopCounterPack::OP;                                    \
  using Sacado::Rad::OP;                                                \
}

SACADO_BINARYFUNC_MACRO(atan2, Atan2Op)
SACADO_BINARYFUNC_MACRO(pow, PowerOp)
SACADO_BINARYFUNC_MACRO(max, MaxOp)
SACADO_BINARYFUNC_MACRO(min, MinOp)

#undef SACADO_BINARYFUNC_MACRO

#if defined(HAVE_SACADO_KOKKOS)

namespace Sacado {
  namespace ELRFad {
    template <typename ValT, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr;
    template <typename T> class DFad;
    template <typename T, int N> class SFad;
    template <typename T, int N> class SLFad;
    template <typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(DFad<T>* dst, const DFad<T>& x);
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SFad<T,N>* dst, const SFad<T,N>& x);
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SLFad<T,N>* dst, const SLFad<T,N>& x);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x);
  }
  namespace CacheFad {
    template <typename ValT, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr;
    template <typename T> class DFad;
    template <typename T, int N> class SFad;
    template <typename T, int N> class SLFad;
    template <typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(DFad<T>* dst, const DFad<T>& x);
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SFad<T,N>* dst, const SFad<T,N>& x);
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SLFad<T,N>* dst, const SLFad<T,N>& x);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x);
  }
  namespace ELRCacheFad {
    template <typename ValT, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr;
    template <typename T> class DFad;
    template <typename T, int N> class SFad;
    template <typename T, int N> class SLFad;
    template <typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(DFad<T>* dst, const DFad<T>& x);
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SFad<T,N>* dst, const SFad<T,N>& x);
    template <typename T, int N>
    SACADO_INLINE_FUNCTION
    void atomic_add(SLFad<T,N>* dst, const SLFad<T,N>& x);
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& x);
  }
}

namespace Kokkos {
  using Sacado::ELRFad::atomic_add;
  using Sacado::CacheFad::atomic_add;
  using Sacado::ELRCacheFad::atomic_add;
}

#endif

#endif // SACADO_MATHFUNCTIONS_HPP
