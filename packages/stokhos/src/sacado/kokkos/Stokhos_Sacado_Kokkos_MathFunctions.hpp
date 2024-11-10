// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SACADO_KOKKOS_MATHFUNCTIONS_HPP
#define STOKHOS_SACADO_KOKKOS_MATHFUNCTIONS_HPP

#include <cmath>

#include "Stokhos_ConfigDefs.h"

#include "Kokkos_Macros.hpp"

#ifdef HAVE_STOKHOS_ENSEMBLE_SCALAR_TYPE

#if STOKHOS_USE_MP_VECTOR_SFS_SPEC
#define UNARYFUNC_MACRO_SFS(OP,FADOP)                                   \
namespace Stokhos {                                                     \
  template <typename O, typename T, int N, typename D>                  \
  class StaticFixedStorage;                                             \
}                                                                       \
namespace Sacado {                                                      \
  namespace MP {                                                        \
    template <typename S> class Vector;                                 \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OP (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >&);         \
  }                                                                     \
}
#else
#define UNARYFUNC_MACRO_SFS(OP,FADOP) /* */
#endif

#define UNARYFUNC_MACRO(OP,FADOP)                                       \
UNARYFUNC_MACRO_SFS(OP,FADOP)                                           \
namespace Sacado {                                                      \
                                                                        \
  namespace MP {                                                        \
    template <typename T> class FADOP;                                  \
    template <typename T> class Expr;                                   \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    FADOP< T >                                                          \
    OP (const Expr<T>&);                                                \
  }                                                                     \
                                                                        \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::MP::OP;                                                 \
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
UNARYFUNC_MACRO(ceil, CeilOp)

#undef UNARYFUNC_MACRO
#undef UNARYFUNC_MACRO_SFS

#if STOKHOS_USE_MP_VECTOR_SFS_SPEC
#define BINARYFUNC_MACRO_SFS(OP,FADOP)                                  \
namespace Stokhos {                                                     \
  template <typename O, typename T, int N, typename D>                  \
  class StaticFixedStorage;                                             \
}                                                                       \
namespace Sacado {                                                      \
  namespace MP {                                                        \
    template <typename S> class Vector;                                 \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OP (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >&,          \
        const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >&);         \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OP (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type&, \
        const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >&);         \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OP (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >&,          \
        const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type&); \
  }                                                                     \
}
#else
#define BINARYFUNC_MACRO_SFS(OP,FADOP) /* */
#endif

#define BINARYFUNC_MACRO(OP,FADOP)                                      \
BINARYFUNC_MACRO_SFS(OP,FADOP)                                          \
namespace Sacado {                                                      \
                                                                        \
  namespace MP {                                                        \
    template <typename T1, typename T2> class FADOP;                    \
    template <typename T> class Expr;                                   \
                                                                        \
    template <typename T1, typename T2>                                 \
    KOKKOS_INLINE_FUNCTION                                              \
    FADOP< T1, T2 >                                                     \
    OP (const Expr<T1>&,                                                \
        const Expr<T2>&);                                               \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    FADOP< typename T::value_type, T >                                  \
    OP (const typename T::value_type&,                                  \
        const Expr<T>&);                                                \
                                                                        \
    template <typename T>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    FADOP< T, typename T::value_type >                                  \
    OP (const Expr<T>&,                                                 \
        const typename T::value_type&);                                 \
  }                                                                     \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::MP::OP;                                                 \
}

BINARYFUNC_MACRO(atan2, Atan2Op)
BINARYFUNC_MACRO(pow, PowerOp)
BINARYFUNC_MACRO(max, MaxOp)
BINARYFUNC_MACRO(min, MinOp)
BINARYFUNC_MACRO(fmax, FMaxOp)
BINARYFUNC_MACRO(fmin, FMinOp)

#undef BINARYFUNC_MACRO
#undef BINARYFUNC_MACRO_SFS

#endif

#ifdef HAVE_STOKHOS_PCE_SCALAR_TYPE

#define UNARYFUNC_MACRO(OP,FADOP)                                       \
namespace Sacado {                                                      \
                                                                        \
  namespace UQ {                                                        \
    template <typename S> class PCE;                                    \
    template <typename S>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    PCE<S> OP (const PCE<S>&);                                          \
  }                                                                     \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::UQ::OP;                                                 \
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
UNARYFUNC_MACRO(ceil, CeilOp)

#undef UNARYFUNC_MACRO

#define BINARYFUNC_MACRO(OP)                                            \
namespace Sacado {                                                      \
                                                                        \
  namespace UQ {                                                        \
    template <typename S> class PCE;                                    \
    template <typename S>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    PCE<S> OP (const PCE<S>&, const PCE<S>&);                           \
    template <typename S>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    PCE<S> OP (const typename PCE<S>::value_type&, const PCE<S>&);      \
    template <typename S>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    PCE<S> OP (const PCE<S>&, const typename PCE<S>::value_type&);      \
  }                                                                     \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::UQ::OP;                                                 \
}

BINARYFUNC_MACRO(atan2)
BINARYFUNC_MACRO(pow)

#undef BINARYFUNC_MACRO

#define BINARYFUNC_MACRO(OP)                                            \
namespace Sacado {                                                      \
                                                                        \
  namespace UQ {                                                        \
    template <typename S> class PCE;                                    \
    template <typename S>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    PCE<S> OP (const typename S::value_type&, const PCE<S>&);           \
    template <typename S>                                               \
    KOKKOS_INLINE_FUNCTION                                              \
    PCE<S> OP (const PCE<S>&, const typename S::value_type&);           \
  }                                                                     \
}                                                                       \
                                                                        \
namespace std {                                                         \
  using Sacado::UQ::OP;                                                 \
}

BINARYFUNC_MACRO(max)
BINARYFUNC_MACRO(min)

#undef BINARYFUNC_MACRO

#endif

#endif // STOKHOS_SACADO_KOKKOS_MATHFUNCTIONS_HPP
