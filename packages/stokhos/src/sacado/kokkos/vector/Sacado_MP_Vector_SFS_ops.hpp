// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Sacado_cmath.hpp"
#include <ostream>      // for std::ostream

#ifdef __CUDACC__
    #include <cuda_runtime_api.h>
    // including math functions via math_functions.h is deprecated in cuda version >= 10.0
    // the deprecation warning indicates to use cuda_runtime_api.h instead
    #if CUDART_VERSION < 10000
        #include <math_functions.h>
    #endif
#endif

#define OPNAME operator+
#define OPER +
#include "Sacado_MP_Vector_SFS_unary_op_tmpl.hpp"
#undef OPNAME
#undef OPER

#define OPNAME operator-
#define OPER -
#include "Sacado_MP_Vector_SFS_unary_op_tmpl.hpp"
#undef OPNAME
#undef OPER

#define OPNAME exp
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME log
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME log10
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME sqrt
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME cbrt
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME cos
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME sin
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME tan
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME acos
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME asin
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME atan
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME cosh
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME sinh
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME tanh
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME acosh
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME asinh
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME atanh
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME abs
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME fabs
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME ceil
#include "Sacado_MP_Vector_SFS_unary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME operator+
#define OPER +
#include "Sacado_MP_Vector_SFS_binary_op_tmpl.hpp"
#undef OPNAME
#undef OPER

#define OPNAME operator-
#define OPER -
#include "Sacado_MP_Vector_SFS_binary_op_tmpl.hpp"
#undef OPNAME
#undef OPER

#define OPNAME operator*
#define OPER *
#include "Sacado_MP_Vector_SFS_binary_op_tmpl.hpp"
#undef OPNAME
#undef OPER

#define OPNAME operator/
#define OPER /
#include "Sacado_MP_Vector_SFS_binary_op_tmpl.hpp"
#undef OPNAME
#undef OPER

#define OPNAME atan2
#include "Sacado_MP_Vector_SFS_binary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME pow
#include "Sacado_MP_Vector_SFS_binary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME max
#include "Sacado_MP_Vector_SFS_binary_func_tmpl.hpp"
#undef OPNAME

#define OPNAME min
#include "Sacado_MP_Vector_SFS_binary_func_tmpl.hpp"
#undef OPNAME

//#ifdef __CUDACC__
//MP_BINARYOP_MACRO(max, ::max)
//MP_BINARYOP_MACRO(min, ::min)
//#else
//MP_BINARYOP_MACRO(max, std::max)
//MP_BINARYOP_MACRO(min, std::min)
//#endif

//-------------------------- Relational Operators -----------------------

#define MP_RELOP_MACRO(OP)                                              \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
                 const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      return a.val() OP b.val();                                        \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,                                            \
                 const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      return a OP b.val();                                              \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
                 const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)                                            \
    {                                                                   \
      return a.val() OP b;                                              \
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

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    bool operator ! (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a)
    {
      return ! a.val();
    }

  } // namespace MP

} // namespace Sacado


//-------------------------- Boolean Operators -----------------------
namespace Sacado {

  namespace MP {

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    bool toBool(const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& x) {
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
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
                 const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      return toBool(a) OP toBool(b);                                    \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,                                            \
                 const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      return a OP toBool(b);                                            \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    bool                                                                \
    operator OP (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
                 const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)                                            \
    {                                                                   \
      return toBool(a) OP b;                                            \
    }                                                                   \
  }                                                                     \
}

PCE_BOOL_MACRO(&&)
PCE_BOOL_MACRO(||)

#undef PCE_BOOL_MACRO
