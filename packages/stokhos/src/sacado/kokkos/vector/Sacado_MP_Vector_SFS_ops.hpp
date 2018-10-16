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
