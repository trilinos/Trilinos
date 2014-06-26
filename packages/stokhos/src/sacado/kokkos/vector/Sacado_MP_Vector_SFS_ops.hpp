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

#define MP_UNARYOP_MACRO(OPNAME,OPER)                                   \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a)    \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) = OPER(a.fastAccessCoeff(i));              \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a) \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) = OPER(a.fastAccessCoeff(i));              \
      return c;                                                         \
    }                                                                   \
                                                                        \
  }                                                                     \
}

MP_UNARYOP_MACRO(operator+, +)
MP_UNARYOP_MACRO(operator-, -)
MP_UNARYOP_MACRO(exp, std::exp)
MP_UNARYOP_MACRO(log, std::log)
MP_UNARYOP_MACRO(log10, std::log10)
MP_UNARYOP_MACRO(sqrt, std::sqrt)
MP_UNARYOP_MACRO(cos, std::cos)
MP_UNARYOP_MACRO(sin, std::sin)
MP_UNARYOP_MACRO(tan, std::tan)
MP_UNARYOP_MACRO(acos, std::acos)
MP_UNARYOP_MACRO(asin, std::asin)
MP_UNARYOP_MACRO(atan, std::atan)
MP_UNARYOP_MACRO(cosh, std::cosh)
MP_UNARYOP_MACRO(sinh, std::sinh)
MP_UNARYOP_MACRO(tanh, std::tanh)
MP_UNARYOP_MACRO(acosh, std::acosh)
MP_UNARYOP_MACRO(asinh, std::asinh)
MP_UNARYOP_MACRO(atanh, std::atanh)
MP_UNARYOP_MACRO(abs, std::abs)
MP_UNARYOP_MACRO(fabs, std::fabs)

#undef MP_UNARYOP_MACRO

#define MP_BINARYOP_MACRO(OPNAME,OPER)                                  \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,    \
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)    \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);               \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);               \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,    \
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);               \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)    \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);               \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,                                                 \
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)    \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<b.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a OPER b.fastAccessCoeff(i);                                  \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,                                                 \
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<b.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a OPER b.fastAccessCoeff(i);                                  \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,    \
            const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)                                                 \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a.fastAccessCoeff(i) OPER b;                                  \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
            const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)                                                 \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          a.fastAccessCoeff(i) OPER b;                                  \
      return c;                                                         \
    }                                                                   \
  }                                                                     \
}

MP_BINARYOP_MACRO(operator+, +)
MP_BINARYOP_MACRO(operator-, -)
MP_BINARYOP_MACRO(operator*, *)
MP_BINARYOP_MACRO(operator/, /)

#undef MP_BINARYOP_MACRO

#define MP_BINARYOP_MACRO(OPNAME,OPER)                                  \
namespace Sacado {                                                      \
  namespace MP {                                                        \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,    \
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)    \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a.fastAccessCoeff(i) ,  b.fastAccessCoeff(i) );         \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a.fastAccessCoeff(i) ,  b.fastAccessCoeff(i) );         \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,    \
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a.fastAccessCoeff(i) ,  b.fastAccessCoeff(i) );         \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)    \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a.fastAccessCoeff(i) ,  b.fastAccessCoeff(i) );         \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,                                                 \
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)    \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<b.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a , b.fastAccessCoeff(i) );                             \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,                                                 \
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b) \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<b.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a , b.fastAccessCoeff(i) );                             \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,    \
            const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)                                                 \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a.fastAccessCoeff(i) , b );                             \
      return c;                                                         \
    }                                                                   \
                                                                        \
    template <typename O, typename T, int N, typename D>                \
    KOKKOS_INLINE_FUNCTION                                              \
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >                      \
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a, \
            const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)                                                 \
    {                                                                   \
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;                 \
      for (O i=0; i<a.size(); ++i)                                      \
        c.fastAccessCoeff(i) =                                          \
          OPER( a.fastAccessCoeff(i) , b );                             \
      return c;                                                         \
    }                                                                   \
  }                                                                     \
}

MP_BINARYOP_MACRO(atan2, std::atan2)
MP_BINARYOP_MACRO(pow, std::pow)
#ifdef __CUDACC__
MP_BINARYOP_MACRO(max, ::max)
MP_BINARYOP_MACRO(min, ::min)
#else
MP_BINARYOP_MACRO(max, std::max)
MP_BINARYOP_MACRO(min, std::min)
#endif

#undef MP_BINARYOP_MACRO

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
