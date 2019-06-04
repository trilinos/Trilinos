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

namespace Sacado {
  namespace MP {

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<a.size(); ++i)
        c.fastAccessCoeff(i) =
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);
      return c;
    }

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<a.size(); ++i)
        c.fastAccessCoeff(i) =
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);
      return c;
    }

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<a.size(); ++i)
        c.fastAccessCoeff(i) =
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);
      return c;
    }

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<a.size(); ++i)
        c.fastAccessCoeff(i) =
          a.fastAccessCoeff(i) OPER b.fastAccessCoeff(i);
      return c;
    }

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,
            const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<b.size(); ++i)
        c.fastAccessCoeff(i) =
          a OPER b.fastAccessCoeff(i);
      return c;
    }

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& a,
            const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<b.size(); ++i)
        c.fastAccessCoeff(i) =
          a OPER b.fastAccessCoeff(i);
      return c;
    }

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,
            const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<a.size(); ++i)
        c.fastAccessCoeff(i) =
          a.fastAccessCoeff(i) OPER b;
      return c;
    }

    template <typename O, typename T, int N, typename D>
    KOKKOS_INLINE_FUNCTION
    Vector< Stokhos::StaticFixedStorage<O,T,N,D> >
    OPNAME (const volatile Vector< Stokhos::StaticFixedStorage<O,T,N,D> >& a,
            const typename Vector< Stokhos::StaticFixedStorage<O,T,N,D> >::value_type& b)
    {
      Vector< Stokhos::StaticFixedStorage<O,T,N,D> > c;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (O i=0; i<a.size(); ++i)
        c.fastAccessCoeff(i) =
          a.fastAccessCoeff(i) OPER b;
      return c;
    }
  }
}
