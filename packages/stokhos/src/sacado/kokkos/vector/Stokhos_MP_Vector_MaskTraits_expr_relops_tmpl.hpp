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

    template <typename V, typename V2>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const Expr<V> &a1,
                     const Expr<V2> &a2)
    {
      const V& v1 = a1.derived();
      const V2& v2 = a2.derived();
      Mask<V> mask;
      if (v2.hasFastAccess(v1.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.fastAccessCoeff(i));
      }
      else{
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.coeff(i));
      }
      return mask;
    }
    
    template <typename V, typename V2>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const volatile Expr<V> &a1,
                     const volatile Expr<V2> &a2)
    {
      const volatile V& v1 = a1.derived();
      const volatile V2& v2 = a2.derived();
      Mask<V> mask;
      if (v2.hasFastAccess(v1.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.fastAccessCoeff(i));
      }
      else{
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.coeff(i));
      }
      return mask;
    }
    
    template <typename V, typename V2>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const Expr<V> &a1,
                     const volatile Expr<V2> &a2)
    {
      const V& v1 = a1.derived();
      const volatile V2& v2 = a2.derived();
      Mask<V> mask;
      if (v2.hasFastAccess(v1.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.fastAccessCoeff(i));
      }
      else{
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.coeff(i));
      }
      return mask;
    }
    
    template <typename V, typename V2>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const volatile Expr<V> &a1,
                     const Expr<V2> &a2)
    {
      const volatile V& v1 = a1.derived();
      const V2& v2 = a2.derived();
      Mask<V> mask;
      if (v2.hasFastAccess(v1.size())) {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.fastAccessCoeff(i));
      }
      else{
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<mask.size; ++i)
          mask.set(i, v1.fastAccessCoeff(i) OPNAME v2.coeff(i));
      }
      return mask;
    }
    
    template <typename V>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const Expr<V> &a1,
                     const typename V::value_type &a2)
    {
      const V& v1 = a1.derived();
      Mask<V> mask;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for(std::size_t i=0; i<mask.size; ++i)
        mask.set(i, v1.fastAccessCoeff(i) OPNAME a2);
      return mask;
    }
    
    template <typename V>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const volatile Expr<V> &a1,
                     const typename V::value_type &a2)
    {
      const volatile V& v1 = a1.derived();
      Mask<V> mask;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for(std::size_t i=0; i<mask.size; ++i)
        mask.set(i, v1.fastAccessCoeff(i) OPNAME a2);
      return mask;
    }
    
    template <typename V>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const typename V::value_type &a1,
                     const Expr<V> &a2)
    {
      const V& v2 = a2.derived();
      Mask<V> mask;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for(std::size_t i=0; i<mask.size; ++i)
        mask.set(i, a1 OPNAME v2.fastAccessCoeff(i));
      return mask;
    }
    
    template <typename V>
    KOKKOS_INLINE_FUNCTION
    Mask<V>
    operator OPNAME (const typename V::value_type &a1,
                     const volatile Expr<V> &a2)
    {
      const volatile V& v2 = a2.derived();
      Mask<V> mask;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for(std::size_t i=0; i<mask.size; ++i)
        mask.set(i, a1 OPNAME v2.fastAccessCoeff(i));
      return mask;
    }
  }
}
