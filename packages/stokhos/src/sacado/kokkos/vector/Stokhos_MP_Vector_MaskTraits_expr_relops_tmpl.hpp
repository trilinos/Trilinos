// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
