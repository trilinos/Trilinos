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
  
    template <typename S>
    KOKKOS_INLINE_FUNCTION
    Mask<Vector<S> >
    operator OPNAME (const Vector<S> &a1,
                     const Vector<S> &a2)
    {      
      Mask<Vector<S> > mask;
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
        mask.set(i, a1.fastAccessCoeff(i) OPNAME a2.fastAccessCoeff(i));
      return mask;
    }

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    Mask<Vector<S> >
    operator OPNAME (const Vector<S> &a1,
                     const typename S::value_type &a2)
    {
      Mask<Vector<S> > mask;
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
        mask.set(i, a1.fastAccessCoeff(i) OPNAME a2);
      return mask;
    }

    template <typename S>
    KOKKOS_INLINE_FUNCTION
    Mask<Vector<S> >
    operator OPNAME (const typename S::value_type &a1,
                     const Vector<S> &a2)
    {
      Mask<Vector<S> > mask;
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
        mask.set(i, a1 OPNAME a2.fastAccessCoeff(i));
      return mask;
    }
  }
}
