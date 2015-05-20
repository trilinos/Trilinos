/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOS_BLAS1_MV_IMPL_FILL_HPP_
#define KOKKOS_BLAS1_MV_IMPL_FILL_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

template<class XMV, class SizeType = typename XMV::size_type>
struct MV_FillFunctor {
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;

  const size_type numCols_;
  const xvalue_type val_;
  XMV X_;

  MV_FillFunctor (const XMV& X, const xvalue_type& val) :
    numCols_ (X.dimension_1 ()), val_ (val), X_ (X)
  {
    static_assert (Kokkos::Impl::is_view<XMV>::value,
                   "KokkosBlas::Impl::MV_FillFunctor: "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename XMV::value_type,
                   typename XMV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_FillFunctor: X is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (XMV::rank == 2, "KokkosBlas::Impl::MV_FillFunctor: "
                   "XMV must have rank 2.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    for (size_type j = 0; j < numCols_; ++j) {
      X_(i,j) = val_;
    }
  }
};

template<class XV, class SizeType = typename XV::size_type>
struct V_FillFunctor {
  typedef typename XV::execution_space            execution_space;
  typedef SizeType                                      size_type;
  typedef typename XV::non_const_value_type           xvalue_type;

  const xvalue_type val_;
  XV x_;

  V_FillFunctor (const XV& x, const xvalue_type& val) :
    val_ (val), x_ (x)
  {
    static_assert (Kokkos::Impl::is_view<XV>::value,
                   "KokkosBlas::Impl::V_FillFunctor: "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename XV::value_type,
                   typename XV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_FillFunctor: X is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (XV::rank == 1, "KokkosBlas::Impl::V_FillFunctor: "
                   "XV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    x_(i) = val_;
  }
};

template<class XV, class SizeType>
void
V_Fill_Invoke (const XV& X, const typename XV::non_const_value_type& val)
{
  typedef typename XV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  typedef V_FillFunctor<XV, SizeType> functor_type;
  functor_type op (X, val);
  Kokkos::parallel_for (policy, op);
}

template<class XMV, class SizeType>
void
MV_Fill_Invoke (const XMV& X, const typename XMV::non_const_value_type& val)
{
  typedef typename XMV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  // Special case for a multivector (2-D View) with a single column.
  // In that case, use the single-vector version of the kernel.
  if (X.dimension_1 () == 1) {
    auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
    typedef decltype (X_0) XV1D;
    V_Fill_Invoke<XV1D, SizeType> (X_0, val);
  }
  else {
    typedef MV_FillFunctor<XMV, SizeType> functor_type;
    functor_type op (X, val);
    Kokkos::parallel_for (policy, op);
  }
}

/// \brief Implementation of KokkosBlas::fill for multivectors and
///   single vectors.
template<class XMV, int rank = XMV::rank>
struct Fill {};

// Specialization for multivectors (2-D Views).
template<class XMV>
struct Fill<XMV, 2> {
  static void fill (const XMV& X, const typename XMV::non_const_value_type& val)
  {
    typedef typename XMV::size_type size_type;
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // The first condition helps avoid overflow with the
    // multiplication in the second condition.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      MV_Fill_Invoke<XMV, int> (X, val);
    }
    else {
      MV_Fill_Invoke<XMV, size_type> (X, val);
    }
  }
};

// Specialization for single vectors (1-D Views).
template<class XV>
struct Fill<XV, 1> {
  static void fill (const XV& X, const typename XV::non_const_value_type& val)
  {
    typedef typename XV::size_type size_type;
    const size_type numRows = X.dimension_0 ();

    if (numRows < static_cast<size_type> (INT_MAX)) {
      V_Fill_Invoke<XV, int> (X, val);
    }
    else {
      V_Fill_Invoke<XV, size_type> (X, val);
    }
  }
};

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Fill for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS_IMPL_MV_FILL_RANK2_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template<> \
struct Fill<Kokkos::View<SCALAR**, \
                         LAYOUT, \
                         Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                         Kokkos::Impl::ViewDefault>, \
            2> \
{ \
  typedef Kokkos::View<SCALAR**, \
                       LAYOUT, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                       Kokkos::Impl::ViewDefault> XMV; \
  static void fill (const XMV& X, const XMV::non_const_value_type& val); \
};

//
// Declarations of full specializations of Impl::Fill for rank == 2.
// Their definitions go in .cpp file(s) in this source directory.
//

#ifdef KOKKOS_HAVE_SERIAL

KOKKOSBLAS_IMPL_MV_FILL_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP

KOKKOSBLAS_IMPL_MV_FILL_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD

KOKKOSBLAS_IMPL_MV_FILL_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA

KOKKOSBLAS_IMPL_MV_FILL_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace )

#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA

KOKKOSBLAS_IMPL_MV_FILL_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )

#endif // KOKKOS_HAVE_CUDA

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Fill for rank == 2.  This is NOT for users!!!  We
// may spread out use of this macro across one or more .cpp files in
// this directory.
//

#define KOKKOSBLAS_IMPL_MV_FILL_RANK2_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
void \
Fill<Kokkos::View<SCALAR**, \
                  LAYOUT, \
                  Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                  Kokkos::Impl::ViewDefault>, \
                  2>:: \
fill (const XMV& X, const XMV::non_const_value_type& val) \
{ \
  typedef XMV::size_type size_type; \
  const size_type numRows = X.dimension_0 (); \
  const size_type numCols = X.dimension_1 (); \
 \
  if (numRows < static_cast<size_type> (INT_MAX) && \
      numRows * numCols < static_cast<size_type> (INT_MAX)) { \
    MV_Fill_Invoke<XMV, int> (X, val); \
  } \
  else { \
    MV_Fill_Invoke<XMV, size_type> (X, val); \
  } \
}


} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_FILL_HPP_
