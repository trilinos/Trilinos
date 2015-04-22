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
#ifndef KOKKOS_BLAS1_MV_IMPL_NRMINF_HPP_
#define KOKKOS_BLAS1_MV_IMPL_NRMINF_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Inf-norm functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class SizeType = typename XV::size_type>
struct V_NrmInf_Functor
{
  typedef typename XV::execution_space              execution_space;
  typedef SizeType                                        size_type;
  typedef typename XV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                         value_type;

  RV m_r;
  XV m_x;

  V_NrmInf_Functor (const RV& r, const XV& x) :
    m_r (r), m_x (x)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "V_NrmInf_Functor: R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value, "KokkosBlas::Impl::"
                   "V_NrmInf_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_NrmInf_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (RV::rank == 0 && XV::rank == 1,
                   "KokkosBlas::Impl::V_NrmInf_Functor: "
                   "RV must have rank 0 and XV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type& update) const
  {
    const typename IPT::mag_type tmp = IPT::norm (m_x(i));
    if (update < tmp) {
      update = tmp;
    }
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = AT::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    if (update < source) {
      update = source;
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void final (const value_type& dst) const
  {
    m_r() = dst;
  }
};

/// \brief Inf-norm functor for multivectors.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class SizeType = typename XMV::size_type>
struct MV_NrmInf_Functor {
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                       value_type[];

  const size_type value_count;
  RV norms_;
  XMV X_;

  MV_NrmInf_Functor (const RV& norms, const XMV& X) :
    value_count (X.dimension_1 ()), norms_ (norms), X_ (X)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value, "KokkosBlas::Impl::"
                   "MV_NrmInf_Functor: R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value, "KokkosBlas::Impl::"
                   "MV_NrmInf_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_NrmInf_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (RV::rank == 1 && XMV::rank == 2,
                   "KokkosBlas::Impl::MV_NrmInf_Functor: "
                   "RV must have rank 1 and XMV must have rank 2.");
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type& i, value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < value_count; ++j) {
      const typename IPT::mag_type tmp = IPT::norm (X_(i,j));
      if (update[j] < tmp) {
        update[j] = tmp;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      update[k] = AT::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      if (update[k] < source[k]) {
        update[k] = source[k];
      }
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      norms_(k) = dst[k];
    }
  }
};


/// \brief Compute the inf-norm of the single vector (1-D View) X, and
///   store the result in the 0-D View r.
template<class RV, class XV, class SizeType>
void
V_NrmInf_Invoke (const RV& r, const XV& X)
{
  typedef typename XV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  typedef V_NrmInf_Functor<RV, XV, SizeType> functor_type;
  functor_type op (r, X);
  Kokkos::parallel_reduce (policy, op);
}


/// \brief Compute the inf-norms of the columns of the multivector
///   (2-D View) X, and store result(s) in the 1-D View r.
template<class RV, class XMV, class SizeType>
void
MV_NrmInf_Invoke (const RV& r, const XMV& X)
{
  typedef typename XMV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  // If the input multivector (2-D View) only has one column, invoke
  // the single-vector version of the kernel.
  if (X.dimension_1 () == 1) {
    auto r_0 = Kokkos::subview (r, 0);
    auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
    typedef decltype (r_0) RV0D;
    typedef decltype (X_0) XV1D;
    V_NrmInf_Invoke<RV0D, XV1D, SizeType> (r_0, X_0);
  }
  else {
    typedef MV_NrmInf_Functor<RV, XMV, SizeType> functor_type;
    functor_type op (r, X);
    Kokkos::parallel_reduce (policy, op);
  }
}


/// \brief Implementation of KokkosBlas::nrmInf for multivectors and
///   single vectors.
template<class RV, class XMV, int rank = XMV::rank>
struct NrmInf_MV {};


//! Special case for multivectors (rank-2 Views).
template<class RV, class XMV>
struct NrmInf_MV<RV, XMV, 2> {
  /// \brief Compute the inf-norm(s) of the column(s) of the
  ///   multivector (2-D View) X, and store result(s) in r.
  static void nrmInf (const RV& r, const XMV& X)
  {
    typedef typename XMV::size_type size_type;
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      MV_NrmInf_Invoke<RV, XMV, int> (r, X);
    }
    else {
      MV_NrmInf_Invoke<RV, XMV, size_type> (r, X);
    }
  }
};


//! Special case for single vectors (rank-1 Views).
template<class RV, class XV>
struct NrmInf_MV<RV, XV, 1> {
  /// \brief Compute the inf-norm of the vector (1-D View) X, and
  ///   store the result in the 0-D View r.
  static void nrmInf (const RV& r, const XV& X)
  {
    typedef typename XV::size_type size_type;
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      V_NrmInf_Invoke<RV, XV, int> (r, X);
    }
    else {
      V_NrmInf_Invoke<RV, XV, size_type> (r, X);
    }
  }
};

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::NrmInf_MV for rank == 2.  This is NOT for
// users!!!  All the declarations of full specializations go in this
// header file.  We may spread out definitions (see _DEF macro below)
// across one or more .cpp files.
//

#define KOKKOSBLAS_IMPL_MV_NRMINF_RANK2_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template<> \
struct NrmInf_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type*, \
                            EXEC_SPACE::array_layout, \
                            Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                            Kokkos::Impl::ViewDefault>, \
               Kokkos::View<const SCALAR**, \
                            LAYOUT, \
                            Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                            Kokkos::Impl::ViewDefault>, \
               2> \
{ \
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type*, \
                       EXEC_SPACE::array_layout, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                       Kokkos::Impl::ViewDefault> RV; \
  typedef Kokkos::View<const SCALAR**, \
                       LAYOUT, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                       Kokkos::Impl::ViewDefault> XMV; \
  static void nrmInf (const RV& r, const XMV& X); \
};

//
// Declarations of full specializations of Impl::NrmInf_MV for rank == 2.
// Their definitions go in .cpp file(s) in this source directory.
//

#ifdef KOKKOS_HAVE_SERIAL

KOKKOSBLAS_IMPL_MV_NRMINF_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP

KOKKOSBLAS_IMPL_MV_NRMINF_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD

KOKKOSBLAS_IMPL_MV_NRMINF_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA

KOKKOSBLAS_IMPL_MV_NRMINF_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace )

#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA

KOKKOSBLAS_IMPL_MV_NRMINF_RANK2_DECL( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )

#endif // KOKKOS_HAVE_CUDA

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::NrmInf_MV for rank == 2.  This is NOT for users!!!
//

#define KOKKOSBLAS_IMPL_MV_NRMINF_RANK2_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
void \
NrmInf_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<SCALAR>::mag_type*, \
                       EXEC_SPACE::array_layout, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                       Kokkos::Impl::ViewDefault>, \
          Kokkos::View<const SCALAR**, \
                       LAYOUT, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                       Kokkos::Impl::ViewDefault>, \
          2>:: \
nrmInf (const RV& r, const XMV& X) \
{ \
  typedef XMV::size_type size_type; \
  const size_type numRows = X.dimension_0 (); \
  const size_type numCols = X.dimension_1 (); \
 \
  if (numRows < static_cast<size_type> (INT_MAX) && \
      numRows * numCols < static_cast<size_type> (INT_MAX)) { \
    MV_NrmInf_Invoke<RV, XMV, int> (r, X); \
  } \
  else { \
    MV_NrmInf_Invoke<RV, XMV, size_type> (r, X); \
  } \
}

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_NRMINF_HPP_
