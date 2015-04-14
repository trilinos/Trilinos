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
#ifndef KOKKOS_BLAS1_MV_IMPL_NRM2_HPP_
#define KOKKOS_BLAS1_MV_IMPL_NRM2_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// nrm2_squared
//

/// \brief 2-norm (squared) functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class SizeType = typename XV::size_type>
struct V_Nrm2Squared_Functor
{
  typedef typename XV::execution_space              execution_space;
  typedef SizeType                                        size_type;
  typedef typename XV::non_const_value_type             xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                         value_type;

  RV m_r;
  typename XV::const_type m_x;

  V_Nrm2Squared_Functor (const RV& r, const XV& x) :
    m_r (r), m_x (x)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RV>::value,
                   "KokkosBlas::Impl::V_Nrm2Squared_Functor: "
                   "R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value,
                   "KokkosBlas::Impl::V_Nrm2Squared_Functor: "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_Nrm2Squared_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (RV::rank == 0 && XV::rank == 1,
                   "KokkosBlas::Impl::V_Nrm2Squared_Functor: "
                   "RV must have rank 0 and XV must have rank 1.");
#endif // KOKKOS_HAVE_CXX11
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type& sum) const
  {
    const typename IPT::mag_type tmp = IPT::norm (m_x(i));
    sum += tmp * tmp;
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& update) const
  {
    update = AT::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source;
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type& dst) const
  {
    m_r() = dst;
  }
};

/// \brief Column-wise 2-norm functor for multivectors; works for
///   any layout, but best performance with LayoutRight.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class SizeType = typename XMV::size_type>
struct MV_Nrm2Squared_Right_FunctorVector
{
  typedef typename XMV::execution_space             execution_space;
  typedef SizeType                                        size_type;
  typedef typename XMV::non_const_value_type            xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type> IPT;
  typedef Kokkos::Details::ArithTraits<typename IPT::mag_type>   AT;
  typedef typename IPT::mag_type                       value_type[];

  size_type value_count;
  RV m_r;
  typename XMV::const_type m_x;

  MV_Nrm2Squared_Right_FunctorVector (const RV& r, const XMV& x) :
    value_count (x.dimension_1 ()), m_r (r), m_x (x)
  {
#ifdef KOKKOS_HAVE_CXX11
    static_assert (Kokkos::Impl::is_view<RV>::value,
                   "KokkosBlas::Impl::MV_Nrm2Squared_Right_FunctorVector: "
                   "R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value,
                   "KokkosBlas::Impl::MV_Nrm2Squared_Right_FunctorVector: "
                   "X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Nrm2Squared_Right_FunctorVector: "
                   "R is const.  It must be nonconst, because it is an output "
                   "argument (we must be able to write to its entries).");
    static_assert (RV::rank == 1 && XMV::rank == 2,
                   "KokkosBlas::Impl::MV_Nrm2Squared_Right_FunctorVector: "
                   "RV must have rank 1 and XMV must have rank 2.");
#endif // KOKKOS_HAVE_CXX11
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type i, value_type sum) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      const typename IPT::mag_type tmp = IPT::norm (m_x(i,j));
      sum[j] += tmp * tmp;
    }
  }

  KOKKOS_INLINE_FUNCTION void
  init (value_type update) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      update[j] = AT::zero ();
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type update,
        const volatile value_type source) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      update[j] += source[j];
    }
  }

  // On device, write the reduction result to the output View.
  KOKKOS_INLINE_FUNCTION void
  final (const value_type dst) const
  {
    const size_type numVecs = value_count;
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type j = 0; j < numVecs; ++j) {
      m_r(j) = dst[j];
    }
  }
};


/// \brief Compute the square of the 2-norm of the single vector (1-D
///   View) X, and store the result in the 0-D View r.
template<class RV, class XV, class SizeType>
void
V_Nrm2Squared_Invoke (const RV& r, const XV& X)
{
  typedef typename XV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  typedef V_Nrm2Squared_Functor<RV, XV, SizeType> functor_type;
  functor_type op (r, X);
  Kokkos::parallel_reduce (policy, op);
}


/// \brief Compute the squares of the 2-norms of the columns of the
///   multivector (2-D View) X, and store result(s) in the 1-D View r.
template<class RV, class XMV, class SizeType>
void
MV_Nrm2Squared_Invoke (const RV& r, const XMV& X)
{
  typedef typename XMV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

#ifdef KOKKOS_HAVE_CXX11
  // We only give you a single-vector special case if you build with
  // C++11 enabled, since we need 'decltype' to ensure that we have
  // the right layouts for RV0D and XV1D.
  if (X.dimension_1 () == 1) {
    auto r_0 = Kokkos::subview (r, 0);
    auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
    typedef decltype (r_0) RV0D;
    typedef decltype (X_0) XV1D;
    V_Nrm2Squared_Invoke<RV0D, XV1D, SizeType> (r_0, X_0);
    return;
  }
#endif // KOKKOS_HAVE_CXX11

  typedef MV_Nrm2Squared_Right_FunctorVector<RV, XMV, SizeType> functor_type;
  functor_type op (r, X);
  Kokkos::parallel_reduce (policy, op);
}


/// \brief Implementation of KokkosBlas::nrm2_squared for multivectors
///   and single vectors.
template<class RV, class XMV, int rank = XMV::rank>
struct Nrm2_MV {};


//! Special case for multivectors (rank-2 Views).
template<class RV, class XMV>
struct Nrm2_MV<RV, XMV, 2> {
  /// \brief Compute the squares of the 2-norm(s) of the column(s) of
  ///   the multivector (2-D View) X, and store result(s) in r.
  static void nrm2_squared (const RV& r, const XMV& X)
  {
    typedef typename XMV::size_type size_type;
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      MV_Nrm2Squared_Invoke<RV, XMV, int> (r, X);
    }
    else {
      MV_Nrm2Squared_Invoke<RV, XMV, size_type> (r, X);
    }
  }
};


//! Special case for single vectors (rank-1 Views).
template<class RV, class XV>
struct Nrm2_MV<RV, XV, 1> {
  /// \brief Compute the square of the 2-norm of the vector (1-D View)
  ///   X, and store the result in the 0-D View r.
  static void nrm2_squared (const RV& r, const XV& X)
  {
    typedef typename XV::size_type size_type;
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      V_Nrm2Squared_Invoke<RV, XV, int> (r, X);
    }
    else {
      V_Nrm2Squared_Invoke<RV, XV, size_type> (r, X);
    }
  }
};


// Full specializations for cases of interest for Tpetra::MultiVector.
//
// Currently, we include specializations for Scalar = double,
// LayoutLeft (which is what Tpetra::MultiVector uses at the moment),
// and all execution spaces.  This may change in the future.  The
// output View _always_ uses the execution space's default array
// layout, which is what Tpetra::MultiVector wants for the output
// argument of norm1().

#ifdef KOKKOS_HAVE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Nrm2_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                            KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                            Kokkos::LayoutLeft,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  static void nrm2_squared (const RV& r, const XMV& X);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Nrm2_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                            KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                            Kokkos::LayoutLeft,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  static void nrm2_squared (const RV& r, const XMV& X);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Nrm2_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                            KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                            Kokkos::LayoutLeft,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  static void nrm2_squared (const RV& r, const XMV& X);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Nrm2_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                            KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                            Kokkos::LayoutLeft,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  static void nrm2_squared (const RV& r, const XMV& X);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double
template<>
struct Nrm2_MV<Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                            KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                            Kokkos::LayoutLeft,
                            Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                            Kokkos::Impl::ViewDefault>,
               2>
{
  typedef Kokkos::View<Kokkos::Details::InnerProductSpaceTraits<KOKKOSBLAS_IMPL_MV_SCALAR>::mag_type*,
                       KOKKOSBLAS_IMPL_MV_EXEC_SPACE::array_layout,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> RV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> XMV;
  static void nrm2_squared (const RV& r, const XMV& X);
};
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_SCALAR
#endif // KOKKOS_HAVE_CUDA


} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_NRM2_HPP_
