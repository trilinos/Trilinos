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
#ifndef KOKKOS_BLAS1_MV_IMPL_MULT_HPP_
#define KOKKOS_BLAS1_MV_IMPL_MULT_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Functor for entry-wise multiply of multivectors.
///
/// \tparam CMV 2-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BMV 2-D Kokkos::View
/// \tparam scalar_ab 0 if ab is zero, else nonzero (preferably 2).
/// \tparam scalar_c 0 if c is zero, else nonzero (preferably 2).
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i,j) = c * C(i,j) + ab * A(i) * B(i,j), subject to the usual
/// BLAS update rules.
template<class CMV, class AV, class BMV,
         int scalar_ab, int scalar_c,
         class SizeType = typename CMV::size_type>
struct MV_MultFunctor
{
  typedef typename CMV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename CMV::non_const_value_type> ATS;

  const size_type m_n;
  typename CMV::const_value_type m_c;
  CMV m_C;
  typename AV::const_value_type m_ab;
  AV m_A;
  BMV m_B;

  MV_MultFunctor (typename CMV::const_value_type& c,
                  const CMV& C,
                  typename AV::const_value_type& ab,
                  const AV& A,
                  const BMV& B) :
    m_n (C.dimension_1 ()),
    m_c (c), m_C (C), m_ab (ab), m_A (A), m_B (B)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i) const
  {
    if (scalar_c == 0) {
      if (scalar_ab == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i,j) = ATS::zero ();
        }
      }
      else { // ab != 0, c == 0
        typename AV::const_value_type Ai = m_A(i);
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i,j) = m_ab * Ai * m_B(i,j);
        }
      }
    } else { // c != 0
      if (scalar_ab == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i,j) = m_c * m_C(i,j);
        }
      }
      else { // m_ab != 0, and m_c != 0
        typename AV::const_value_type Ai = m_A(i);
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
        for (size_type j = 0; j < m_n; ++j) {
          m_C(i,j) = m_c * m_C(i,j) + m_ab * Ai * m_B(i,j);
        }
      }
    }
  }
};

/// \brief Functor for entry-wise multiply of vectors.
///
/// \tparam CV 1-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BV 1-D Kokkos::View
/// \tparam scalar_ab 0 if ab is zero, else nonzero (preferably 2).
/// \tparam scalar_c 0 if c is zero, else nonzero (preferably 2).
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i) = c * C(i) + ab * A(i) * B(i), subject to the usual
/// BLAS update rules.
template<class CV, class AV, class BV,
         int scalar_ab, int scalar_c,
         class SizeType = typename CV::size_type>
struct V_MultFunctor
{
  typedef typename CV::execution_space execution_space;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<typename CV::non_const_value_type> ATS;

  typename CV::const_value_type m_c;
  CV m_C;
  typename AV::const_value_type m_ab;
  AV m_A;
  BV m_B;

  V_MultFunctor (typename CV::const_value_type& c,
                 const CV& C,
                 typename AV::const_value_type& ab,
                 const AV& A,
                 const BV& B) :
    m_c (c), m_C (C), m_ab (ab), m_A (A), m_B (B)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const size_type& i) const
  {
    if (scalar_c == 0) {
      if (scalar_ab == 0) {
        m_C(i) = ATS::zero ();
      }
      else { // ab != 0, c == 0
        m_C(i) = m_ab * m_A(i) * m_B(i);
      }
    } else { // c != 0
      if (scalar_ab == 0) {
        m_C(i) = m_c * m_C(i);
      }
      else { // m_ab != 0, and m_c != 0
        m_C(i) = m_c * m_C(i) + m_ab * m_A(i) * m_B(i);
      }
    }
  }
};

/// \brief Implementation of entry-wise multiply of vectors, that
///   dispatches to the right functor invocation.
///
/// \tparam CV 1-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BV 1-D Kokkos::View
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i) = c * C(i) + ab * A(i) * B(i), subject to the usual BLAS
/// update rules.
template<class CV, class AV, class BV, class SizeType>
void
V_Mult_Generic (typename CV::const_value_type& c,
                const CV& C,
                typename AV::const_value_type& ab,
                const AV& A,
                const BV& B)
{
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::Details::ArithTraits<typename AV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename CV::non_const_value_type> ATC;
  typedef typename CV::execution_space execution_space;

  const SizeType numRows = C.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (c == ATC::zero ()) {
    if (ab == ATA::zero ()) {
      typedef V_MultFunctor<CV, AV, BV, 0, 0, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
    else {
      typedef V_MultFunctor<CV, AV, BV, 2, 0, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
  }
  else { // c != 0
    if (ab == ATA::zero ()) {
      typedef V_MultFunctor<CV, AV, BV, 0, 2, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
    else {
      typedef V_MultFunctor<CV, AV, BV, 2, 2, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
  }
}

/// \brief Implementation of entry-wise multiply of multivectors, that
///   dispatches to the right functor invocation (or calls
///   V_Mult_Generic if C and B have one column).
///
/// \tparam CMV 2-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BMV 2-D Kokkos::View
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i,j) = c * C(i,j) + ab * A(i) * B(i,j), subject to the usual
/// BLAS update rules.
template<class CMV, class AV, class BMV, class SizeType>
void
MV_Mult_Generic (typename CMV::const_value_type& c,
                 const CMV& C,
                 typename AV::const_value_type& ab,
                 const AV& A,
                 const BMV& B)
{
  using Kokkos::ALL;
  using Kokkos::subview;
  typedef Kokkos::Details::ArithTraits<typename AV::non_const_value_type> ATA;
  typedef Kokkos::Details::ArithTraits<typename CMV::non_const_value_type> ATC;
  typedef typename CMV::execution_space execution_space;

  if (C.dimension_1 () == 1) {
    // It's better to use decltype if we have C++11, because that will
    // always work, no matter the layout of CMV and BMV.  If either
    // has LayoutRight, then the correct layout of a single column is
    // LayoutStride, which the non-C++11 branch below won't get right.
#ifdef KOKKOS_HAVE_CXX11
    auto C_0 = subview (C, ALL (), 0);
    auto B_0 = subview (B, ALL (), 0);
    typedef decltype (C_0) CV;
    typedef decltype (B_0) BV;
#else
    typedef Kokkos::View<typename CMV::value_type*,
      typename CMV::array_layout, typename CMV::device_type,
      typename CMV::memory_traits, typename CMV::specialize> CV;
    typedef Kokkos::View<typename BMV::value_type*,
      typename BMV::array_layout, typename BMV::device_type,
      typename BMV::memory_traits, typename BMV::specialize> BV;
    CV C_0 = subview (C, ALL (), 0);
    BV B_0 = subview (B, ALL (), 0);
#endif // KOKKOS_HAVE_CXX11
    V_Mult_Generic<CV, AV, BV, SizeType> (c, C_0, ab, A, B_0);
    return;
  }

  const SizeType numRows = C.dimension_0 ();
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  if (c == ATC::zero ()) {
    if (ab == ATA::zero ()) {
      typedef MV_MultFunctor<CMV, AV, BMV, 0, 0, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
    else {
      typedef MV_MultFunctor<CMV, AV, BMV, 2, 0, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
  }
  else { // c != 0
    if (ab == ATA::zero ()) {
      typedef MV_MultFunctor<CMV, AV, BMV, 0, 2, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
    else {
      typedef MV_MultFunctor<CMV, AV, BMV, 2, 2, SizeType> functor_type;
      functor_type op (c, C, ab, A, B);
      Kokkos::parallel_for (policy, op);
    }
  }
}

/// \brief Implementation of entry-wise multiply of multivectors or
///   single vectors (depending on the rank template parameter).
///
/// This struct is how the front-end interface KokkosBlas::mult talks
/// to the implementation.  The specializations for rank == 1 and rank
/// == 2 have meaningful content.
template<class CMV, class AV, class BMV, int rank = CMV::rank>
struct Mult {};

/// \brief Implementation of entry-wise multiply of multivectors, that
///   dispatches to MV_Mult_Generic.
///
/// \tparam CMV 2-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BMV 2-D Kokkos::View
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i,j) = c * C(i,j) + ab * A(i) * B(i,j), subject to the usual
/// BLAS update rules.
template<class CMV, class AV, class BMV>
struct Mult<CMV, AV, BMV, 2> {
  static void
  mult (typename CMV::const_value_type& c,
        const CMV& C,
        typename AV::const_value_type& ab,
        const AV& A,
        const BMV& B)
  {
    typedef typename CMV::size_type size_type;

    const size_type numRows = C.dimension_0 ();
    const size_type numCols = C.dimension_1 ();
    if (numRows < static_cast<int> (INT_MAX) &&
        numRows * numCols < static_cast<int> (INT_MAX)) {
      MV_Mult_Generic<CMV, AV, BMV, int> (c, C, ab, A, B);
    }
    else {
      MV_Mult_Generic<CMV, AV, BMV, size_type> (c, C, ab, A, B);
    }
  }
};

/// \brief Implementation of entry-wise multiply of vectors, that
///   dispatches to V_Mult_Generic.
///
/// \tparam CV 1-D Kokkos::View
/// \tparam AV 1-D Kokkos::View
/// \tparam BV 1-D Kokkos::View
/// \tparam SizeType Index type for iterating over rows.
///
/// C(i) = c * C(i) + ab * A(i) * B(i), subject to the usual
/// BLAS update rules.
template<class CV, class AV, class BV>
struct Mult<CV, AV, BV, 1> {
  static void
  mult (typename CV::const_value_type& c,
        const CV& C,
        typename AV::const_value_type& ab,
        const AV& A,
        const BV& B)
  {
    typedef typename CV::size_type size_type;

    const size_type numRows = C.dimension_0 ();
    if (numRows < static_cast<int> (INT_MAX)) {
      V_Mult_Generic<CV, AV, BV, int> (c, C, ab, A, B);
    }
    else {
      V_Mult_Generic<CV, AV, BV, size_type> (c, C, ab, A, B);
    }
  }
};

// Full specializations for cases of interest for Tpetra::MultiVector.
//
// Currently, we include specializations for Scalar = double,
// LayoutLeft (which is what Tpetra::MultiVector uses at the moment),
// and all execution spaces.  This may change in the future.

#ifdef KOKKOS_HAVE_SERIAL
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Serial
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Mult<Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                         Kokkos::LayoutLeft,
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
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> CMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> BMV;

  static void
  mult (CMV::const_value_type& c,
        const CMV& C,
        AV::const_value_type& ab,
        const AV& A,
        const BMV& B);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::OpenMP
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Mult<Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                         Kokkos::LayoutLeft,
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
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> CMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> BMV;

  static void
  mult (CMV::const_value_type& c,
        const CMV& C,
        AV::const_value_type& ab,
        const AV& A,
        const BMV& B);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Threads
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::HostSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Mult<Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                         Kokkos::LayoutLeft,
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
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> CMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> BMV;

  static void
  mult (CMV::const_value_type& c,
        const CMV& C,
        AV::const_value_type& ab,
        const AV& A,
        const BMV& B);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Mult<Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                         Kokkos::LayoutLeft,
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
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> CMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> BMV;

  static void
  mult (CMV::const_value_type& c,
        const CMV& C,
        AV::const_value_type& ab,
        const AV& A,
        const BMV& B);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA
#define KOKKOSBLAS_IMPL_MV_EXEC_SPACE Kokkos::Cuda
#define KOKKOSBLAS_IMPL_MV_MEM_SPACE Kokkos::CudaUVMSpace
#define KOKKOSBLAS_IMPL_MV_SCALAR double

template<>
struct Mult<Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                         Kokkos::LayoutLeft,
                         Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                         Kokkos::Impl::ViewDefault>,
            Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                         Kokkos::LayoutLeft,
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
  typedef Kokkos::View<KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> CMV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR*,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> AV;
  typedef Kokkos::View<const KOKKOSBLAS_IMPL_MV_SCALAR**,
                       Kokkos::LayoutLeft,
                       Kokkos::Device<KOKKOSBLAS_IMPL_MV_EXEC_SPACE, KOKKOSBLAS_IMPL_MV_MEM_SPACE>,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                       Kokkos::Impl::ViewDefault> BMV;

  static void
  mult (CMV::const_value_type& c,
        const CMV& C,
        AV::const_value_type& ab,
        const AV& A,
        const BMV& B);
};

#undef KOKKOSBLAS_IMPL_MV_SCALAR
#undef KOKKOSBLAS_IMPL_MV_MEM_SPACE
#undef KOKKOSBLAS_IMPL_MV_EXEC_SPACE
#endif // KOKKOS_HAVE_CUDA

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_MV_IMPL_MULT_HPP_
