/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOS_BLAS1_MV_IMPL_FILL_HPP_
#define KOKKOS_BLAS1_MV_IMPL_FILL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <type_traits>
#include <cstring> // for memset (see KokkosBlas::fill specializations below)

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
///
/// \tparam XMV Type of the Kokkos::View to fill.
/// \tparam rank (Integer) rank of XMV.  1 means a 1-D View, 2 means a
///   2-D View, etc.
/// \tparam useMemset Whether to use memset if possible (if the View
///   occupies a contiguous chunk of memory).  This depends on the
///   value type in the View, as well as the View's layout, execution
///   space (we only use memset for Kokkos::Serial, else no
///   parallelism), and memory space (it must be host-reachable
///   memory).
///
template<class XMV,
         const int rank = XMV::rank,
         const bool useMemset =
#ifdef KOKKOS_HAVE_SERIAL
           std::is_fundamental<typename XMV::non_const_value_type>::value &&
           (std::is_integral<typename XMV::non_const_value_type>::value ||
            std::is_floating_point<typename XMV::non_const_value_type>::value) &&
           (std::is_same<typename XMV::array_layout, Kokkos::LayoutLeft>::value ||
            std::is_same<typename XMV::array_layout, Kokkos::LayoutRight>::value) &&
           std::is_same<typename XMV::memory_space, Kokkos::HostSpace>::value &&
           std::is_same<typename XMV::execution_space, Kokkos::Serial>::value
#else // NOT KOKKOS_HAVE_SERIAL
         false
#endif // KOKKOS_HAVE_SERIAL
         >
struct Fill;

// Specialization for multivectors (2-D Views), not using memset at all.
template<class XMV>
struct Fill<XMV, 2, false>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
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
}
#endif
;

// Specialization for multivectors (2-D Views), using memset if possible.
template<class XMV>
struct Fill<XMV, 2, true>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  static void fill (const XMV& X, const typename XMV::non_const_value_type& val)
  {
    typedef typename XMV::size_type size_type;
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // Don't call one of the special cases (memset or 1-D fill) unless
    // the memory in X is contiguous.
    if (X.capacity () == numRows * numCols) {
      typedef typename XMV::non_const_value_type SC;

      if (val == Kokkos::Details::ArithTraits<SC>::zero ()) {
        // It might not necessarily be true for ALL Scalar types that
        // memset with 0 does the right thing, but it certainly works
        // here.
        memset (X.ptr_on_device (), 0, numRows * numCols * sizeof (SC));
      }
      else {
        typedef Kokkos::View<SC*,
          typename XMV::array_layout,
          typename XMV::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV1D;

        XV1D X1D (X.ptr_on_device (), X.capacity ());
        Kokkos::deep_copy (X1D, val);

        // mfh 14 Apr 2015: The code below didn't actually help
        // performance over using ViewFill (mfh 21 Mar 2017: or
        // Kokkos::deep_copy?) on the 1-D View.  The key thing is
        // using the right ViewFill specialization: 1-D is faster (.18
        // s / 26 s, for example).

        // if (numRows < static_cast<size_type> (INT_MAX)) {
        //   V_Fill_Invoke<XV1D, int> (X1D, val);
        // }
        // else {
        //   V_Fill_Invoke<XV1D, size_type> (X1D, val);
        // }
      }
    }
    else {
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
  }
}
#endif
;

// Specialization for single vectors (1-D Views), not using memset at all.
template<class XV>
struct Fill<XV, 1, false>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
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
}
#endif
;

// Specialization for single vectors (1-D Views), using memset if possible.
template<class XV>
struct Fill<XV, 1, true>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  static void fill (const XV& X, const typename XV::non_const_value_type& val)
  {
    typedef typename XV::non_const_value_type SC;
    typedef typename XV::size_type size_type;
    const size_type numRows = X.dimension_0 ();

    if (val == Kokkos::Details::ArithTraits<SC>::zero ()) {
      memset (X.ptr_on_device (), 0, numRows * sizeof (SC));
    }
    else {
      Kokkos::deep_copy (X, val);
    }
  }
}
#endif
;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Fill for rank == 2.  This is NOT for users!!!
// All the declarations of full specializations go in this header
// file.  We may spread out definitions (see _DEF macro below) across
// one or more .cpp files.
//

#define KOKKOSBLAS1_IMPL_MV_FILL_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Fill<Kokkos::View<SCALAR**, \
                              LAYOUT, \
                              Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                            2>;

//
// KokkosBlas::Impl::Fill for rank == 2.  This is NOT for users!!!  We
// may spread out use of this macro across one or more .cpp files in
// this directory.
//

#define KOKKOSBLAS1_IMPL_MV_FILL_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Fill<Kokkos::View<SCALAR**, \
                       LAYOUT, \
                       Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                     2>;

} // namespace Impl
} // namespace KokkosBlas


#include<generated_specializations_hpp/KokkosBlas1_impl_MV_fill_decl_specializations.hpp>
#endif // KOKKOS_BLAS1_MV_IMPL_FILL_HPP_
