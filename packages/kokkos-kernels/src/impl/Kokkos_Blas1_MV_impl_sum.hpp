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
#ifndef KOKKOS_BLAS1_MV_IMPL_SUM_HPP_
#define KOKKOS_BLAS1_MV_IMPL_SUM_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace KokkosBlas {
namespace Impl {

//
// sum
//

/// \brief Sum functor for single vectors.
///
/// \tparam RV 0-D output View
/// \tparam XV 1-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XV, class SizeType = typename XV::size_type>
struct V_Sum_Functor
{
  typedef typename XV::execution_space   execution_space;
  typedef SizeType                             size_type;
  typedef typename XV::non_const_value_type   value_type;
  typedef Kokkos::Details::ArithTraits<value_type>    AT;

  RV m_r;
  XV m_x;

  V_Sum_Functor (const RV& r, const XV& x) :
    m_r (r), m_x (x)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value,
                   "KokkosBlas::Impl::V_Sum_Functor: R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XV>::value,
                   "KokkosBlas::Impl::V_Sum_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                     typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::V_Sum_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (RV::rank == 0 && XV::rank == 1,
                   "KokkosBlas::Impl::V_Sum_Functor: "
                   "RV must have rank 0 and XV must have rank 1.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i, value_type& sum) const
  {
    sum += m_x(i);
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
  KOKKOS_INLINE_FUNCTION void final (const value_type& dst) const
  {
    m_r() = dst;
  }
};

/// \brief Sum functor for multivectors.
///
/// \tparam RV 1-D output View
/// \tparam XMV 2-D input View
/// \tparam SizeType Index type.  Use int (32 bits) if possible.
template<class RV, class XMV, class SizeType = typename XMV::size_type>
struct MV_Sum_Functor {
  typedef typename XMV::execution_space                       execution_space;
  typedef SizeType                                                  size_type;
  typedef typename XMV::non_const_value_type                     value_type[];
  typedef Kokkos::Details::ArithTraits<typename XMV::non_const_value_type> AT;

  const size_type value_count;
  RV sums_;
  XMV X_;

  MV_Sum_Functor (const RV& sums, const XMV& X) :
    value_count (X.dimension_1 ()), sums_ (sums), X_ (X)
  {
    static_assert (Kokkos::Impl::is_view<RV>::value,
                   "KokkosBlas::Impl::MV_Sum_Functor: R is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XMV>::value,
                   "KokkosBlas::Impl::MV_Sum_Functor: X is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_same<typename RV::value_type,
                   typename RV::non_const_value_type>::value,
                   "KokkosBlas::Impl::MV_Sum_Functor: R is const.  "
                   "It must be nonconst, because it is an output argument "
                   "(we have to be able to write to its entries).");
    static_assert (RV::rank == 1 && XMV::rank == 2,
                   "KokkosBlas::Impl::MV_Sum_Functor: "
                   "RV must have rank 1 and XMV must have rank 2.");
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const size_type& i, value_type sum) const
  {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_VECTOR
#pragma vector always
#endif
    for (size_type k = 0; k < value_count; ++k) {
      sum[k] += X_(i,k);
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
      update[k] += source[k];
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
      sums_(k) = dst[k];
    }
  }
};

/// \brief Compute the sum of the entries of the single vector (1-D
///   View) X, and store the result in the 0-D View r.
template<class RV, class XV, class SizeType>
void
V_Sum_Invoke (const RV& r, const XV& X)
{
  typedef typename XV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  typedef V_Sum_Functor<RV, XV, SizeType> functor_type;
  functor_type op (r, X);
  Kokkos::parallel_reduce (policy, op);
}


/// \brief Compute the sums of the entries of the columns of the
///   multivector (2-D View) X, and store result(s) in the 1-D View r.
template<class RV, class XMV, class SizeType>
void
MV_Sum_Invoke (const RV& r, const XMV& X)
{
  typedef typename XMV::execution_space execution_space;
  const SizeType numRows = static_cast<SizeType> (X.dimension_0 ());
  Kokkos::RangePolicy<execution_space, SizeType> policy (0, numRows);

  // For the single-vector special case, we use 'decltype' to ensure
  // that we have the right layouts for RV0D and XV1D.
  if (X.dimension_1 () == 1) {
    auto r_0 = Kokkos::subview (r, 0);
    auto X_0 = Kokkos::subview (X, Kokkos::ALL (), 0);
    typedef decltype (r_0) RV0D;
    typedef decltype (X_0) XV1D;
    V_Sum_Invoke<RV0D, XV1D, SizeType> (r_0, X_0);
    return;
  }

  typedef MV_Sum_Functor<RV, XMV, SizeType> functor_type;
  functor_type op (r, X);
  Kokkos::parallel_reduce (policy, op);
}


/// \brief Implementation of KokkosBlas::sum for multivectors and
///   single vectors.
template<class RV, class XMV, int rank = XMV::rank>
struct Sum;


//! Special case for multivectors (rank-2 Views).
template<class RV, class XMV>
struct Sum<RV, XMV, 2>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  /// \brief Compute the sums of the entries of the column(s) of the
  ///   multivector (2-D View) X, and store result(s) in r.
  static void sum (const RV& r, const XMV& X)
  {
    typedef typename XMV::size_type size_type;
    const size_type numRows = X.dimension_0 ();
    const size_type numCols = X.dimension_1 ();

    // int is generally faster than size_t, but check for overflow first.
    if (numRows < static_cast<size_type> (INT_MAX) &&
        numRows * numCols < static_cast<size_type> (INT_MAX)) {
      MV_Sum_Invoke<RV, XMV, int> (r, X);
    }
    else {
      MV_Sum_Invoke<RV, XMV, size_type> (r, X);
    }
  }
}
#endif
;


//! Special case for single vectors (rank-1 Views).
template<class RV, class XV>
struct Sum<RV, XV, 1>
#ifndef KOKKOSKERNELS_ETI_ONLY
{
  /// \brief Compute the sum of the entries of the vector (1-D View)
  ///   X, and store the result in the 0-D View r.
  static void sum (const RV& r, const XV& X)
  {
    typedef typename XV::size_type size_type;
    // int is generally faster than size_t, but check for overflow first.
    if (X.dimension_0 () < static_cast<size_type> (INT_MAX)) {
      V_Sum_Invoke<RV, XV, int> (r, X);
    }
    else {
      V_Sum_Invoke<RV, XV, size_type> (r, X);
    }
  }
}
#endif
;

//
// Macro for declaration of full specialization of
// KokkosBlas::Impl::Sum for rank == 2.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//

#define KOKKOSBLAS1_IMPL_MV_SUM_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
extern template struct Sum<Kokkos::View<SCALAR*, \
                        EXEC_SPACE::array_layout, \
                        Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                        Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
           Kokkos::View<const SCALAR**, \
                        LAYOUT, \
                        Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                        Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                           2>;

//
// Macro for definition of full specialization of
// KokkosBlas::Impl::Sum for rank == 2.  This is NOT for users!!!
//

#define KOKKOSBLAS1_IMPL_MV_SUM_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template struct Sum<Kokkos::View<SCALAR*, \
                 EXEC_SPACE::array_layout, \
                 Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
    Kokkos::View<const SCALAR**, \
                 LAYOUT, \
                 Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
                    2>;

} // namespace Impl
} // namespace KokkosBlas

#include<generated_specializations_hpp/KokkosBlas1_impl_MV_sum_decl_specializations.hpp>
#endif // KOKKOS_BLAS1_MV_IMPL_SUM_HPP_
