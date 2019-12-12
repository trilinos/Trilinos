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

/// \file KokkosSparse_sptrsv.hpp
/// \brief Parallel sparse triangular solve
///
/// This file provides KokkosSparse::sptrsv.  This function performs a
/// local (no MPI) sparse triangular solve on matrices stored in
/// compressed row sparse ("Crs") format.

#ifndef KOKKOSSPARSE_SPTRSV_HPP_
#define KOKKOSSPARSE_SPTRSV_HPP_

#include <type_traits>

//#include "KokkosSparse_sptrsv_handle.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosSparse_sptrsv_symbolic_spec.hpp"
#include "KokkosSparse_sptrsv_solve_spec.hpp"

namespace KokkosSparse {
namespace Experimental {

#define KOKKOSKERNELS_SPTRSV_SAME_TYPE(A, B) std::is_same<typename std::remove_const<A>::type, typename std::remove_const<B>::type>::value

  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_>
  void sptrsv_symbolic(
      KernelHandle *handle, 
      lno_row_view_t_ rowmap,
      lno_nnz_view_t_ entries)
  {
    typedef typename KernelHandle::size_type size_type;
    typedef typename KernelHandle::nnz_lno_t ordinal_type;

    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
        "sptrsv_symbolic: A size_type must match KernelHandle size_type (const doesn't matter)");

    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
        "sptrsv_symbolic: A entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");


    typedef typename KernelHandle::const_size_type c_size_t;
    typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
    typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

    typedef typename KernelHandle::HandleExecSpace c_exec_t;
    typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
    typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

    typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
    const_handle_type tmp_handle (*handle);

    typedef Kokkos::View<
          typename lno_row_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
          typename lno_row_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > RowMap_Internal;

    typedef Kokkos::View<
          typename lno_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
          typename lno_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Entries_Internal;


    RowMap_Internal rowmap_i = rowmap;
    Entries_Internal entries_i = entries;

    KokkosSparse::Impl::SPTRSV_SYMBOLIC<const_handle_type, RowMap_Internal, Entries_Internal>::sptrsv_symbolic (&tmp_handle, rowmap_i, entries_i);

  } // sptrsv_symbolic



  template <typename KernelHandle,
            typename lno_row_view_t_,
            typename lno_nnz_view_t_,
            typename scalar_nnz_view_t_,
            class BType,
            class XType>
  void sptrsv_solve(
      KernelHandle *handle, 
      lno_row_view_t_ rowmap,
      lno_nnz_view_t_ entries,
      scalar_nnz_view_t_ values,
      BType b,
      XType x)
  {
    typedef typename KernelHandle::size_type size_type;
    typedef typename KernelHandle::nnz_lno_t ordinal_type;
    typedef typename KernelHandle::nnz_scalar_t scalar_type;
    
    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_row_view_t_::non_const_value_type, size_type),
        "sptrsv_solve: A size_type must match KernelHandle size_type (const doesn't matter)");
    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename lno_nnz_view_t_::non_const_value_type, ordinal_type),
        "sptrsv_solve: A entry type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");
    static_assert(KOKKOSKERNELS_SPTRSV_SAME_TYPE(typename scalar_nnz_view_t_::value_type, scalar_type),
        "sptrsv_solve: A scalar type must match KernelHandle entry type (aka nnz_lno_t, and const doesn't matter)");

    static_assert (Kokkos::Impl::is_view<BType>::value,
        "sptrsv: b is not a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<XType>::value,
        "sptrsv: x is not a Kokkos::View.");
    static_assert ((int) BType::rank == (int) XType::rank,
        "sptrsv: The ranks of b and x do not match.");
    static_assert (BType::rank == 1,
        "sptrsv: b and x must both either have rank 1.");
    static_assert (std::is_same<typename XType::value_type,
        typename XType::non_const_value_type>::value,
        "sptrsv: The output x must be nonconst.");
    static_assert (std::is_same<typename BType::device_type, typename XType::device_type>::value,
        "sptrsv: Views BType and XType have different device_types.");
    static_assert (std::is_same<typename BType::device_type::execution_space, typename KernelHandle::SPTRSVHandleType::execution_space>::value,
        "sptrsv: KernelHandle and Views have different execution spaces.");
    static_assert (std::is_same<typename lno_row_view_t_::device_type, typename lno_nnz_view_t_::device_type>::value,
        "sptrsv: rowmap and entries have different device types.");
    static_assert (std::is_same<typename lno_row_view_t_::device_type, typename scalar_nnz_view_t_::device_type>::value,
        "sptrsv: rowmap and values have different device types.");


    typedef typename KernelHandle::const_size_type c_size_t;
    typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
    typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

    typedef typename KernelHandle::HandleExecSpace c_exec_t;
    typedef typename KernelHandle::HandleTempMemorySpace c_temp_t;
    typedef typename KernelHandle::HandlePersistentMemorySpace c_persist_t;

    typedef typename  KokkosKernels::Experimental::KokkosKernelsHandle<c_size_t, c_lno_t, c_scalar_t, c_exec_t, c_temp_t, c_persist_t> const_handle_type;
    const_handle_type tmp_handle (*handle);

    typedef Kokkos::View<
          typename lno_row_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_row_view_t_>::array_layout,
          typename lno_row_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > RowMap_Internal;

    typedef Kokkos::View<
          typename lno_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<lno_nnz_view_t_>::array_layout,
          typename lno_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Entries_Internal;

    typedef Kokkos::View<
          typename scalar_nnz_view_t_::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<scalar_nnz_view_t_>::array_layout,
          typename scalar_nnz_view_t_::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > Values_Internal;


    typedef Kokkos::View<
          typename BType::const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<BType>::array_layout,
          typename BType::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > BType_Internal;

    typedef Kokkos::View<
          typename XType::non_const_value_type*,
          typename KokkosKernels::Impl::GetUnifiedLayout<XType>::array_layout,
          typename XType::device_type,
          Kokkos::MemoryTraits<Kokkos::Unmanaged> > XType_Internal;


    RowMap_Internal rowmap_i = rowmap;
    Entries_Internal entries_i = entries;
    Values_Internal values_i = values;

    BType_Internal b_i = b;
    XType_Internal x_i = x;

    KokkosSparse::Impl::SPTRSV_SOLVE<const_handle_type, RowMap_Internal, Entries_Internal, Values_Internal, BType_Internal, XType_Internal>::sptrsv_solve (&tmp_handle, rowmap_i, entries_i, values_i, b_i, x_i);

  } // sptrsv_solve

} // namespace Experimental
} // namespace KokkosSparse

#undef KOKKOSKERNELS_SPTRSV_SAME_TYPE

#endif // KOKKOSSPARSE_SPTRSV_HPP_

