/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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
#ifndef _KOKKOSGRAPH_DISTANCE1_COLOR_HPP
#define _KOKKOSGRAPH_DISTANCE1_COLOR_HPP

#include "KokkosGraph_color_d1_spec.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Utils.hpp"

namespace KokkosGraph {

namespace Experimental {

template <class KernelHandle, typename lno_row_view_t_,
          typename lno_nnz_view_t_>
void graph_color_symbolic(KernelHandle *handle,
                          typename KernelHandle::nnz_lno_t num_rows,
                          typename KernelHandle::nnz_lno_t /* num_cols */,
                          lno_row_view_t_ row_map, lno_nnz_view_t_ entries,
                          bool /* is_symmetric */ = true) {
  typedef typename KernelHandle::HandleExecSpace ExecSpace;
  typedef typename KernelHandle::HandleTempMemorySpace MemSpace;
  typedef typename KernelHandle::HandlePersistentMemorySpace PersistentMemSpace;
  typedef typename Kokkos::Device<ExecSpace, MemSpace> DeviceType;

  typedef typename KernelHandle::const_size_type c_size_t;
  typedef typename KernelHandle::const_nnz_lno_t c_lno_t;
  typedef typename KernelHandle::const_nnz_scalar_t c_scalar_t;

  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<
      c_size_t, c_lno_t, c_scalar_t, ExecSpace, MemSpace, PersistentMemSpace>
      ConstKernelHandle;
  ConstKernelHandle tmp_handle(*handle);

  typedef Kokkos::View<typename lno_row_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           lno_row_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_rowmap;
  typedef Kokkos::View<typename lno_nnz_view_t_::const_value_type *,
                       typename KokkosKernels::Impl::GetUnifiedLayout<
                           lno_nnz_view_t_>::array_layout,
                       DeviceType, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      Internal_entries;
  KokkosGraph::Impl::
      COLOR_D1<ConstKernelHandle, Internal_rowmap, Internal_entries>::color_d1(
          &tmp_handle, num_rows,
          Internal_rowmap(row_map.data(), row_map.extent(0)),
          Internal_entries(entries.data(), entries.extent(0)));
}

template <class KernelHandle, typename lno_row_view_t_,
          typename lno_nnz_view_t_>
void graph_color(KernelHandle *handle,
                 typename KernelHandle::nnz_lno_t num_rows,
                 typename KernelHandle::nnz_lno_t num_cols,
                 lno_row_view_t_ row_map, lno_nnz_view_t_ entries,
                 bool is_symmetric = true) {
  graph_color_symbolic(handle, num_rows, num_cols, row_map, entries,
                       is_symmetric);
}

}  // end namespace Experimental
}  // end namespace KokkosGraph

#endif  // _KOKKOSGRAPH_DISTANCE1_COLOR_HPP
