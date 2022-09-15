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
#ifndef KOKKOSSPARSE_IMPL_SPADD_SYMBOLIC_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPADD_SYMBOLIC_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosSparse_spadd_symbolic_impl.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t,
          class b_size_view_t_, class b_lno_view_t, class c_size_view_t_>
struct spadd_symbolic_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPADD_SYMBOLIC_ETI_SPEC_AVAIL(                       \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    MEM_SPACE_TYPE)                                                       \
  template <>                                                             \
  struct spadd_symbolic_eti_spec_avail<                                   \
      KokkosKernels::Experimental::KokkosKernelsHandle<                   \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,               \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {          \
    enum : bool { value = true };                                         \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spadd_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spadd_symbolic_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosBlas::spadd (sparse-sparse matrix addition)

template <class KernelHandle, class a_size_view_t, class a_lno_view_t,
          class b_size_view_t, class b_lno_view_t, class c_size_view_t,
          bool tpl_spec_avail = spadd_symbolic_tpl_spec_avail<
              KernelHandle, a_size_view_t, a_lno_view_t, b_size_view_t,
              b_lno_view_t, c_size_view_t>::value,
          bool eti_spec_avail = spadd_symbolic_eti_spec_avail<
              KernelHandle, a_size_view_t, a_lno_view_t, b_size_view_t,
              b_lno_view_t, c_size_view_t>::value>
struct SPADD_SYMBOLIC {
  static void spadd_symbolic(KernelHandle *handle, a_size_view_t row_mapA,
                             a_lno_view_t entriesA, b_size_view_t row_mapB,
                             b_lno_view_t entriesB, c_size_view_t row_mapC);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

template <class KernelHandle, class a_size_view_t, class a_lno_view_t,
          class b_size_view_t, class b_lno_view_t, class c_size_view_t>
struct SPADD_SYMBOLIC<KernelHandle, a_size_view_t, a_lno_view_t, b_size_view_t,
                      b_lno_view_t, c_size_view_t, false,
                      KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void spadd_symbolic(KernelHandle *handle, a_size_view_t row_mapA,
                             a_lno_view_t entriesA, b_size_view_t row_mapB,
                             b_lno_view_t entriesB, c_size_view_t row_mapC) {
    spadd_symbolic_impl(handle, row_mapA, entriesA, row_mapB, entriesB,
                        row_mapC);
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPADD_SYMBOLIC_ETI_SPEC_DECL(                        \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    MEM_SPACE_TYPE)                                                       \
  extern template struct SPADD_SYMBOLIC<                                  \
      typename KokkosKernels::Experimental::KokkosKernelsHandle<          \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,               \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      false, true>;

#define KOKKOSSPARSE_SPADD_SYMBOLIC_ETI_SPEC_INST(                        \
    SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, \
    MEM_SPACE_TYPE)                                                       \
  template struct SPADD_SYMBOLIC<                                         \
      KokkosKernels::Experimental::KokkosKernelsHandle<                   \
          const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE,       \
          EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,               \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE,                      \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE,                     \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE,                            \
                   Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,             \
      false, true>;

#include <KokkosSparse_spadd_tpl_spec_decl.hpp>
#include <generated_specializations_hpp/KokkosSparse_spadd_symbolic_eti_spec_decl.hpp>

#endif
