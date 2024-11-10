//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSPARSE_SPADD_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPADD_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
//
template <class ExecSpace, class KernelHandle, class a_size_view_t, class a_lno_view_t, class b_size_view_t,
          class b_lno_view_t, class c_size_view_t>
struct spadd_symbolic_tpl_spec_avail {
  enum : bool { value = false };
};

template <class ExecSpace, class KernelHandle, class a_size_view_t, class a_lno_view_t, class a_scalar_view_t,
          class b_size_view_t, class b_lno_view_t, class b_scalar_view_t, class c_size_view_t, class c_lno_view_t,
          class c_scalar_view_t>
struct spadd_numeric_tpl_spec_avail {
  enum : bool { value = false };
};

#define KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,          \
                                                   EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                              \
  template <>                                                                                                    \
  struct spadd_symbolic_tpl_spec_avail<                                                                          \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

#define KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,           \
                                                  EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                               \
  template <>                                                                                                    \
  struct spadd_numeric_tpl_spec_avail<                                                                           \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

#define KOKKOSSPARSE_SPADD_TPL_SPEC_AVAIL(ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_AVAIL(float, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,       \
                                             MEM_SPACE_TYPE)                                                       \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_AVAIL(double, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,      \
                                             MEM_SPACE_TYPE)                                                       \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_AVAIL(Kokkos::complex<float>, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,       \
                                             EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                                      \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_AVAIL(Kokkos::complex<double>, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,      \
                                             EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                                      \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_AVAIL(float, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,        \
                                            MEM_SPACE_TYPE)                                                        \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_AVAIL(double, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE,       \
                                            MEM_SPACE_TYPE)                                                        \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_AVAIL(Kokkos::complex<float>, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,        \
                                            EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                                       \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_AVAIL(Kokkos::complex<double>, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,       \
                                            EXEC_SPACE_TYPE, MEM_SPACE_TYPE)

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
KOKKOSSPARSE_SPADD_TPL_SPEC_AVAIL(int, int, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
KOKKOSSPARSE_SPADD_TPL_SPEC_AVAIL(rocsparse_int, rocsparse_int, Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace)
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#endif
