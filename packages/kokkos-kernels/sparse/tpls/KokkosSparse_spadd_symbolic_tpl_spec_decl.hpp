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

#ifndef KOKKOSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_HPP_

namespace KokkosSparse {
namespace Impl {

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#define KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE(TOKEN, KOKKOS_SCALAR_TYPE, TPL_SCALAR_TYPE, ORDINAL_TYPE,   \
                                                           OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,  \
                                                           ETI_SPEC_AVAIL)                                             \
  template <>                                                                                                          \
  struct SPADD_SYMBOLIC<                                                                                               \
      EXEC_SPACE_TYPE,                                                                                                 \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,                          \
                                                       const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,      \
                                                       MEM_SPACE_TYPE>,                                                \
      Kokkos::View<const OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using kernelhandle_t = KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,     \
                                                                            const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, \
                                                                            MEM_SPACE_TYPE, MEM_SPACE_TYPE>;           \
    using rowmap_view_t =                                                                                              \
        Kokkos::View<const OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                                        \
    using non_const_rowmap_view_t =                                                                                    \
        Kokkos::View<OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                                        \
    using colidx_view_t =                                                                                              \
        Kokkos::View<const ORDINAL_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                                        \
    static void spadd_symbolic(const EXEC_SPACE_TYPE& exec, kernelhandle_t* handle, const ORDINAL_TYPE m,              \
                               const ORDINAL_TYPE n, rowmap_view_t rowmapA, colidx_view_t colidxA,                     \
                               rowmap_view_t rowmapB, colidx_view_t colidxB, non_const_rowmap_view_t rowmapC) {        \
      Kokkos::Profiling::pushRegion("KokkosSparse::spadd_symbolic[TPL_CUSPARSE," +                                     \
                                    Kokkos::ArithTraits<KOKKOS_SCALAR_TYPE>::name() + "]");                            \
                                                                                                                       \
      auto addHandle   = handle->get_spadd_handle();                                                                   \
      auto& cuspData   = addHandle->cusparseData;                                                                      \
      auto& cuspHandle = KokkosKernels::Impl::CusparseSingleton::singleton().cusparseHandle;                           \
                                                                                                                       \
      /* Not easy to init 'one' for cuda complex, so we don't init it. Anyway,                                         \
       * the uninit'ed var won't affect C's pattern.                                                                   \
       */                                                                                                              \
      TPL_SCALAR_TYPE one;                                                                                             \
      size_t nbytes;                                                                                                   \
      OFFSET_TYPE nnzA = colidxA.extent(0);                                                                            \
      OFFSET_TYPE nnzB = colidxB.extent(0);                                                                            \
      OFFSET_TYPE nnzC = 0;                                                                                            \
                                                                                                                       \
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(cuspHandle, exec.cuda_stream()));                                    \
                                                                                                                       \
      /* https://docs.nvidia.com/cuda/cusparse/index.html#cusparsecreatematdescr                                       \
       It sets the fields MatrixType and IndexBase to the default values                                               \
       CUSPARSE_MATRIX_TYPE_GENERAL and CUSPARSE_INDEX_BASE_ZERO,                                                      \
       respectively, while leaving other fields uninitialized. */                                                      \
                                                                                                                       \
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&cuspData.descrA));                                             \
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&cuspData.descrB));                                             \
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&cuspData.descrC));                                             \
      KOKKOS_CUSPARSE_SAFE_CALL(cusparse##TOKEN##csrgeam2_bufferSizeExt(                                               \
          cuspHandle, m, n, &one, cuspData.descrA, nnzA, NULL, rowmapA.data(), colidxA.data(), &one, cuspData.descrB,  \
          nnzB, NULL, rowmapB.data(), colidxB.data(), cuspData.descrC, NULL, rowmapC.data(), NULL, &nbytes));          \
      cuspData.nbytes    = nbytes;                                                                                     \
      cuspData.workspace = Kokkos::kokkos_malloc<MEM_SPACE_TYPE>(nbytes);                                              \
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseXcsrgeam2Nnz(                                                                  \
          cuspHandle, m, n, cuspData.descrA, nnzA, rowmapA.data(), colidxA.data(), cuspData.descrB, nnzB,              \
          rowmapB.data(), colidxB.data(), cuspData.descrC, rowmapC.data(), &nnzC, cuspData.workspace));                \
      addHandle->set_c_nnz(nnzC);                                                                                      \
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetStream(cuspHandle, NULL));                                                  \
                                                                                                                       \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE_EXT(ETI_SPEC_AVAIL)                                      \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE(S, float, float, int, int, Kokkos::LayoutLeft, Kokkos::Cuda,   \
                                                     Kokkos::CudaSpace, ETI_SPEC_AVAIL)                             \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE(D, double, double, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, \
                                                     Kokkos::CudaSpace, ETI_SPEC_AVAIL)                             \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE(C, Kokkos::complex<float>, cuComplex, int, int,                \
                                                     Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace,           \
                                                     ETI_SPEC_AVAIL)                                                \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE(Z, Kokkos::complex<double>, cuDoubleComplex, int, int,         \
                                                     Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace,           \
                                                     ETI_SPEC_AVAIL)

KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE_EXT(true)
KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_CUSPARSE_EXT(false)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#define KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE(                                                           \
    KOKKOS_SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE, ETI_SPEC_AVAIL)       \
  template <>                                                                                                          \
  struct SPADD_SYMBOLIC<                                                                                               \
      EXEC_SPACE_TYPE,                                                                                                 \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,                          \
                                                       const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,      \
                                                       MEM_SPACE_TYPE>,                                                \
      Kokkos::View<const OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                   \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<const ORDINAL_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      Kokkos::View<OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                         \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                          \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using kernelhandle_t = KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,     \
                                                                            const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, \
                                                                            MEM_SPACE_TYPE, MEM_SPACE_TYPE>;           \
    using rowmap_view_t =                                                                                              \
        Kokkos::View<const OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                                        \
    using non_const_rowmap_view_t =                                                                                    \
        Kokkos::View<OFFSET_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                                        \
    using colidx_view_t =                                                                                              \
        Kokkos::View<const ORDINAL_TYPE*, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;                                                        \
    static void spadd_symbolic(const EXEC_SPACE_TYPE& exec, kernelhandle_t* handle, const ORDINAL_TYPE m,              \
                               const ORDINAL_TYPE n, rowmap_view_t rowmapA, colidx_view_t colidxA,                     \
                               rowmap_view_t rowmapB, colidx_view_t colidxB, non_const_rowmap_view_t rowmapC) {        \
      Kokkos::Profiling::pushRegion("KokkosSparse::spadd_symbolic[TPL_ROCSPARSE," +                                    \
                                    Kokkos::ArithTraits<KOKKOS_SCALAR_TYPE>::name() + "]");                            \
                                                                                                                       \
      auto addHandle    = handle->get_spadd_handle();                                                                  \
      auto& rocData     = addHandle->rocsparseData;                                                                    \
      auto& rocspHandle = KokkosKernels::Impl::RocsparseSingleton::singleton().rocsparseHandle;                        \
      OFFSET_TYPE nnzA  = colidxA.extent(0);                                                                           \
      OFFSET_TYPE nnzB  = colidxB.extent(0);                                                                           \
      OFFSET_TYPE nnzC  = 0;                                                                                           \
                                                                                                                       \
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_stream(rocspHandle, exec.hip_stream()));                           \
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&rocData.descrA));                                    \
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&rocData.descrB));                                    \
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&rocData.descrC));                                    \
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_csrgeam_nnz(rocspHandle, m, n, rocData.descrA, nnzA, rowmapA.data(),   \
                                                            colidxA.data(), rocData.descrB, nnzB, rowmapB.data(),      \
                                                            colidxB.data(), rocData.descrC, rowmapC.data(), &nnzC));   \
      addHandle->set_c_nnz(nnzC);                                                                                      \
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_set_stream(rocspHandle, NULL));                                        \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE_EXT(ETI_SPEC_AVAIL)                                 \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE(float, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft,  \
                                                      Kokkos::HIP, Kokkos::HIPSpace, ETI_SPEC_AVAIL)            \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE(double, rocsparse_int, rocsparse_int, Kokkos::LayoutLeft, \
                                                      Kokkos::HIP, Kokkos::HIPSpace, ETI_SPEC_AVAIL)            \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE(Kokkos::complex<float>, rocsparse_int, rocsparse_int,     \
                                                      Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace,        \
                                                      ETI_SPEC_AVAIL)                                           \
  KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE(Kokkos::complex<double>, rocsparse_int, rocsparse_int,    \
                                                      Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace,        \
                                                      ETI_SPEC_AVAIL)

KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE_EXT(true)
KOKKOSSPARSE_SPADD_SYMBOLIC_TPL_SPEC_DECL_ROCSPARSE_EXT(false)
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#endif
