// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_HPP_
#define KOKKOSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_HPP_

namespace KokkosSparse {
namespace Impl {

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#define KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE(TOKEN, KOKKOS_SCALAR_TYPE, TPL_SCALAR_TYPE, ORDINAL_TYPE,    \
                                                          OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,   \
                                                          ETI_SPEC_AVAIL)                                              \
  template <>                                                                                                          \
  struct SPADD_NUMERIC<                                                                                                \
      EXEC_SPACE_TYPE,                                                                                                 \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,                          \
                                                       const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,      \
                                                       MEM_SPACE_TYPE>,                                                \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using kernelhandle_t = KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,     \
                                                                            const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, \
                                                                            MEM_SPACE_TYPE, MEM_SPACE_TYPE>;           \
    using rowmap_view_t =                                                                                              \
        Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using non_const_rowmap_view_t =                                                                                    \
        Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using colidx_view_t =                                                                                              \
        Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,               \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using non_const_colidx_view_t =                                                                                    \
        Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using scalar_view_t =                                                                                              \
        Kokkos::View<const KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,         \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using non_const_scalar_view_t =                                                                                    \
        Kokkos::View<KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,               \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    static void spadd_numeric(const EXEC_SPACE_TYPE &exec, kernelhandle_t *handle, ORDINAL_TYPE m, ORDINAL_TYPE n,     \
                              const KOKKOS_SCALAR_TYPE alpha, rowmap_view_t rowmapA, colidx_view_t colidxA,            \
                              scalar_view_t valuesA, const KOKKOS_SCALAR_TYPE beta, rowmap_view_t rowmapB,             \
                              colidx_view_t colidxB, scalar_view_t valuesB, rowmap_view_t rowmapC,                     \
                              non_const_colidx_view_t colidxC, non_const_scalar_view_t valuesC) {                      \
      Kokkos::Profiling::pushRegion("KokkosSparse::spadd_numeric[TPL_CUSPARSE," +                                      \
                                    KokkosKernels::ArithTraits<KOKKOS_SCALAR_TYPE>::name() + "]");                     \
                                                                                                                       \
      auto addHandle   = handle->get_spadd_handle();                                                                   \
      auto &cuspData   = addHandle->cusparseData;                                                                      \
      auto &cuspHandle = KokkosKernels::Impl::CusparseSingleton::singleton().cusparseHandle;                           \
      cusparsePointerMode_t oldPtrMode;                                                                                \
                                                                                                                       \
      KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseSetStream(cuspHandle, exec.cuda_stream()));                         \
      KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseGetPointerMode(cuspHandle, &oldPtrMode));                           \
      KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(                                                                            \
          cusparseSetPointerMode(cuspHandle, CUSPARSE_POINTER_MODE_HOST)); /* alpha, beta on host*/                    \
      OFFSET_TYPE nnzA = colidxA.extent(0);                                                                            \
      OFFSET_TYPE nnzB = colidxB.extent(0);                                                                            \
      KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparse##TOKEN##csrgeam2(                                                  \
          cuspHandle, m, n, reinterpret_cast<const TPL_SCALAR_TYPE *>(&alpha), cuspData.descrA, nnzA,                  \
          reinterpret_cast<const TPL_SCALAR_TYPE *>(valuesA.data()), rowmapA.data(), colidxA.data(),                   \
          reinterpret_cast<const TPL_SCALAR_TYPE *>(&beta), cuspData.descrB, nnzB,                                     \
          reinterpret_cast<const TPL_SCALAR_TYPE *>(valuesB.data()), rowmapB.data(), colidxB.data(), cuspData.descrC,  \
          reinterpret_cast<TPL_SCALAR_TYPE *>(valuesC.data()), const_cast<OFFSET_TYPE *>(rowmapC.data()),              \
          colidxC.data(), cuspData.workspace));                                                                        \
      KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(cuspHandle, oldPtrMode));                            \
      KOKKOSSPARSE_IMPL_CUSPARSE_SAFE_CALL(cusparseSetStream(cuspHandle, NULL));                                       \
                                                                                                                       \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE_EXT(ETI_SPEC_AVAIL)                                      \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE(S, float, float, int, int, Kokkos::LayoutLeft, Kokkos::Cuda,   \
                                                    Kokkos::CudaSpace, ETI_SPEC_AVAIL)                             \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE(D, double, double, int, int, Kokkos::LayoutLeft, Kokkos::Cuda, \
                                                    Kokkos::CudaSpace, ETI_SPEC_AVAIL)                             \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE(C, Kokkos::complex<float>, cuComplex, int, int,                \
                                                    Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace,           \
                                                    ETI_SPEC_AVAIL)                                                \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE(Z, Kokkos::complex<double>, cuDoubleComplex, int, int,         \
                                                    Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace,           \
                                                    ETI_SPEC_AVAIL)

KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE_EXT(true)
KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_CUSPARSE_EXT(false)
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE

#define KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE(TOKEN, KOKKOS_SCALAR_TYPE, TPL_SCALAR_TYPE, ORDINAL_TYPE,   \
                                                           OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,  \
                                                           ETI_SPEC_AVAIL)                                             \
  template <>                                                                                                          \
  struct SPADD_NUMERIC<                                                                                                \
      EXEC_SPACE_TYPE,                                                                                                 \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,                          \
                                                       const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE,      \
                                                       MEM_SPACE_TYPE>,                                                \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                  \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                       \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      Kokkos::View<KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>,                                                           \
      true, ETI_SPEC_AVAIL> {                                                                                          \
    using kernelhandle_t = KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE,     \
                                                                            const KOKKOS_SCALAR_TYPE, EXEC_SPACE_TYPE, \
                                                                            MEM_SPACE_TYPE, MEM_SPACE_TYPE>;           \
    using rowmap_view_t =                                                                                              \
        Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using non_const_rowmap_view_t =                                                                                    \
        Kokkos::View<OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                      \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using colidx_view_t =                                                                                              \
        Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,               \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using non_const_colidx_view_t =                                                                                    \
        Kokkos::View<ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                     \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using scalar_view_t =                                                                                              \
        Kokkos::View<const KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,         \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    using non_const_scalar_view_t =                                                                                    \
        Kokkos::View<KOKKOS_SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,               \
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;                                                         \
    static void spadd_numeric(const EXEC_SPACE_TYPE &exec, kernelhandle_t *handle, ORDINAL_TYPE m, ORDINAL_TYPE n,     \
                              const KOKKOS_SCALAR_TYPE alpha, rowmap_view_t rowmapA, colidx_view_t colidxA,            \
                              scalar_view_t valuesA, const KOKKOS_SCALAR_TYPE beta, rowmap_view_t rowmapB,             \
                              colidx_view_t colidxB, scalar_view_t valuesB, rowmap_view_t rowmapC,                     \
                              non_const_colidx_view_t colidxC, non_const_scalar_view_t valuesC) {                      \
      Kokkos::Profiling::pushRegion("KokkosSparse::spadd_numeric[TPL_ROCSPARSE," +                                     \
                                    KokkosKernels::ArithTraits<KOKKOS_SCALAR_TYPE>::name() + "]");                     \
                                                                                                                       \
      auto addHandle    = handle->get_spadd_handle();                                                                  \
      auto &rocData     = addHandle->rocsparseData;                                                                    \
      auto &rocspHandle = KokkosKernels::Impl::RocsparseSingleton::singleton().rocsparseHandle;                        \
      rocsparse_pointer_mode oldPtrMode;                                                                               \
                                                                                                                       \
      KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(rocspHandle, exec.hip_stream()));                     \
      KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_get_pointer_mode(rocspHandle, &oldPtrMode));                     \
      KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(                                                                           \
          rocsparse_set_pointer_mode(rocspHandle, rocsparse_pointer_mode_host)); /* alpha, beta on host*/              \
      OFFSET_TYPE nnzA = colidxA.extent(0);                                                                            \
      OFFSET_TYPE nnzB = colidxB.extent(0);                                                                            \
      KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_##TOKEN##csrgeam(                                                \
          rocspHandle, m, n, reinterpret_cast<const TPL_SCALAR_TYPE *>(&alpha), rocData.descrA, nnzA,                  \
          reinterpret_cast<const TPL_SCALAR_TYPE *>(valuesA.data()), rowmapA.data(), colidxA.data(),                   \
          reinterpret_cast<const TPL_SCALAR_TYPE *>(&beta), rocData.descrB, nnzB,                                      \
          reinterpret_cast<const TPL_SCALAR_TYPE *>(valuesB.data()), rowmapB.data(), colidxB.data(), rocData.descrC,   \
          reinterpret_cast<TPL_SCALAR_TYPE *>(valuesC.data()), const_cast<OFFSET_TYPE *>(rowmapC.data()),              \
          colidxC.data()));                                                                                            \
      KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_set_pointer_mode(rocspHandle, oldPtrMode));                      \
      KOKKOSSPARSE_IMPL_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(rocspHandle, NULL));                                  \
                                                                                                                       \
      Kokkos::Profiling::popRegion();                                                                                  \
    }                                                                                                                  \
  };

#define KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE_EXT(ETI_SPEC_AVAIL)                                       \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE(s, float, float, int, int, Kokkos::LayoutLeft, Kokkos::HIP,     \
                                                     Kokkos::HIPSpace, ETI_SPEC_AVAIL)                               \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE(d, double, double, int, int, Kokkos::LayoutLeft, Kokkos::HIP,   \
                                                     Kokkos::HIPSpace, ETI_SPEC_AVAIL)                               \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE(c, Kokkos::complex<float>, rocsparse_float_complex, int, int,   \
                                                     Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace,              \
                                                     ETI_SPEC_AVAIL)                                                 \
  KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE(z, Kokkos::complex<double>, rocsparse_double_complex, int, int, \
                                                     Kokkos::LayoutLeft, Kokkos::HIP, Kokkos::HIPSpace,              \
                                                     ETI_SPEC_AVAIL)

KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE_EXT(true)
KOKKOSSPARSE_SPADD_NUMERIC_TPL_SPEC_DECL_ROCSPARSE_EXT(false)
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#endif
