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

#ifndef KOKKOSPARSE_SPMV_TPL_SPEC_AVAIL_HPP_
#define KOKKOSPARSE_SPMV_TPL_SPEC_AVAIL_HPP_

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM>
struct spmv_tpl_spec_avail {
  enum : bool { value = false };
};

// cuSPARSE
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

//These versions of cuSPARSE require the ordinal and offset types to be the same.
//For KokkosKernels, this means int/int only.

#define KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, ORDINAL, OFFSET, XL, YL, MEMSPACE) \
template <> \
struct spmv_tpl_spec_avail<const SCALAR, const ORDINAL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET, \
			   const SCALAR*,       XL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, \
			   SCALAR*,             YL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> > { \
  enum : bool { value = true }; \
};

#if (9000 <= CUDA_VERSION)

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int, int, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

//CUDA_VERSION by itself cannot determine whether the generic cuSPARSE API is available:
//cuSPARSE version 10.1.105 does not have the generic API, but it comes with the same CUDA_VERSION (10010) as 10.1.243 which does.
#if defined(CUSPARSE_VERSION) && (CUSPARSE_VERSION >= 10300)

//Can enable int64/size_t.
//TODO: if Nvidia ever supports int/size_t, add that too.

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<float>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  && defined (KOKKOSKERNELS_INST_OFFSET_SIZE_T)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(Kokkos::complex<double>, int64_t, size_t, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaUVMSpace)
#endif

#endif  // CUSPARSE >= 10.3 (nested, implies >= 9.0)
#endif  // CUDA/CUSPARSE >= 9.0?
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

// Specialization struct which defines whether a specialization exists
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM,
         const bool integerScalarType =
           std::is_integral<typename std::decay<AT>::type>::value>
struct spmv_mv_tpl_spec_avail {
  enum : bool { value = false };
};

}}  // KokkosSparse::Impl

#endif // KOKKOSPARSE_SPMV_TPL_SPEC_AVAIL_HPP_
