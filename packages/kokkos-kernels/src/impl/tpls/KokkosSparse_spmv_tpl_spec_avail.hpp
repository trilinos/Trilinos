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
#if defined (KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && (9000 <= CUDA_VERSION)

#define KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(SCALAR, XL, YL, MEMSPACE) \
template <> \
struct spmv_tpl_spec_avail<const SCALAR, const int, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged>, const int, \
			   const SCALAR*,       XL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>, \
			   SCALAR*,             YL, Kokkos::Device<Kokkos::Cuda, MEMSPACE>, Kokkos::MemoryTraits<Kokkos::Unmanaged> > { \
  enum : bool { value = true }; \
};

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTLEFT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_FLOAT) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(float, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
  && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT) \
  && defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  && defined (KOKKOSKERNELS_INST_OFFSET_INT)
  KOKKOSSPARSE_SPMV_TPL_SPEC_AVAIL_CUSPARSE(double, Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif

#endif

// Specialization struct which defines whether a specialization exists
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM,
         const bool integerScalarType =
           std::is_integral<typename std::decay<AT>::type>::value>
struct spmv_mv_tpl_spec_avail {
  enum : bool { value = false };
};

}
}

#endif // KOKKOSPARSE_SPMV_TPL_SPEC_AVAIL_HPP_
