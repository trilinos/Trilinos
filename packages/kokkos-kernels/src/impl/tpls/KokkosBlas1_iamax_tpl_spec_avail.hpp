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

#ifndef KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class RV, class XMV, int Xrank = XMV::Rank>
struct iamax_tpl_spec_avail {
  enum : bool { value = false };
};
}
}

namespace KokkosBlas {
namespace Impl {

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
// double
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS( INDEX_TYPE, SCALAR, LAYOUT, MEMSPACE ) \
template<class ExecSpace> \
struct iamax_tpl_spec_avail< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1> { enum : bool { value = true }; };

#if defined (KOKKOSKERNELS_INST_DOUBLE)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS( unsigned long, double,                  Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_FLOAT)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS( unsigned long, float,                   Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS( unsigned long, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_BLAS( unsigned long, Kokkos::complex<float>,  Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
// double
#define KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( INDEX_TYPE, SCALAR, LAYOUT, MEMSPACE ) \
template<class ExecSpace> \
struct iamax_tpl_spec_avail< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::HostSpace, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1> { enum : bool { value = true }; }; \
template<class ExecSpace> \
struct iamax_tpl_spec_avail< \
Kokkos::View<INDEX_TYPE, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::View<const SCALAR*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1> { enum : bool { value = true }; };

#if defined (KOKKOSKERNELS_INST_DOUBLE)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned long, double,                 Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned int, double,                  Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_FLOAT)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned long, float,                  Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned int, float,                   Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned long, Kokkos::complex<double>,Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned int, Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned long, Kokkos::complex<float>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
KOKKOSBLAS1_IAMAX_TPL_SPEC_AVAIL_CUBLAS( unsigned int, Kokkos::complex<float>,  Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

#endif

}
}
#endif
