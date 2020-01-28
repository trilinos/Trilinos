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

#ifndef KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_HPP_
#define KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_HPP_

namespace KokkosBlas {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AMatrix,
         class BXMV>
struct gesv_tpl_spec_avail {
  enum : bool { value = false };
};

// Generic Host side BLAS (could be MKL or whatever)
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS

#define KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( SCALAR, LAYOUT, MEMSPACE ) \
template<class ExecSpace> \
struct gesv_tpl_spec_avail< \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
     >  { enum : bool { value = true }; };

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( double,                  Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( float,                   Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( Kokkos::complex<float>,  Kokkos::LayoutLeft, Kokkos::HostSpace)
#endif
/*
#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( double,                  Kokkos::LayoutRight, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( float,                   Kokkos::LayoutRight, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( Kokkos::complex<double>, Kokkos::LayoutRight, Kokkos::HostSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_BLAS( Kokkos::complex<float>,  Kokkos::LayoutRight, Kokkos::HostSpace)
#endif
*/
#endif

// MAGMA
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA

#define KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( SCALAR, LAYOUT, MEMSPACE ) \
template<class ExecSpace> \
struct gesv_tpl_spec_avail< \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
     Kokkos::View<SCALAR**, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
     >  { enum : bool { value = true }; };

#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( double,                  Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( float,                   Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( Kokkos::complex<double>, Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTLEFT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( Kokkos::complex<float>,  Kokkos::LayoutLeft, Kokkos::CudaSpace)
#endif

/*
#if defined (KOKKOSKERNELS_INST_DOUBLE) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( double,                 Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_FLOAT) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( float,                  Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_DOUBLE_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( Kokkos::complex<double>,Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif
#if defined (KOKKOSKERNELS_INST_KOKKOS_COMPLEX_FLOAT_) \
 && defined (KOKKOSKERNELS_INST_LAYOUTRIGHT)
 KOKKOSBLAS_GESV_TPL_SPEC_AVAIL_MAGMA( Kokkos::complex<float>, Kokkos::LayoutRight, Kokkos::CudaSpace)
#endif
*/
#endif

}
}

#endif
