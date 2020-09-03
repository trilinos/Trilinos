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

#ifndef KOKKOSBLAS1_SCAL_TPL_SPEC_DECL_HPP_
#define KOKKOSBLAS1_SCAL_TPL_SPEC_DECL_HPP_

namespace KokkosBlas {
namespace Impl {

namespace {
  template<class RV, class AV, class XV>
  inline void scal_print_specialization() {
      #ifdef KOKKOSKERNELS_ENABLE_CHECK_SPECIALIZATION
        printf("KokkosBlas1::scal<> TPL Blas specialization for < %s , %s , %s >\n",typeid(RV).name(),typeid(AV).name(),typeid(XV).name());
      #endif
  }
}
}
}

#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include "KokkosBlas_Host_tpl.hpp"

namespace KokkosBlas {
namespace Impl {


#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
double, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef double AV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const double& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_BLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      int N = numElems; \
      int one = 1; \
      HostBlas<double>::scal(N,alpha,R.data(),one);     \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
float, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef float AV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const float& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_BLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      int N = numElems; \
      int one = 1; \
      HostBlas<float>::scal(N,alpha,R.data(),one);     \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::complex<double>, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::complex<double> AV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const Kokkos::complex<double>& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_BLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      int N = numElems; \
      int one = 1; \
      const std::complex<double> alpha_val = alpha;     \
      HostBlas<std::complex<double> >::scal\
        (N,alpha_val,              \
         reinterpret_cast<std::complex<double>*>(R.data()), one);       \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::complex<float>, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::complex<float> AV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const Kokkos::complex<float>& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_BLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      int N = numElems; \
      int one = 1; \
      const std::complex<float> alpha_val = alpha;     \
      HostBlas<std::complex<float> >::scal\
        (N,alpha_val,               \
         reinterpret_cast<std::complex<float>*>(R.data()), one);        \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_BLAS( Kokkos::LayoutLeft, Kokkos::HostSpace, false)

}
}

#endif

// cuBLAS
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include<KokkosBlas_tpl_spec.hpp>

namespace KokkosBlas {
namespace Impl {

#define KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
double, \
Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef double AV; \
  typedef Kokkos::View<const double*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const double& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_CUBLAS,double]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasDscal(s.handle, N, &alpha, R.data(), one); \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
float, \
Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef float AV; \
  typedef Kokkos::View<const float*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const float& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_CUBLAS,float]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasSscal(s.handle, N, &alpha, R.data(), one); \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::complex<double>, \
Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::complex<double> AV; \
  typedef Kokkos::View<const Kokkos::complex<double>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const Kokkos::complex<double>& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_CUBLAS,complex<double>]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasZscal(s.handle, N, reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<cuDoubleComplex*>(R.data()), one); \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

#define KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS( LAYOUT, MEMSPACE, ETI_SPEC_AVAIL ) \
template<class ExecSpace> \
struct Scal< \
Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
Kokkos::complex<float>, \
Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
             Kokkos::MemoryTraits<Kokkos::Unmanaged> >, \
1,true, ETI_SPEC_AVAIL > { \
  \
  typedef Kokkos::View<Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > RV; \
  typedef Kokkos::complex<float> AV; \
  typedef Kokkos::View<const Kokkos::complex<float>*, LAYOUT, Kokkos::Device<ExecSpace, MEMSPACE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > XV; \
  typedef typename XV::size_type size_type; \
  \
  static void scal (const RV& R, const Kokkos::complex<float>& alpha, const XV& X) \
  { \
    Kokkos::Profiling::pushRegion("KokkosBlas::scal[TPL_CUBLAS,complex<float>]"); \
    const size_type numElems = X.extent(0); \
    if ((numElems < static_cast<size_type> (INT_MAX)) && (R.data() == X.data())) { \
      scal_print_specialization<RV,AV,XV>(); \
      const int N = static_cast<int> (numElems); \
      constexpr int one = 1; \
      KokkosBlas::Impl::CudaBlasSingleton & s = KokkosBlas::Impl::CudaBlasSingleton::singleton(); \
      cublasCscal(s.handle, N, reinterpret_cast<const cuComplex*>(&alpha), reinterpret_cast<cuComplex*>(R.data()), one); \
    } else { \
      Scal<RV,AV,XV,1,false,ETI_SPEC_AVAIL>::scal(R,alpha,X); \
    } \
    Kokkos::Profiling::popRegion(); \
  } \
};

KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaSpace, false)

#if defined (KOKKOS_ENABLE_CUDA_UVM)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_DSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_SSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_ZSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)

KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, true)
KOKKOSBLAS1_CSCAL_TPL_SPEC_DECL_CUBLAS( Kokkos::LayoutLeft, Kokkos::CudaUVMSpace, false)
#endif

}
}

#endif

#endif
