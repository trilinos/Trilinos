/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels: Linear Algebra and Graph Kernels
//                 Copyright 2016 Sandia Corporation
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





/// \author Kyungjoo Kim (kyukim@sandia.gov)


#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_Vector.hpp"

namespace KokkosKernels {

  namespace Test {
    
    enum { TEST_ADD = 0,
           TEST_MINUS = 1,
           TEST_MULT = 2,
           TEST_DIV = 3 };
    
    template<typename ViewType, int TestID>
    struct Functor {
      ViewType _a, _b, _c;

      KOKKOS_INLINE_FUNCTION
      Functor(ViewType a_, ViewType b_, ViewType c_) 
        : _a(a_), _b(b_), _c(c_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const int i) const {
        switch (TestID) {
        case 0: _c(i) = _a(i) + _b(i); break;
        case 1: _c(i) = _a(i) - _b(i); break;
        case 2: _c(i) = _a(i) * _b(i); break;
        case 3: _c(i) = _a(i) / _b(i); break;
        }
      }
      
    };

    template<typename DeviceSpaceType, typename VectorTagType, int TestID>
    void VectorArithmatic() {
      constexpr int N = 100000;
      
      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      const int iter_begin = -100, iter_end = 100;
      Kokkos::Impl::Timer timer;

      {
        typedef Kokkos::View<ValueType*,HostSpaceType> ViewType;
        ViewType a_host("a_host", N*VectorLength), b_host("b_host", N*VectorLength), c_host("c_host", N*VectorLength);
        for (int k=0;k<N*VectorLength;++k) {
          const int 
            i = k/VectorLength,
            j = k%VectorLength;
          a_host(k) = j + 1;
          b_host(k) = i + 1;
          c_host(k) = 0;
        }
        
        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);
        auto c = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), c_host);
        
        Kokkos::deep_copy(a, a_host);
        Kokkos::deep_copy(b, b_host);
        Kokkos::deep_copy(c, c_host);
        
        {
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, N*VectorLength);
            Kokkos::parallel_for( policy, Functor<ViewType,TestID>(a, b, c) );
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          std::cout << "test = " << TestID << " time(plain) = " << (t/iter_end) << std::endl;
        }
      }

      {
        typedef Vector<VectorTagType> VectorType;
        typedef Kokkos::View<VectorType*,HostSpaceType> ViewType;
        ViewType a_host("a_host", N), b_host("b_host", N), c_host("c_host", N);
        for (int i=0;i<N;++i)
          for (int j=0;j<VectorLength;++j) {
            a_host(i)[j] = j + 1;
            b_host(i)[j] = i + 1;
            c_host(i)[j] = 0;
          }
        
        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);
        auto c = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), c_host);
        
        Kokkos::deep_copy(a, a_host);
        Kokkos::deep_copy(b, b_host);
        Kokkos::deep_copy(c, c_host);
        
        {
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, N);
            Kokkos::parallel_for( policy, Functor<ViewType,TestID>(a, b, c) );
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          std::cout << "test = " << TestID << " time(simd ) = " << (t/iter_end) << std::endl;
        }
      }
    }
  }
}

int main(int argc, char *argv[]) {

  Kokkos::initialize();
  
  std::cout << " Testing SIMD4 double \n";

  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,4>,0>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,4>,1>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,4>,2>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,4>,3>();

  std::cout << " Testing SIMD2 complex<double> \n";
  
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<Kokkos::complex<double> >,2>,0>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<Kokkos::complex<double> >,2>,1>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<Kokkos::complex<double> >,2>,2>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<Kokkos::complex<double> >,2>,3>();

#if defined(__AVX__) || defined(__AVX2__)
  std::cout << " Testing AVX256 double \n";

  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,4>,0>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,4>,1>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,4>,2>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,4>,3>();

  std::cout << " Testing AVX256 complex<double> \n";

  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<Kokkos::complex<double> >,2>,0>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<Kokkos::complex<double> >,2>,1>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<Kokkos::complex<double> >,2>,2>();
  
  // division is not yet implemented
  // KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<Kokkos::complex<double> >,2>,3>();
#endif

  std::cout << " Testing SIMD8 \n";

  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,8>,0>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,8>,1>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,8>,2>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::SIMD<double>,8>,3>();

#if defined(__AVX512F__)
  std::cout << " Testing AVX512 \n";

  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,8>,0>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,8>,1>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,8>,2>();
  KokkosKernels::Test::VectorArithmatic<Kokkos::Serial,KokkosKernels::VectorTag<KokkosKernels::AVX<double>,8>,3>();
#endif
  
  Kokkos::finalize();

  return 0;
}

