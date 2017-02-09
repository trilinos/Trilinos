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

#include <iomanip>

#if defined(__KOKKOSKERNELS_INTEL_MKL__)
#include "mkl.h"
#endif

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_Vector.hpp"

#include "KokkosKernels_Gemm_Decl.hpp"
#include "KokkosKernels_Gemm_Serial_Impl.hpp"
#include "KokkosKernels_Gemm_Team_Impl.hpp"

namespace KokkosKernels {

  namespace Test {

#define FLOP_MUL 1.0
#define FLOP_ADD 1.0

    double FlopCount(int mm, int nn, int kk) {
      double m = (double)mm;    double n = (double)nn;    double k = (double)kk;
      return (FLOP_MUL*(m*n*k) +
              FLOP_ADD*(m*n*k));
    }

    template<int BlkSize, typename DeviceSpaceType, typename VectorTagType, typename AlgoTagType>
    void Gemm(const int N) {
      typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;
      //typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;

      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      double flop = (N*VectorLength)*FlopCount(BlkSize,BlkSize,BlkSize);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> cref;

      ///
      /// Reference version using MKL DGEMM
      ///
#if defined(__KOKKOSKERNELS_INTEL_MKL__)
      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N*VectorLength, BlkSize, BlkSize),
          b("b", N*VectorLength, BlkSize, BlkSize),
          c("c", N*VectorLength, BlkSize, BlkSize);

        for (int k=0;k<N*VectorLength;++k) {
          const int
            k0 = k/VectorLength,
            k1 = k%VectorLength;
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              a(k, i, j) = k1 + 1;
              b(k, i, j) = k0 + 1;
              c(k, i, j) = 0;
            }
        }

        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);

          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());
                
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            BlkSize, BlkSize, BlkSize,
                            1.0,
                            (double*)aa.data(), aa.stride_0(),
                            (double*)bb.data(), bb.stride_0(),
                            1.0,
                            (double*)cc.data(), cc.stride_0());
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          std::cout << std::setw(12) << "MKL DGEMM"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << std::endl;

          cref = c;
        }
      }

#if defined(__KOKKOSKERNELS_INTEL_MKL_BATCHED__)
      {
        typedef Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> ViewType;
        ViewType
          a("a", N*VectorLength, BlkSize, BlkSize),
          b("b", N*VectorLength, BlkSize, BlkSize),
          c("c", N*VectorLength, BlkSize, BlkSize);

        for (int k=0;k<N*VectorLength;++k) {
          const int
            k0 = k/VectorLength,
            k1 = k%VectorLength;
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              a(k, i, j) = k1 + 1;
              b(k, i, j) = k0 + 1;
              c(k, i, j) = 0;
            }
        }

        ValueType
          *aa[N*VectorLength],
          *bb[N*VectorLength],
          *cc[N*VectorLength];

        for (int k=0;k<N*VectorLength;++k) {
          aa[k] = &a(k, 0, 0);
          bb[k] = &b(k, 0, 0);
          cc[k] = &c(k, 0, 0);
        }

        {
          double t = 0;

          MKL_INT blksize[1] = { BlkSize };
          MKL_INT lda[1] = { a.stride_1() };
          MKL_INT ldb[1] = { b.stride_1() };
          MKL_INT ldc[1] = { c.stride_1() };

          CBLAS_TRANSPOSE transA[1] = { CblasNoTrans };
          CBLAS_TRANSPOSE transB[1] = { CblasNoTrans };

          double one[1] = { 1.0 };
          MKL_INT size_per_grp[1] = { N*VectorLength };

          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            cblas_dgemm_batch(CblasRowMajor, 
                              transA,
                              transB,
                              blksize, blksize, blksize, 
                              one,
                              (const double**)aa, lda,
                              (const double**)bb, ldb,
                              one,
                              cc, ldc,
                              1, size_per_grp);

            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c(i,j,k));

          std::cout << std::setw(12) << "MKL Batch"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }
#endif
#endif

      ///
      /// Plain version (comparable to micro BLAS version)
      ///
      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a_host("a_host", N*VectorLength, BlkSize, BlkSize),
          b_host("b_host", N*VectorLength, BlkSize, BlkSize),
          c_host("c_host", N*VectorLength, BlkSize, BlkSize);

        for (int k=0;k<N*VectorLength;++k) {
          const int
            k0 = k/VectorLength,
            k1 = k%VectorLength;
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              a_host(k, i, j) = k1 + 1;
              b_host(k, i, j) = k0 + 1;
              c_host(k, i, j) = 0;
            }
        }

        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);
        auto c = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), c_host);

        Kokkos::deep_copy(a, a_host);
        Kokkos::deep_copy(b, b_host);
        Kokkos::deep_copy(c, c_host);

        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
          
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());
                
                KokkosKernels::Serial::
                  Gemm<Trans::NoTranspose,Trans::NoTranspose,AlgoTagType>::
                  invoke(1.0, aa, bb, 1.0, cc);
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          Kokkos::deep_copy(c_host, c);
          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c_host(i,j,k));

          std::cout << std::setw(12) << "Plain"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }

      ///
      /// Serial SIMD with appropriate data layout
      ///
      {
        typedef Vector<VectorTagType> VectorType;
        Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
          a_host("a_host", N, BlkSize, BlkSize),
          b_host("b_host", N, BlkSize, BlkSize),
          c_host("c_host", N, BlkSize, BlkSize);

        for (int k0=0;k0<N;++k0)
          for (int k1=0;k1<VectorLength;++k1)
            for (int i=0;i<BlkSize;++i)
              for (int j=0;j<BlkSize;++j) {
                a_host(k0, i, j)[k1] = k1 + 1;
                b_host(k0, i, j)[k1] = k0 + 1;
                c_host(k0, i, j)[k1] = 0;
              }

        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);
        auto c = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), c_host);

        Kokkos::deep_copy(a, a_host);
        Kokkos::deep_copy(b, b_host);
        Kokkos::deep_copy(c, c_host);

        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N);
          
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());
                
                KokkosKernels::Serial::
                  Gemm<Trans::NoTranspose,Trans::NoTranspose,AlgoTagType>::
                  invoke(1.0, aa, bb, 1.0, cc);
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          Kokkos::deep_copy(c_host, c);
          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c_host(i/VectorLength,j,k)[i%VectorLength]);

          std::cout << std::setw(12) << "Serial SIMD"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }

      ///
      /// Team SIMD with appropriate data layout
      ///

      {
        typedef Vector<VectorTagType> VectorType;
        Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
          a_host("a_host", N, BlkSize, BlkSize),
          b_host("b_host", N, BlkSize, BlkSize),
          c_host("c_host", N, BlkSize, BlkSize);

        for (int k0=0;k0<N;++k0)
          for (int k1=0;k1<VectorLength;++k1)
            for (int i=0;i<BlkSize;++i)
              for (int j=0;j<BlkSize;++j) {
                a_host(k0, i, j)[k1] = k1 + 1;
                b_host(k0, i, j)[k1] = k0 + 1;
                c_host(k0, i, j)[k1] = 0;
              }

        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);
        auto c = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), c_host);

        Kokkos::deep_copy(a, a_host);
        Kokkos::deep_copy(b, b_host);
        Kokkos::deep_copy(c, c_host);

        {
          double t = 0;
          typedef Kokkos::TeamPolicy<DeviceSpaceType,ScheduleType> policy_type;
          typedef typename policy_type::member_type member_type;
          const policy_type policy(N, Kokkos::AUTO, VectorLength);
          
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const member_type &member) {
                const int k = member.league_rank();

                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());
                
                KokkosKernels::Team::
                  Gemm<member_type,Trans::NoTranspose,Trans::NoTranspose,AlgoTagType>::
                  invoke(member, 1.0, aa, bb, 1.0, cc);
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          Kokkos::deep_copy(c_host, c);
          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c_host(i/VectorLength,j,k)[i%VectorLength]);

          std::cout << std::setw(12) << "Team SIMD"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " diff to ref = " << diff
                    << std::endl << std::endl;
        }
      }
    }
  
  }
}

using namespace KokkosKernels;

template<typename VectorType,
         typename AlgoTagType>
void run(const int N) {
  typedef typename VectorType::exec_space ExecSpace;

  std::cout << "ExecSpace::  "; 
  if (std::is_same<ExecSpace,Kokkos::Serial>::value) 
    std::cout << "Kokkos::Serial " << std::endl;
  else 
    ExecSpace::print_configuration(std::cout, false);

  // Test::Gemm< 4, ExecSpace,VectorType,AlgoTagType>();
  // Test::Gemm< 8, ExecSpace,VectorType,AlgoTagType>();
  // Test::Gemm<16, ExecSpace,VectorType,AlgoTagType>();
  // Test::Gemm<20, ExecSpace,VectorType,AlgoTagType>();
  // Test::Gemm<32, ExecSpace,VectorType,AlgoTagType>();
  // Test::Gemm<64, ExecSpace,VectorType,AlgoTagType>();

  Test::Gemm< 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemm< 9, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemm<15, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemm<20, ExecSpace,VectorType,AlgoTagType>(N);
}

int main(int argc, char *argv[]) {
  
  Kokkos::initialize(argc, argv);

  const int ntest = 6;
  const int N[6] = { 256, 512, 768, 1024, 1280, 1536 };
  //const int N[1] = { 512 };

  {    
    /// 256 test
#if (defined(__AVX2__) || defined(__AVX__)) && !defined(__AVX512F__)
    for (int i=0;i<ntest;++i) {
      std::cout << " N = " << N[i] << std::endl;
      
      // std::cout << "\n Testing SIMD4 and Algo::Gemm::Triple\n";
      // run<VectorTag<SIMD<double>,4>,Algo::Gemm::Triple>();
      
      // std::cout << "\n Testing AVX256 and Algo::Gemm::Triple\n";
      // run<VectorTag<AVX<double>,4>,Algo::Gemm::Triple>();
      
      //std::cout << "\n Testing SIMD4 and Algo::Gemm::Blocked\n";
      //run<VectorTag<SIMD<double>,4>,Algo::Gemm::Blocked>(N[i]/4);
      
      std::cout << "\n Testing AVX256 and Algo::Gemm::Blocked\n";
      run<VectorTag<AVX<double>,4>,Algo::Gemm::Blocked>(N[i]/4);
    }
#endif
    
    /// 512 test
#if defined(__AVX512F__)
    for (int i=0;i<ntest;++i) {
      std::cout << " N = " << N[i] << std::endl;

      //  std::cout << "\n Testing SIMD8 and Algo::Gemm::Triple\n";
      //run<VectorTag<SIMD<double>,8>,Algo::Gemm::Triple>();
      
      //std::cout << "\n Testing AVX512 and Algo::Gemm::Triple\n";
      //run<VectorTag<AVX<double>,8>,Algo::Gemm::Triple>();
      
      //std::cout << "\n Testing SIMD8 and Algo::Gemm::Blocked\n";
      //run<VectorTag<SIMD<double>,8>,Algo::Gemm::Blocked>(N[i]/8);
      
      std::cout << "\n Testing AVX512 and Algo::Gemm::Blocked\n";
      run<VectorTag<AVX<double>,8>,Algo::Gemm::Blocked>(N[i]/8);
    }
#endif
  }

  Kokkos::finalize();

  return 0;
}
