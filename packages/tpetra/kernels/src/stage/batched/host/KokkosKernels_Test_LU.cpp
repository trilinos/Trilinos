/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>
#include "mkl.h"

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_Vector.hpp"

#include "KokkosKernels_LU_Serial_Decl.hpp"
#include "KokkosKernels_LU_Serial_Impl.hpp"

namespace KokkosKernels {

  namespace Test {

#define FLOP_MUL 1.0
#define FLOP_ADD 1.0

    double FlopCount(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      if (m > n)
        return (FLOP_MUL*(0.5*m*n*n-(1.0/6.0)*n*n*n+0.5*m*n-0.5*n*n+(2.0/3.0)*n) +
                FLOP_ADD*(0.5*m*n*n-(1.0/6.0)*n*n*n-0.5*m*n+        (1.0/6.0)*n));
      else
        return (FLOP_MUL*(0.5*n*m*m-(1.0/6.0)*m*m*m+0.5*n*m-0.5*m*m+(2.0/3.0)*m) +
                FLOP_ADD*(0.5*n*m*m-(1.0/6.0)*m*m*m-0.5*n*m+        (1.0/6.0)*m));
    }

    template<int BlkSize, typename DeviceSpaceType, typename VectorTagType, typename AlgoTagType>
    void LU(const int N) {
      //constexpr int N = 100;

      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      double flop = (N*VectorLength)*FlopCount(BlkSize,BlkSize);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      ///
      /// Reference version using MKL DGETRF
      ///
      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> aref;
      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> a("a", N*VectorLength, BlkSize, BlkSize);
        Kokkos::View<int**,Kokkos::LayoutRight,HostSpaceType> p("p", N*VectorLength, BlkSize);

        for (int k=0;k<N*VectorLength;++k) {
          const int
            k0 = k/VectorLength,
            k1 = k%VectorLength;
          for (int i=0;i<BlkSize;++i) {
            a(k, i,   i  ) = 100*(k0+1);
            if ((i+1) < BlkSize) {
              a(k, i,   i+1) = 1*(k1+1);
              a(k, i+1, i  ) = 1*(k1+1);
            }
          }
        }

        {
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Dynamic> > policy(0, N*VectorLength);
            Kokkos::parallel_for
              (policy,
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                auto pp = Kokkos::subview(p, k, Kokkos::ALL());
                
                LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
                               BlkSize, BlkSize,
                               (double*)aa.data(), aa.stride_0(),
                               (int*)pp.data());
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;
          std::cout << std::setw(10) << "MKL LU"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << std::endl;
        }

        aref = a;
      }

      ///
      /// Plain version (comparable to micro BLAS version)
      ///

      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType>
          a_host("a_host", N*VectorLength, BlkSize, BlkSize);
        
        for (int k=0;k<N*VectorLength;++k) {
          const int
            k0 = k/VectorLength,
            k1 = k%VectorLength;
          for (int i=0;i<BlkSize;++i) {
            a_host(k, i, i) = 100*(k0+1);
            if ((i+1) < BlkSize) {
              a_host(k, i,   i+1) = 1*(k1+1);
              a_host(k, i+1, i  ) = 1*(k1+1);
            }
          }
        }

        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);

        Kokkos::deep_copy(a, a_host);

        {
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Dynamic> > policy(0, N*VectorLength);
            Kokkos::parallel_for
              (policy,
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                
                Serial::LU<AlgoTagType>::invoke(aa);
              });

            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          Kokkos::deep_copy(a_host, a);
          double diff = 0;
          for (int i=0;i<aref.dimension(0);++i)
            for (int j=0;j<aref.dimension(1);++j)
              for (int k=0;k<aref.dimension(2);++k)
                diff += std::abs(aref(i,j,k) - a_host(i,j,k));

          std::cout << std::setw(10) << "Plain"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }

      ///
      /// SIMD with appropriate data layout
      ///

      {
        typedef Vector<VectorTagType> VectorType;
        Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
          a_host("a_host", N, BlkSize, BlkSize);

        for (int k0=0;k0<N;++k0)
          for (int k1=0;k1<VectorLength;++k1) {
            for (int i=0;i<BlkSize;++i) 
              for (int j=0;j<BlkSize;++j) 
                a_host(k0, i,j)[k1] = 0;
            for (int i=0;i<BlkSize;++i) {
              a_host(k0, i, i)[k1] = 100*(k0+1);
              if ((i+1) < BlkSize) {
                a_host(k0, i,   i+1)[k1] = 1*(k1+1);
                a_host(k0, i+1, i  )[k1] = 1*(k1+1);
              }
            }
          }
        
        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);

        Kokkos::deep_copy(a, a_host);

        {
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Dynamic> > policy(0, N);
            Kokkos::parallel_for
              (policy,
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                
                Serial::LU<AlgoTagType>::invoke(aa);
              });

            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          Kokkos::deep_copy(a_host, a);
          double diff = 0;
          for (int i=0;i<aref.dimension(0);++i)
            for (int j=0;j<aref.dimension(1);++j)
              for (int k=0;k<aref.dimension(2);++k) 
                diff += std::abs(aref(i,j,k) - a_host(i/VectorLength,j,k)[i%VectorLength]);
          std::cout << std::setw(10) << "SIMD"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }
    }
  }
}

using namespace KokkosKernels;

template<typename VectorType,
         typename AlgoTagType>
void run(const int N) {
  typedef Kokkos::OpenMP ExecSpace;

  std::cout << "ExecSpace::  ";
  if (std::is_same<ExecSpace,Kokkos::Serial>::value)
    std::cout << "Kokkos::Serial " << std::endl;
  else
    ExecSpace::print_configuration(std::cout, false);

  Test::LU< 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::LU< 9, ExecSpace,VectorType,AlgoTagType>(N);
  Test::LU<15, ExecSpace,VectorType,AlgoTagType>(N);
  Test::LU<20, ExecSpace,VectorType,AlgoTagType>(N);

  // Test::LU< 4, ExecSpace,VectorType,AlgoTagType>();
  // Test::LU< 8, ExecSpace,VectorType,AlgoTagType>();
  // Test::LU<16, ExecSpace,VectorType,AlgoTagType>();
  // Test::LU<20, ExecSpace,VectorType,AlgoTagType>();
}

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  const int N = 512;

#if defined(__AVX__) || defined(__AVX2__)
  // std::cout << "\n Testing SIMD4 and Algo::LU::Unblocked\n";
  // run<VectorTag<SIMD<double>,4>,Algo::LU::Unblocked>();

  // std::cout << "\n Testing AVX256 and Algo::LU::Unblocked\n";
  // run<VectorTag<AVX<double>,4>,Algo::LU::Unblocked>();

  // std::cout << "\n Testing SIMD4 and Algo::LU::Blocked\n";
  // run<VectorTag<SIMD<double>,4>,Algo::LU::Blocked>();

  std::cout << "\n Testing AVX256 and Algo::LU::Blocked\n";
  run<VectorTag<AVX<double>,4>,Algo::LU::Blocked>(N/4);
#endif 

#if defined(__AVX512F__)
  // std::cout << "\n Testing SIMD8 and Algo::LU::Unblocked\n";
  // run<VectorTag<SIMD<double>,8>,Algo::LU::Unblocked>();

  // std::cout << "\n Testing AVX512 and Algo::LU::Unblocked\n";
  // run<VectorTag<AVX<double>,8>,Algo::LU::Unblocked>();

  // std::cout << "\n Testing SIMD8 and Algo::LU::Blocked\n";
  // run<VectorTag<SIMD<double>,8>,Algo::LU::Blocked>();

  std::cout << "\n Testing AVX512 and Algo::LU::Blocked\n";
  run<VectorTag<AVX<double>,8>,Algo::LU::Blocked>(N/8);
#endif 


  Kokkos::finalize();

  return 0;
}


