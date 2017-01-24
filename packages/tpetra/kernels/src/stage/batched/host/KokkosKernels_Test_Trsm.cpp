/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>
#include "mkl.h"

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_Vector.hpp"

#include "KokkosKernels_Trsm_Serial_Decl.hpp"
#include "KokkosKernels_Trsm_Serial_Impl.hpp"

namespace KokkosKernels {

  namespace Test {

#define FLOP_MUL 1.0 
#define FLOP_ADD 1.0

    double FlopCountLower(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      return (FLOP_MUL*(0.5*m*n*(n+1.0)) +
              FLOP_ADD*(0.5*m*n*(n-1.0)));
    }

    double FlopCountUpper(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      return (FLOP_MUL*(0.5*m*n*(n+1.0)) +
              FLOP_ADD*(0.5*m*n*(n-1.0)));
    }

    template<int test, int BlkSize, typename DeviceSpaceType, typename VectorTagType, typename AlgoTagType>
    void Trsm(const int N) {
      switch (test) {
      case 0: std::cout << "TestID = Left,  Lower, NoTrans,    UnitDiag\n"; break;
      case 1: std::cout << "TestID = Left,  Lower, NoTrans, NonUnitDiag\n"; break;
      case 2: std::cout << "TestID = Right, Upper, NoTrans,    UnitDiag\n"; break;
      case 3: std::cout << "TestID = Right, Upper, NoTrans, NonUnitDiag\n"; break;
      }

      //constexpr int N = 100;

      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      // when m == n, lower upper does not matter (unit and nonunit)
      double flop = 0;
      switch (test) {
      case 0:
      case 1:
        flop = FlopCountLower(BlkSize,BlkSize);
        break;
      case 2:
      case 3:
        flop = FlopCountUpper(BlkSize,BlkSize);        
        break;
      }
      flop *= (N*VectorLength);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      ///
      /// Reference version using MKL DTRSM
      ///
      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> bref;
      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
        amat("amat", N*VectorLength, BlkSize, BlkSize),
        bmat("bmat", N*VectorLength, BlkSize, BlkSize);
      
      for (int k=0;k<N*VectorLength;++k) {
        const int
          k0 = k/VectorLength,
          k1 = k%VectorLength;
        for (int i=0;i<BlkSize;++i)
          for (int j=0;j<BlkSize;++j) {
            amat(k, i, j) = (rand()/(double)RAND_MAX)+4.0*(i==j);
            bmat(k, i, j) = (rand()/(double)RAND_MAX)+1.0;
          }
      }
        
      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N*VectorLength, BlkSize, BlkSize),
          b("b", N*VectorLength, BlkSize, BlkSize);
        
        for (int k=0;k<N*VectorLength;++k) {
          const int
            k0 = k/VectorLength,
            k1 = k%VectorLength;
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              a(k, i, j) = amat(k, i, j);
              b(k, i, j) = bmat(k, i, j);
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
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                
                switch (test) {
                case 0:
                  cblas_dtrsm(CblasRowMajor, 
                              CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                              BlkSize, BlkSize, 
                              1.0,
                              (double*)aa.data(), aa.stride_0(),
                              (double*)bb.data(), bb.stride_0());
                  break;
                case 1:
                  cblas_dtrsm(CblasRowMajor, 
                              CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                              BlkSize, BlkSize, 
                              1.0,
                              (double*)aa.data(), aa.stride_0(),
                              (double*)bb.data(), bb.stride_0());
                  break;
                case 2:
                  cblas_dtrsm(CblasRowMajor, 
                              CblasRight, CblasUpper, CblasNoTrans, CblasUnit,
                              BlkSize, BlkSize, 
                              1.0,
                              (double*)aa.data(), aa.stride_0(),
                              (double*)bb.data(), bb.stride_0());
                  break;
                case 3:
                  cblas_dtrsm(CblasRowMajor, 
                              CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
                              BlkSize, BlkSize, 
                              1.0,
                              (double*)aa.data(), aa.stride_0(),
                              (double*)bb.data(), bb.stride_0());
                  break;
                }
              });

            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          double sum = 0;
          for (int i=0;i<b.dimension(0);++i)
            for (int j=0;j<b.dimension(1);++j)
              for (int k=0;k<b.dimension(2);++k)
                sum += std::abs(b(i,j,k));
          
          std::cout << std::setw(10) << "MKL TRSM"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " sum abs(B)  = " << sum
                    << std::endl;

          bref = b;
        }
      }

      ///
      /// Plain version (comparable to micro BLAS version)
      ///

      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a_host("a_host", N*VectorLength, BlkSize, BlkSize),
          b_host("b_host", N*VectorLength, BlkSize, BlkSize);

        for (int k=0;k<N*VectorLength;++k) {
          const int
            k0 = k/VectorLength,
            k1 = k%VectorLength;
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              a_host(k, i, j) = amat(k, i, j);
              b_host(k, i, j) = bmat(k, i, j);
            }
        }

        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);

        Kokkos::deep_copy(a, a_host);
        Kokkos::deep_copy(b, b_host);

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
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                
                switch (test) {
                case 0: 
                  Serial::Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                case 1:
                  Serial::Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NonUnit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                case 2:
                  Serial::Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::Unit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                case 3:
                  Serial::Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                }
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;
          
          Kokkos::deep_copy(b_host, b);
          double diff = 0;
          for (int i=0;i<bref.dimension(0);++i)
            for (int j=0;j<bref.dimension(1);++j)
              for (int k=0;k<bref.dimension(2);++k)
                diff += std::abs(bref(i,j,k) - b_host(i,j,k));

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
          a_host("a_host", N, BlkSize, BlkSize),
          b_host("b_host", N, BlkSize, BlkSize);

        for (int k0=0;k0<N;++k0) 
          for (int k1=0;k1<VectorLength;++k1) 
            for (int i=0;i<BlkSize;++i)
              for (int j=0;j<BlkSize;++j) {
                a_host(k0, i, j)[k1] = amat(k0*VectorLength+k1, i, j);
                b_host(k0, i, j)[k1] = bmat(k0*VectorLength+k1, i, j);
              }
        
        auto a = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);
        
        Kokkos::deep_copy(a, a_host);
        Kokkos::deep_copy(b, b_host);
        
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
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                
                switch (test) {
                case 0:
                  Serial::Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                case 1:
                  Serial::Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NonUnit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                case 2:
                  Serial::Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::Unit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                case 3:
                  Serial::Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                }
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;
          
          Kokkos::deep_copy(b_host, b);
          double diff = 0;
          for (int i=0;i<bref.dimension(0);++i)
            for (int j=0;j<bref.dimension(1);++j)
              for (int k=0;k<bref.dimension(2);++k)
                diff += std::abs(bref(i,j,k) - b_host(i/VectorLength,j,k)[i%VectorLength]);

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

  /// UnitDiag ( division does not involve )
  ///
  // Test::Trsm<0, 5, ExecSpace,VectorType,AlgoTagType>();
  // Test::Trsm<0, 9, ExecSpace,VectorType,AlgoTagType>();
  // Test::Trsm<0,15, ExecSpace,VectorType,AlgoTagType>();
  // Test::Trsm<0,20, ExecSpace,VectorType,AlgoTagType>();

  // Test::Trsm<2, 5, ExecSpace,VectorType,AlgoTagType>();
  // Test::Trsm<2, 9, ExecSpace,VectorType,AlgoTagType>();
  // Test::Trsm<2,15, ExecSpace,VectorType,AlgoTagType>();
  // Test::Trsm<2,20, ExecSpace,VectorType,AlgoTagType>();

  /// NonUnitDiag ( division seems to be important )

  Test::Trsm<1, 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<1, 9, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<1,15, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<1,20, ExecSpace,VectorType,AlgoTagType>(N);

  // Test::Trsm<3, 5, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Trsm<3, 9, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Trsm<3,15, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Trsm<3,20, ExecSpace,VectorType,AlgoTagType>(N);
}

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  const int N = 512;

#if defined(__AVX__) || defined(__AVX2__)
  // std::cout << "\n Testing SIMD4 and Algo::Trsm::Unblocked\n";
  // run<VectorTag<SIMD<double>,4>,Algo::Trsm::Unblocked>(N/4);

  // std::cout << "\n Testing AVX256 and Algo::Trsm::Unblocked\n";
  // run<VectorTag<AVX<double>,4>,Algo::Trsm::Unblocked>(N/4);

  // std::cout << "\n Testing SIMD4 and Algo::Trsm::Blocked\n";
  // run<VectorTag<SIMD<double>,4>,Algo::Trsm::Blocked>(N/4);

  std::cout << "\n Testing AVX256 and Algo::Trsm::Blocked\n";
  run<VectorTag<AVX<double>,4>,Algo::Trsm::Blocked>(N/4);
#endif
  
#if defined(__AVX512F__)
  // std::cout << "\n Testing SIMD8 and Algo::Trsm::Unblocked\n";
  // run<VectorTag<SIMD<double>,8>,Algo::Trsm::Unblocked>(N/8);

  // std::cout << "\n Testing AVX512 and Algo::Trsm::Unblocked\n";
  // run<VectorTag<AVX<double>,8>,Algo::Trsm::Unblocked>(N/8);

  // std::cout << "\n Testing SIMD8 and Algo::Trsm::Blocked\n";
  // run<VectorTag<SIMD<double>,8>,Algo::Trsm::Blocked>(N/8);

  std::cout << "\n Testing AVX512 and Algo::Trsm::Blocked\n";
  run<VectorTag<AVX<double>,8>,Algo::Trsm::Blocked>(N/8);
#endif

  Kokkos::finalize();

  return 0;
}
