/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>

#if defined(__KOKKOSBATCHED_INTEL_MKL__)
#include "mkl.h"
#endif

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Gemv_Decl.hpp"
#include "KokkosBatched_Gemv_Serial_Impl.hpp"

namespace KokkosKernels {

  namespace Test {

#define FLOP_MUL 1.0
#define FLOP_ADD 1.0

    double FlopCount(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;    
      return (FLOP_MUL*(m*n) +
              FLOP_ADD*(m*n));
    }
    
    template<int BlkSize, int NumVecs, typename DeviceSpaceType, typename VectorTagType, typename AlgoTagType>
    void Gemv(const int N) {
      typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;
      //typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;
      
      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      double flop = (N*VectorLength)*FlopCount(BlkSize,BlkSize)*NumVecs;

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;
      
      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> yref;
      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
        amat("amat", N*VectorLength, BlkSize, BlkSize);
      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
        xvec("xvec", N*VectorLength, NumVecs, BlkSize);

      {
        Random random;
        for (int k=0;k<N*VectorLength;++k) 
          for (int i=0;i<BlkSize;++i) {
            for (int j=0;j<NumVecs;++j)
              xvec(k, j, i) = random.value();
            for (int j=0;j<BlkSize;++j) 
              amat(k, i, j) = random.value();
          }
      }
      
      // for KNL
      constexpr size_t LLC_CAPACITY = 34*1024*1024;
      Flush<LLC_CAPACITY> flush;
      
      ///
      /// Reference version using MKL DGEMM
      ///
#if defined(__KOKKOSBATCHED_INTEL_MKL__)
      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N*VectorLength, BlkSize, BlkSize),
          x("x", N*VectorLength, NumVecs, BlkSize),
          y("y", N*VectorLength, NumVecs, BlkSize);

        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
          
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();
            
            // initialize matrices
            Kokkos::deep_copy(a, amat);
            Kokkos::deep_copy(x, xvec);
            Kokkos::deep_copy(y, 0);
            
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                for (int j=0;j<NumVecs;++j) {
                  auto xx = Kokkos::subview(x, k, j, Kokkos::ALL());
                  auto yy = Kokkos::subview(y, k, j, Kokkos::ALL());
                  
                  cblas_dgemv(CblasRowMajor, CblasNoTrans, 
                              BlkSize, BlkSize, 
                              1.0,
                              (double*)aa.data(), aa.stride_0(),
                              (double*)xx.data(), xx.stride_0(),
                              1.0,
                              (double*)yy.data(), yy.stride_0());
                }
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          std::cout << std::setw(12) << "MKL DGEMV"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " NumVecs = " << std::setw(3) << NumVecs
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << std::endl;

          yref = y;
        }
      }
#endif

      ///
      /// Plain version (comparable to micro BLAS version)
      ///
      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N*VectorLength, BlkSize, BlkSize),
          x("x", N*VectorLength, NumVecs, BlkSize),
          y("y", N*VectorLength, NumVecs, BlkSize);

        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
          
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();

            // initialize matrices
            Kokkos::deep_copy(a, amat);
            Kokkos::deep_copy(x, xvec);
            Kokkos::deep_copy(y, 0);

            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());

                for (int j=0;j<NumVecs;++j) {                
                  auto xx = Kokkos::subview(x, k, j, Kokkos::ALL());
                  auto yy = Kokkos::subview(y, k, j, Kokkos::ALL());
                
                  KokkosKernels::Serial::
                    Gemv<Trans::NoTranspose,AlgoTagType>::
                    invoke(1.0, aa, xx, 1.0, yy);
                }
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;

          double diff = 0;
          for (int i=0;i<yref.dimension(0);++i)
            for (int j=0;j<yref.dimension(1);++j)
              for (int k=0;k<yref.dimension(2);++k)
                diff += std::abs(yref(i,j,k) - y(i,j,k));

          std::cout << std::setw(12) << "Plain"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " NumVecs = " << std::setw(3) << NumVecs
                    << " time = " << std::scientific << t
                    << " flop/s = " << (flop/t)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }
      
      typedef Vector<VectorTagType> VectorType;
      Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
        amat_simd("amat_simd", N, BlkSize, BlkSize),
        xvec_simd("xvec_simd", N, NumVecs, BlkSize);
        
      for (int k0=0;k0<N;++k0)
        for (int k1=0;k1<VectorLength;++k1) 
          for (int i=0;i<BlkSize;++i) {
            for (int j=0;j<NumVecs;++j) 
              xvec_simd(k0, j, i)[k1] = xvec(k0*VectorLength+k1, j, i);
            for (int j=0;j<BlkSize;++j) 
              amat_simd(k0, i, j)[k1] = amat(k0*VectorLength+k1, i, j);
          }
      
      
      ///
      /// Serial SIMD with appropriate data layout
      ///
      {
        typedef Vector<VectorTagType> VectorType;
        Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N, BlkSize, BlkSize),
          x("x", N, NumVecs, BlkSize),
          y("y", N, NumVecs, BlkSize);
        
        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N);
          
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();

            // initialize matrices
            Kokkos::deep_copy(a, amat_simd);
            Kokkos::deep_copy(x, xvec_simd);
            Kokkos::deep_copy(y, 0);

            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());

                for (int j=0;j<NumVecs;++j) {                
                  auto xx = Kokkos::subview(x, k, j, Kokkos::ALL());
                  auto yy = Kokkos::subview(y, k, j, Kokkos::ALL());
                
                  KokkosKernels::Serial::
                    Gemv<Trans::NoTranspose,AlgoTagType>::
                    invoke(1.0, aa, xx, 1.0, yy);
                }
              });
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          t /= iter_end;
          
          double diff = 0;
          for (int i=0;i<yref.dimension(0);++i)
            for (int j=0;j<yref.dimension(1);++j)
              for (int k=0;k<yref.dimension(2);++k)
                diff += std::abs(yref(i,j,k) - y(i/VectorLength,j,k)[i%VectorLength]);

          std::cout << std::setw(12) << "Serial SIMD"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " NumVecs = " << std::setw(3) << NumVecs
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
  typedef typename VectorType::exec_space ExecSpace;

  std::cout << "ExecSpace::  "; 
  if (std::is_same<ExecSpace,Kokkos::Serial>::value) 
    std::cout << "Kokkos::Serial " << std::endl;
  else 
    ExecSpace::print_configuration(std::cout, false);

  // Test::Gemv< 4, 1, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemv< 8, 1, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemv<16, 1, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemv<20, 1, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemv<32, 1, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemv<64, 1, ExecSpace,VectorType,AlgoTagType>(N);

  Test::Gemv< 3, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemv< 5, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemv<10, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemv<15, 1, ExecSpace,VectorType,AlgoTagType>(N);
}

int main(int argc, char *argv[]) {
  
  Kokkos::initialize(argc, argv);

  const int ntest = 1;
  //const int N[6] = { 256, 512, 768, 1024, 1280, 1536 };
  const int N[1] = { 128*128 };

#if defined(__AVX512F__)
  constexpr int VectorLength = 8;
#elif defined(__AVX2__) || defined(__AVX__)
  constexpr int VectorLength = 4;
#else
  static_assert(false, "AVX is not supported");
#endif

  {        
    for (int i=0;i<ntest;++i) {
      std::cout << " N = " << N[i] << std::endl;
      
      // std::cout << "\n Testing SIMD-" << VectorLength << " and Algo::Gemv::Unblocked\n";
      // run<VectorTag<SIMD<double>,VectorLength>,Algo::Gemv::Unblocked>(N[i]/VectorLength);
      
      std::cout << "\n Testing AVX-" << VectorLength << " and Algo::Gemv::Unblocked\n";
      run<VectorTag<AVX<double>,VectorLength>,Algo::Gemv::Unblocked>(N[i]/VectorLength);

      // std::cout << "\n Testing SIMD-" << VectorLength << " and Algo::Gemv::Blocked\n";
      // run<VectorTag<SIMD<double>,VectorLength>,Algo::Gemv::Blocked>(N[i]/VectorLength);
      
      std::cout << "\n Testing AVX-" << VectorLength << " and Algo::Gemv::Blocked\n";
      run<VectorTag<AVX<double>,VectorLength>,Algo::Gemv::Blocked>(N[i]/VectorLength);

      // std::cout << "\n Testing AVX-" << VectorLength << " and Algo::Gemv::CompactMKL\n";
      // run<VectorTag<AVX<double>,VectorLength>,Algo::Gemv::CompactMKL>(N[i]/VectorLength);
    }
  }

  Kokkos::finalize();
  
  return 0;
}
