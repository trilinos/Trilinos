/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>

#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cublas_api.h"
#endif

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_Vector.hpp"

#include "KokkosKernels_Gemm_Decl.hpp"
#include "KokkosKernels_Gemm_Serial_Impl.hpp"
//#include "KokkosKernels_Gemm_Team_Impl.hpp"

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

      const double flop = (N*VectorLength)*FlopCount(BlkSize,BlkSize,BlkSize);
      const double tmax = 1.0e15;

      typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      Kokkos::View<ValueType***,Kokkos::LayoutLeft,HostSpaceType>
        amat("amat", N*VectorLength, BlkSize, BlkSize),
        bmat("bmat", N*VectorLength, BlkSize, BlkSize),
        cref("cref", N*VectorLength, BlkSize, BlkSize);

      {
        Random random;
        for (int k=0;k<N*VectorLength;++k)
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              amat(k, i, j) = random.value();
              bmat(k, i, j) = random.value();
            }
      }

      // P100 L2 cache 4MB per core
      constexpr size_t LLC_CAPACITY = 56*64*4*1024*1024;
      Flush<LLC_CAPACITY,DeviceSpaceType> flush;

#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
      {
        ///
        /// CUBLAS Strided version
        ///

        const Kokkos::LayoutStride stride(N*VectorLength, BlkSize*BlkSize,
                                          BlkSize, 1,
                                          BlkSize, BlkSize);

        Kokkos::View<ValueType***,Kokkos::LayoutStride,DeviceSpaceType>
          a("a", stride),
          b("b", stride),
          c("c", stride);

        double tavg = 0, tmin = tmax;

        cublasStatus_t stat;
        cublasHandle_t handle;

        stat = cublasCreate(&handle);
        if (stat != CUBLAS_STATUS_SUCCESS)
          Kokkos::abort("CUBLAS initialization failed\n");

        auto amat_device = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), amat);
        auto bmat_device = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), bmat);

        Kokkos::deep_copy(amat_device, amat);
        Kokkos::deep_copy(bmat_device, bmat);

        DeviceSpaceType::fence();

        const double one(1.0), zero(0.0);
        {
          tavg = 0; tmin = tmax;

          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();

            // initialize matrices
            Kokkos::deep_copy(a, amat_device);
            Kokkos::deep_copy(b, bmat_device);
            Kokkos::deep_copy(c, 0);

            DeviceSpaceType::fence();
            timer.reset();

            stat = cublasDgemmStridedBatched(handle,
                                             CUBLAS_OP_N,
                                             CUBLAS_OP_N,
                                             BlkSize, BlkSize, BlkSize,
                                             &one,
                                             (const ValueType*)a.data(), BlkSize, BlkSize*BlkSize,
                                             (const ValueType*)b.data(), BlkSize, BlkSize*BlkSize,
                                             &zero,
                                             (ValueType*)c.data(), BlkSize, BlkSize*BlkSize,
                                             N*VectorLength);

            DeviceSpaceType::fence();
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
          Kokkos::deep_copy(csol, c);
          Kokkos::deep_copy(cref, csol);

          std::cout << std::setw(8) << "CUBLAS"
                    << std::setw(8) << "Strided"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
                    << std::endl;
        }
        cublasDestroy(handle);
      }
#endif

      {
        ///
        /// Plain version
        ///
        Kokkos::View<ValueType***,Kokkos::LayoutLeft,DeviceSpaceType>
          a("a", N*VectorLength, BlkSize, BlkSize),
          b("b", N*VectorLength, BlkSize, BlkSize),
          c("c", N*VectorLength, BlkSize, BlkSize);

        double tavg = 0, tmin = tmax;
        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);

          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();

            // initialize matrices
            Kokkos::deep_copy(a, amat);
            Kokkos::deep_copy(b, bmat);
            Kokkos::deep_copy(c, 0);

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
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
          Kokkos::deep_copy(csol, c);

          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - csol(i,j,k));

          std::cout << std::setw(8) << "Kokkos"
                    << std::setw(8) << "Range"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
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

  // Test::Gemm< 4, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemm< 8, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemm<16, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemm<20, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemm<32, ExecSpace,VectorType,AlgoTagType>(N);
  // Test::Gemm<64, ExecSpace,VectorType,AlgoTagType>(N);

  Test::Gemm< 3, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemm< 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemm<10, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Gemm<15, ExecSpace,VectorType,AlgoTagType>(N);
}

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  const int ntest = 1;
  //const int N[6] = { 256, 512, 768, 1024, 1280, 1536 };
  int N[1] = { 128*128 };

  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (token == std::string("-N")) N[0] = std::atoi(argv[++i]);
  }

  constexpr int VectorLength = 4;

  {
    for (int i=0;i<ntest;++i) {
      std::cout << " N = " << N[i] << std::endl;

      std::cout << "\n Testing SIMT-" << VectorLength << " and Algo::Gemm::Unblocked\n";
      run<VectorTag<SIMT<double>,VectorLength>,Algo::Gemm::Unblocked>(N[i]/VectorLength);

      std::cout << "\n Testing SIMT-" << VectorLength << " and Algo::Gemm::Blocked\n";
      run<VectorTag<SIMT<double>,VectorLength>,Algo::Gemm::Blocked>(N[i]/VectorLength);
    }
  }

  Kokkos::finalize();

  return 0;
}
