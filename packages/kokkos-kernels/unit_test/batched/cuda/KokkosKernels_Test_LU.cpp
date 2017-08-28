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

#include "KokkosKernels_LU_Decl.hpp"
#include "KokkosKernels_LU_Serial_Impl.hpp"
//#include "KokkosKernels_LU_Team_Impl.hpp"

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
      typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;
      //typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;
      
      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;
      
      const double flop = (N*VectorLength)*FlopCount(BlkSize,BlkSize);
      const double tmax = 1.0e15;

      typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      Kokkos::View<ValueType***,Kokkos::LayoutLeft,HostSpaceType>
        amat("amat", N*VectorLength, BlkSize, BlkSize),
        aref("aref", N*VectorLength, BlkSize, BlkSize);

      {
        Random random;
        for (int k=0;k<N*VectorLength;++k) {
          // use tridiagonal matrices; for now we just check elementwise l/u factors
          // do not allow pivots
          for (int i=0;i<BlkSize;++i) {
            amat(k, i, i) = random.value() + 10.0;
            if ((i+1) < BlkSize) {
              amat(k, i, i+1) = random.value() + 1.0;
              amat(k, i+1, i) = random.value() + 1.0;
            }
          }
          
          // ValueType d[BlkSize], v[BlkSize][BlkSize];
          // for (int i=0;i<BlkSize;++i) {
          //   d[i] = random.value() + 1.0; // positive value
          //   for (int j=0;j<BlkSize;++j)
          //     v[i][j] = random.value();
          // }
          // for (int i=0;i<BlkSize;++i)
          //   for (int j=0;j<BlkSize;++j)
          //     for (int l=0;l<BlkSize;++l)
          //       amat(k, i, j) = v[i][l]*d[l]*v[l][j];
        }
      }

      constexpr size_t LLC_CAPACITY = 56*64*4*1024*1024;
      Flush<LLC_CAPACITY> flush;
      
#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
      {
        ///
        /// CUBLAS Batch version
        ///
        const Kokkos::LayoutStride stride(N*VectorLength, BlkSize*BlkSize,
                                          BlkSize, 1,
                                          BlkSize, BlkSize);
        
        Kokkos::View<ValueType***,Kokkos::LayoutStride,DeviceSpaceType> a("a", stride);        
        Kokkos::View<int*,DeviceSpaceType>                           info("info", N*VectorLength);

        cublasStatus_t stat;
        cublasHandle_t handle;

        stat = cublasCreate(&handle);
        if (stat != CUBLAS_STATUS_SUCCESS)
          Kokkos::abort("CUBLAS initialization failed\n");

        auto amat_device = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), amat);
        Kokkos::deep_copy(amat_device, amat);

        DeviceSpaceType::fence();
        {
          double tavg = 0, tmin = tmax;
          ValueType *aa[N*VectorLength];

          for (int k=0;k<N*VectorLength;++k) {
            aa[k] = a.data() + k*a.stride_0();
          }
          ValueType **aa_device;
          if (cudaMalloc(&aa_device, N*VectorLength*sizeof(ValueType*)) != cudaSuccess) {
            Kokkos::abort("CUDA memory allocation failed\n");
          }
          if (cudaMemcpy(aa_device, aa, sizeof(ValueType*)*N*VectorLength, cudaMemcpyHostToDevice) != cudaSuccess) {
            Kokkos::abort("CUDA memcpy failed\n");
          }
          DeviceSpaceType::fence();          
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();
            
            // initialize matrix
            Kokkos::deep_copy(a, amat_device);

            DeviceSpaceType::fence();
            timer.reset();

            stat = cublasDgetrfBatched(handle, 
                                       BlkSize, 
                                       (ValueType**)aa_device, BlkSize, 
                                       NULL, 
                                       (int*)info.data(), 
                                       N*VectorLength);
            if (stat != CUBLAS_STATUS_SUCCESS) {
              Kokkos::abort("CUBLAS LU Batched failed\n");
            }
            
            DeviceSpaceType::fence();
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;
          
          auto asol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a);
          Kokkos::deep_copy(asol, a);
          Kokkos::deep_copy(aref, asol);
          
          if (cudaFree(aa_device) != cudaSuccess) {
            Kokkos::abort("CUDA memory free failed\n");
          }
          
          std::cout << std::setw(8) << "CUBLAS"
                    << std::setw(8) << "Batch"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
                    << std::endl;
        }
      }
#endif

      {
        ///
        /// Plain version (comparable to micro BLAS version)
        ///
        
        Kokkos::View<ValueType***,Kokkos::LayoutLeft,DeviceSpaceType>
          a("a", N*VectorLength, BlkSize, BlkSize);

        {
          double tavg = 0, tmin = tmax;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();

            // initialize matrix
            Kokkos::deep_copy(a, amat);

            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
            Kokkos::parallel_for
              (policy,
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());

                Serial::LU<AlgoTagType>::invoke(aa);
              });

            DeviceSpaceType::fence();
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;
          
          auto asol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a);
          Kokkos::deep_copy(asol, a);
          
          double diff = 0;
          for (int i=0;i<aref.dimension(0);++i)
            for (int j=0;j<aref.dimension(1);++j)
              for (int k=0;k<aref.dimension(2);++k)
                diff += std::abs(aref(i,j,k) - asol(i,j,k));

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

  Test::LU< 3, ExecSpace,VectorType,AlgoTagType>(N);
  Test::LU< 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::LU<10, ExecSpace,VectorType,AlgoTagType>(N);
  Test::LU<15, ExecSpace,VectorType,AlgoTagType>(N);
}

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  int N = 128*128;

  for (int i=1;i<argc;++i) {
    const std::string& token = argv[i];
    if (token == std::string("-N")) N = std::atoi(argv[++i]);
  }

  constexpr int VectorLength = 4;

  {
    std::cout << " N = " << N << std::endl;
    
    std::cout << "\n Testing SIMT-" << VectorLength << " and Algo::LU::Unblocked\n";
    run<VectorTag<SIMT<double>,VectorLength>,Algo::LU::Unblocked>(N/VectorLength);
    
    std::cout << "\n Testing SIMT-" << VectorLength << " and Algo::LU::Blocked\n";
    run<VectorTag<SIMT<double>,VectorLength>,Algo::LU::Blocked>(N/VectorLength);
  }

  Kokkos::finalize();

  return 0;
}


