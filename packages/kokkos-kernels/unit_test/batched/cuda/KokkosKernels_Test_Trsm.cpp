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

#include "KokkosKernels_Trsm_Decl.hpp"
#include "KokkosKernels_Trsm_Serial_Impl.hpp"
//#include "KokkosKernels_Trsm_Team_Impl.hpp"

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

    template<int test, int BlkSize, int NumCols, typename DeviceSpaceType, typename VectorTagType, typename AlgoTagType>
    void Trsm(const int N) {
      typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;
      //typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;

      switch (test) {
      case 0: std::cout << "TestID = Left,  Lower, NoTrans,    UnitDiag\n"; break;
      case 1: std::cout << "TestID = Left,  Lower, NoTrans, NonUnitDiag\n"; break;
      case 2: std::cout << "TestID = Right, Upper, NoTrans,    UnitDiag\n"; break;
      case 3: std::cout << "TestID = Right, Upper, NoTrans, NonUnitDiag\n"; break;
      case 4: std::cout << "TestID = Left,  Upper, NoTrans, NonUnitDiag\n"; break;
      }

      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      // when m == n, lower upper does not matter (unit and nonunit)
      double flop = 0;
      switch (test) {
      case 0:
      case 1:
        flop = FlopCountLower(BlkSize,NumCols);
        break;
      case 2:
      case 3:
      case 4:
        flop = FlopCountUpper(BlkSize,NumCols);
        break;
      }
      flop *= (N*VectorLength);

      const double tmax = 1.0e15;

      typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      Kokkos::View<ValueType***,Kokkos::LayoutLeft,HostSpaceType>
        amat("amat", N*VectorLength, BlkSize, BlkSize),
        bmat("bmat", N*VectorLength, BlkSize, NumCols),
        bref("bmat", N*VectorLength, BlkSize, NumCols);

      {
        Random random;
        for (int k=0;k<N*VectorLength;++k) {
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j)
              amat(k, i, j) = random.value() + 4.0*(i==j);
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<NumCols;++j)
              bmat(k, i, j) = random.value();
        }
      }

      constexpr size_t LLC_CAPACITY = 56*64*4*1024*1024;
      Flush<LLC_CAPACITY,DeviceSpaceType> flush;

#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
      {
        ///                                                                                                           
        /// CUBLAS Batch version                                                                                    
        ///  

        const Kokkos::LayoutStride stride(N*VectorLength, BlkSize*BlkSize,
                                          BlkSize, 1,
                                          BlkSize, BlkSize);

        Kokkos::View<ValueType***,Kokkos::LayoutStride,DeviceSpaceType>
          a("a", stride),
          b("b", stride);       

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

        const double one(1.0); //, zero(0.0);
        {
          double tavg = 0, tmin = tmax;
          ValueType 
            *aa[N*VectorLength],
            *bb[N*VectorLength];
          for (int k=0;k<N*VectorLength;++k) {
            aa[k] = a.data() + k*a.stride_0();
            bb[k] = b.data() + k*b.stride_0();
          }
          ValueType 
            **aa_device,
            **bb_device;
          if (cudaMalloc(&aa_device, N*VectorLength*sizeof(ValueType*)) != cudaSuccess || 
              cudaMalloc(&bb_device, N*VectorLength*sizeof(ValueType*)) != cudaSuccess) {
            Kokkos::abort("CUDA memory allocation failed\n"); 
          }
          if (cudaMemcpy(aa_device, aa, sizeof(ValueType*)*N*VectorLength, cudaMemcpyHostToDevice) != cudaSuccess ||
              cudaMemcpy(bb_device, bb, sizeof(ValueType*)*N*VectorLength, cudaMemcpyHostToDevice) != cudaSuccess) {
            Kokkos::abort("CUDA memcpy failed\n");
          }
          DeviceSpaceType::fence();
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();

            // initialize matrices
            Kokkos::deep_copy(a, amat_device);
            Kokkos::deep_copy(b, bmat_device);
            
            timer.reset();
            switch (test) {
            case 0: {
              // Left,  Lower, NoTrans,    UnitDiag 
              stat = cublasDtrsmBatched(handle, 
                                        CUBLAS_SIDE_LEFT,
                                        CUBLAS_FILL_MODE_LOWER,
                                        CUBLAS_OP_N,
                                        CUBLAS_DIAG_UNIT, 
                                        BlkSize, NumCols,
                                        &one, 
                                        (const ValueType**)aa_device, BlkSize, 
                                        (ValueType**)bb_device, BlkSize, 
                                        N*VectorLength);
              break;
            }
            case 1: {
              // Left,  Lower, NoTrans, NonUnitDiag 
              stat = cublasDtrsmBatched(handle, 
                                        CUBLAS_SIDE_LEFT,
                                        CUBLAS_FILL_MODE_LOWER,
                                        CUBLAS_OP_N,
                                        CUBLAS_DIAG_NON_UNIT,
                                        BlkSize, NumCols,
                                        &one, 
                                        (const ValueType**)aa_device, BlkSize, 
                                        (ValueType**)bb_device, BlkSize, 
                                        N*VectorLength);
              break;
            }
            case 2: {
              // Right, Upper, NoTrans,    UnitDiag
              stat = cublasDtrsmBatched(handle, 
                                        CUBLAS_SIDE_RIGHT,
                                        CUBLAS_FILL_MODE_UPPER,
                                        CUBLAS_OP_N,
                                        CUBLAS_DIAG_UNIT,
                                        BlkSize, NumCols,
                                        &one, 
                                        (const ValueType**)aa_device, BlkSize, 
                                        (ValueType**)bb_device, BlkSize, 
                                        N*VectorLength);
              break;             
            }
            case 3: { 
              // Right, Upper, NoTrans, NonUnitDiag
              stat = cublasDtrsmBatched(handle, 
                                        CUBLAS_SIDE_RIGHT,
                                        CUBLAS_FILL_MODE_UPPER,
                                        CUBLAS_OP_N,
                                        CUBLAS_DIAG_NON_UNIT,
                                        BlkSize, NumCols,
                                        &one, 
                                        (const ValueType**)aa_device, BlkSize, 
                                        (ValueType**)bb_device, BlkSize, 
                                        N*VectorLength);
              break;             
            }
            case 4: {
              // Left,  Upper, NoTrans, NonUnitDiag
              stat = cublasDtrsmBatched(handle, 
                                        CUBLAS_SIDE_LEFT,
                                        CUBLAS_FILL_MODE_UPPER,
                                        CUBLAS_OP_N,
                                        CUBLAS_DIAG_NON_UNIT,
                                        BlkSize, NumCols,
                                        &one, 
                                        (const ValueType**)aa_device, BlkSize, 
                                        (ValueType**)bb_device, BlkSize, 
                                        N*VectorLength);
              break;                           
            }
            }

            if (stat != CUBLAS_STATUS_SUCCESS) {
              Kokkos::abort("CUBLAS Trsm Batched failed\n");
            }
            DeviceSpaceType::fence();
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          auto bsol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b);
          Kokkos::deep_copy(bsol, b);
          Kokkos::deep_copy(bref, bsol);

          if (cudaFree(aa_device) != cudaSuccess || 
              cudaFree(bb_device) != cudaSuccess) {
            Kokkos::abort("CUDA memory free failed\n"); 
          }
          
          std::cout << std::setw(8) << "CUBLAS"
                    << std::setw(8) << "Batch"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " NumCols = " << std::setw(3) << NumCols
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
        /// Plain version (comparable to micro BLAS version)
        ///
        
        Kokkos::View<ValueType***,Kokkos::LayoutLeft,DeviceSpaceType>
          a("a", N*VectorLength, BlkSize, BlkSize),
          b("b", N*VectorLength, BlkSize, NumCols);
        
        {
          double tavg = 0, tmin = tmax;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            flush.run();
            
            // initialize matrices
            Kokkos::deep_copy(a, amat);
            Kokkos::deep_copy(b, bmat);
            
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
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
                case 4:
                  Serial::Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,AlgoTagType>::
                    invoke(1.0, aa, bb);
                  break;
                }
              });

            DeviceSpaceType::fence();
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          auto bsol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b);
          Kokkos::deep_copy(bsol, b);
          
          double diff = 0;
          for (int i=0;i<bref.dimension(0);++i)
            for (int j=0;j<bref.dimension(1);++j)
              for (int k=0;k<bref.dimension(2);++k)
                diff += std::abs(bref(i,j,k) - bsol(i,j,k));

          std::cout << std::setw(8) << "Kokkos"
                    << std::setw(8) << "Range"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " NumCols = " << std::setw(3) << NumCols
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }
      std::cout << "\n\n";
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

  std::cout << "\n\n Used for Factorization \n\n";

  /// Left, Lower, NoTrans, UnitDiag (used in LU factorization and LU solve)

  Test::Trsm<0, 3, 3, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<0, 5, 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<0,10,10, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<0,15,15, ExecSpace,VectorType,AlgoTagType>(N);

  /// Left, Lower, NoTrans, NonUnitDiag

  Test::Trsm<1, 3, 3, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<1, 5, 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<1,10,10, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<1,15,15, ExecSpace,VectorType,AlgoTagType>(N);

  /// Right, Upper, NoTrans, UnitDiag

  Test::Trsm<2, 3, 3, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<2, 5, 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<2,10,10, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<2,15,15, ExecSpace,VectorType,AlgoTagType>(N);

  /// Right, Upper, NoTrans, NonUnitDiag (used in LU factorization)

  Test::Trsm<3, 3, 3, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<3, 5, 5, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<3,10,10, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<3,15,15, ExecSpace,VectorType,AlgoTagType>(N);

  std::cout << "\n\n Used for Solve \n\n";

  /// Left, Lower, NoTrans, UnitDiag (used in LU solve)

  Test::Trsm<0, 3, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<0, 5, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<0,10, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<0,15, 1, ExecSpace,VectorType,AlgoTagType>(N);

  /// Left, Upper, Notrans, NonUnitDiag (user in LU solve)

  Test::Trsm<4, 3, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<4, 5, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<4,10, 1, ExecSpace,VectorType,AlgoTagType>(N);
  Test::Trsm<4,15, 1, ExecSpace,VectorType,AlgoTagType>(N);
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

    std::cout << "\n Testing SIMT-" << VectorLength << " and Algo::Trsm::Unblocked\n";
    run<VectorTag<SIMT<double>,VectorLength>,Algo::Trsm::Unblocked>(N/VectorLength);

    std::cout << "\n Testing SIMT-" << VectorLength << " and Algo::Trsm::Blocked\n";
    run<VectorTag<SIMT<double>,VectorLength>,Algo::Trsm::Blocked>(N/VectorLength);
  }

  Kokkos::finalize();

  return 0;
}
