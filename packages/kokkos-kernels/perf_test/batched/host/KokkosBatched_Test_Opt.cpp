/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "mkl.h"

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"

#include "KokkosBatched_Trsm_Serial_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"

#include "KokkosBatched_LU_Serial_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"

namespace KokkosKernels {

  namespace Test {

#define FLOP_MUL 1.0 // 6.0
#define FLOP_ADD 1.0 // 2.0

    double FlopCountLowerTrsm(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      return (FLOP_MUL*(0.5*m*n*(n+1.0)) +
              FLOP_ADD*(0.5*m*n*(n-1.0)));
    }

    double FlopCountUpperTrsm(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      return (FLOP_MUL*(0.5*m*n*(n+1.0)) +
              FLOP_ADD*(0.5*m*n*(n-1.0)));
    }

    double FlopCountGemm(int mm, int nn, int kk) {
      double m = (double)mm;    double n = (double)nn;    double k = (double)kk;
      return (FLOP_MUL*(m*n*k) +
              FLOP_ADD*(m*n*k));
    }

    double FlopCountLU(int mm, int nn) {
      double m = (double)mm;    double n = (double)nn;
      if (m > n)
        return (FLOP_MUL*(0.5*m*n*n-(1.0/6.0)*n*n*n+0.5*m*n-0.5*n*n+(2.0/3.0)*n) +
                FLOP_ADD*(0.5*m*n*n-(1.0/6.0)*n*n*n-0.5*m*n+        (1.0/6.0)*n));
      else
        return (FLOP_MUL*(0.5*n*m*m-(1.0/6.0)*m*m*m+0.5*n*m-0.5*m*m+(2.0/3.0)*m) +
                FLOP_ADD*(0.5*n*m*m-(1.0/6.0)*m*m*m-0.5*n*m+        (1.0/6.0)*m));
    }

    template<int BlkSize, 
             typename DeviceSpaceType, 
             typename VectorTagType, 
             typename AlgoLUTagType = void,
             typename AlgoTrsmTagType = void,
             typename AlgoGemmTagType = void>
    void Opt() {
      constexpr int N = 100;

      typedef typename VectorTagType::value_type ValueType;
      constexpr int VectorLength = VectorTagType::length;

      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      typedef Vector<VectorTagType> VectorType;
      Kokkos::View<VectorType***,Kokkos::LayoutRight,DeviceSpaceType> 
        ATL("ATL", N, BlkSize, BlkSize),
        ATR("ATR", N, BlkSize, BlkSize),
        ABL("ABL", N, BlkSize, BlkSize),
        ABR("ABR", N, BlkSize, BlkSize);

      std::srand(0);
      auto RandVal = [=](const ValueType low, const ValueType hi) {
        const ValueType val = (ValueType)std::rand()/RAND_MAX;
        return low + val*(hi - low); 
      };
      
      for (int k=0;k<N;++k)
        for (int i=0;i<BlkSize;++i) 
          for (int j=0;j<BlkSize;++j) {
            ATL(k, i, j) = VectorType(RandVal(-1,1) + (i==j)*100);
            ATR(k, i, j) = VectorType(RandVal(-1,1));
            ABL(k, i, j) = VectorType(RandVal(-1,1));
            ABR(k, i, j) = VectorType(RandVal(-1,1));
          }

      double t = 0;
      for (int iter=iter_begin;iter<iter_end;++iter) {
        DeviceSpaceType::fence();
        timer.reset();

        for (int k=0;k<N;++k) {
          
          if (!std::is_same<AlgoLUTagType,void>::value) {
            auto atl = Kokkos::subview(ATL, k, Kokkos::ALL(), Kokkos::ALL());
            Serial::LU<AlgoLUTagType>::invoke(atl);
          }

          if (!std::is_same<AlgoTrsmTagType,void>::value) {
            auto atl = Kokkos::subview(ATL, k, Kokkos::ALL(), Kokkos::ALL());

            auto atr = Kokkos::subview(ATR, k, Kokkos::ALL(), Kokkos::ALL());            
            Serial::Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NonUnit,AlgoTrsmTagType>::
              invoke(1.0, atl, atr);

            // auto abl = Kokkos::subview(ABL, k, Kokkos::ALL(), Kokkos::ALL());
            // Serial::Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,AlgoTrsmTagType>::
            //   invoke(1.0, atl, abl);
          }

          if (!std::is_same<AlgoGemmTagType,void>::value) {
            auto atr = Kokkos::subview(ATR, k, Kokkos::ALL(), Kokkos::ALL());   
            auto abl = Kokkos::subview(ABL, k, Kokkos::ALL(), Kokkos::ALL());                     
            auto abr = Kokkos::subview(ABR, k, Kokkos::ALL(), Kokkos::ALL());                     

            Serial::Gemm<Trans::NoTranspose,Trans::NoTranspose,AlgoGemmTagType>::
              invoke(1.0, atr, abl, 1.0, abr);
          }
        }
        DeviceSpaceType::fence();
        t += (iter >= 0)*timer.seconds();
      }
      t /= iter_end;
      
      
      const double flop = (N*VectorLength); //*
        (!std::is_same<AlgoLUTagType,void>::value)*FlopCountLU(BlkSize,BlkSize) + 
        (!std::is_same<AlgoTrsmTagType,void>::value)*FlopCountLowerTrsm(BlkSize,BlkSize) + 
        //(!std::is_same<AlgoTrsmTagType,void>::value)*FlopCountUpperTrsm(BlkSize,BlkSize) + 
        (!std::is_same<AlgoGemmTagType,void>::value)*FlopCountGemm(BlkSize,BlkSize,BlkSize); 

      std::cout << " BlkSize = " << std::setw(3) << BlkSize
                << " time = " << std::scientific << t
                << " flop/s = " << (flop/t)
                << std::endl;
    }
  }
}

using namespace KokkosKernels;

template<typename VectorType,
         typename AlgoLUTagType,
         typename AlgoTrsmTagType,
         typename AlgoGemmTagType>
void run() {
  Test::Opt< 5, Kokkos::Serial, VectorType, AlgoLUTagType, AlgoTrsmTagType, AlgoGemmTagType>();
}

int main(int argc, char *argv[]) {

  Kokkos::initialize();

#if 1
#if defined(__AVX__) || defined(__AVX2__) 
  typedef VectorTag<AVX<double>,4> VectorTagType;

  // typedef Algo::LU::Unblocked AlgoLUType;
  // typedef Algo::Trsm::Unblocked AlgoTrsmType;
  // typedef Algo::Gemm::Unblocked AlgoGemmType;

  typedef Algo::LU::Blocked AlgoLUType;
  typedef Algo::Trsm::Blocked AlgoTrsmType;
  typedef Algo::Gemm::Blocked AlgoGemmType;
#endif
#else
#if defined(__AVX512F__)
  typedef VectorTag<AVX<double>,8> VectorTagType;

  typedef Algo::LU::Unblocked AlgoLUType;
  typedef Algo::Trsm::Unblocked AlgoTrsmType;
  typedef Algo::Gemm::Unblocked AlgoGemmType;

  //typedef Algo::LU::Blocked AlgoLUType;
  //typedef Algo::Trsm::Blocked AlgoTrsmType;
  //typedef Algo::Gemm::Blocked AlgoGemmType;
#endif
#endif

  run<VectorTagType,void,void,AlgoGemmType>();
  //run<VectorTagType,void,AlgoTrsmType,void>();
  //run<VectorTagType,AlgoLUType,void,void>();


  Kokkos::finalize();

  return 0;
}


