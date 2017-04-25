/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>
#if defined(__KOKKOSKERNELS_LIBXSMM__)
#include "libxsmm.h"
#endif

#if defined(__KOKKOSKERNELS_INTEL_MKL__)
#include "mkl.h"
#endif

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_Vector.hpp"

#include "KokkosKernels_Gemm_Decl.hpp"
#include "KokkosKernels_Gemm_Serial_Impl.hpp"
#include "KokkosKernels_Gemm_Team_Impl.hpp"

bool hot = false;

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

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> cref;
      Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
        amat("amat", N*VectorLength, BlkSize, BlkSize),
        bmat("bmat", N*VectorLength, BlkSize, BlkSize);

      {
        Random random;
        for (int k=0;k<N*VectorLength;++k) 
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              amat(k, i, j) = random.value();
              bmat(k, i, j) = random.value();
            }
      }
      
      typedef Vector<VectorTagType> VectorType;
      Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
        amat_simd("amat_simd", N, BlkSize, BlkSize),
        bmat_simd("bmat_simd", N, BlkSize, BlkSize);
        
      for (int k0=0;k0<N;++k0)
        for (int k1=0;k1<VectorLength;++k1) 
          for (int i=0;i<BlkSize;++i)
            for (int j=0;j<BlkSize;++j) {
              amat_simd(k0, i, j)[k1] = amat(k0*VectorLength+k1, i, j);
              bmat_simd(k0, i, j)[k1] = bmat(k0*VectorLength+k1, i, j);
            }
      
      // for KNL
      constexpr size_t LLC_CAPACITY = 34*1024*1024;
      Flush<LLC_CAPACITY> flush;
      
      ///
      /// Reference version using MKL DGEMM
      ///
#if defined(__KOKKOSKERNELS_INTEL_MKL__)
      {
        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N*VectorLength, BlkSize, BlkSize),
          b("b", N*VectorLength, BlkSize, BlkSize),
          c("c", N*VectorLength, BlkSize, BlkSize);

        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
          
          double tavg = 0, tmin = tmax;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            if (!hot)
              flush.run();

            // initialize matrices
            if (!hot && iter == iter_begin) {
              Kokkos::deep_copy(a, amat);
              Kokkos::deep_copy(b, bmat);
            }
            Kokkos::deep_copy(c, 0);

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
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          std::cout << std::setw(12) << "MKL DGEMM"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
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
          double tavg = 0, tmin = tmax;

          MKL_INT blksize[1] = { BlkSize };
          MKL_INT lda[1] = { a.stride_1() };
          MKL_INT ldb[1] = { b.stride_1() };
          MKL_INT ldc[1] = { c.stride_1() };

          CBLAS_TRANSPOSE transA[1] = { CblasNoTrans };
          CBLAS_TRANSPOSE transB[1] = { CblasNoTrans };

          double one[1] = { 1.0 };
          MKL_INT size_per_grp[1] = { N*VectorLength };

          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            if (!hot)
              flush.run();

            // initialize matrices
            if (!hot && iter == iter_begin) {
              Kokkos::deep_copy(a, amat);
              Kokkos::deep_copy(b, bmat);
            }
            Kokkos::deep_copy(c, 0);

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
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c(i,j,k));

          std::cout << std::setw(12) << "MKL Batch"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }
#endif
#if defined(__KOKKOSKERNELS_INTEL_MKL_COMPACT_BATCHED__)
      {
        Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N, BlkSize, BlkSize),
          b("b", N, BlkSize, BlkSize),
          c("c", N, BlkSize, BlkSize);

        {
          double tavg = 0, tmin = tmax;

          MKL_INT blksize[1] = { BlkSize };
          MKL_INT lda[1] = { a.stride_1() };
          MKL_INT ldb[1] = { b.stride_1() };
          MKL_INT ldc[1] = { c.stride_1() };

          CBLAS_TRANSPOSE transA[1] = { CblasNoTrans };
          CBLAS_TRANSPOSE transB[1] = { CblasNoTrans };

          double one[1] = { 1.0 };
          MKL_INT size_per_grp[1] = { N*VectorLength };

          compact_t A_p, B_p, C_p;
          A_p.layout = CblasRowMajor;
          A_p.rows = blksize;
          A_p.cols = blksize;
          A_p.stride = lda;
          A_p.group_count = 1;
          A_p.size_per_group = size_per_grp;
          A_p.format = VectorLength;
          A_p.mat = (double*)a.data();

          B_p.layout = CblasRowMajor;
          B_p.rows = blksize;
          B_p.cols = blksize;
          B_p.stride = ldb;
          B_p.group_count = 1;
          B_p.size_per_group = size_per_grp;
          B_p.format = VectorLength;
          B_p.mat = (double*)b.data();

          C_p.layout = CblasRowMajor;
          C_p.rows = blksize;
          C_p.cols = blksize;
          C_p.stride = ldc;
          C_p.group_count = 1;
          C_p.size_per_group = size_per_grp;
          C_p.format = VectorLength;
          C_p.mat = (double*)c.data();

          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            if (!hot)
              flush.run();

            // initialize matrices
            if (!hot && iter == iter_begin) {
              Kokkos::deep_copy(a, amat_simd);
              Kokkos::deep_copy(b, bmat_simd);
            }
            Kokkos::deep_copy(c, 0);

            DeviceSpaceType::fence();
            timer.reset();

            cblas_dgemm_compute_batch(transA, 
                                      transB, 
                                      one, 
                                      &A_p, 
                                      &B_p, 
                                      one, 
                                      &C_p);

            DeviceSpaceType::fence();
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c(i/VectorLength,j,k)[i%VectorLength]);
          
          std::cout << std::setw(12) << "MKL Cmpct"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }
#endif
#endif

#if defined(__KOKKOSKERNELS_LIBXSMM__)
      {
        libxsmm_init();

        Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N*VectorLength, BlkSize, BlkSize),
          b("b", N*VectorLength, BlkSize, BlkSize),
          c("c", N*VectorLength, BlkSize, BlkSize);

        libxsmm_blasint 
          lda = a.stride_1(),
          ldb = b.stride_1(),
          ldc = c.stride_1();

        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
          
          double tavg = 0, tmin = tmax;

          // adjust column major order in xsmm
          char transA = 'N',  transB = 'N';
          libxsmm_blasint blksize = BlkSize;
          double one = 1.0;

          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            if (!hot)
              flush.run();

            // initialize matrices
            if (!hot && iter == iter_begin) {
              Kokkos::deep_copy(a, amat);
              Kokkos::deep_copy(b, bmat);
            }
            Kokkos::deep_copy(c, 0);

            DeviceSpaceType::fence();
            timer.reset();

            Kokkos::parallel_for
              (policy, 
               KOKKOS_LAMBDA(const int k) {
                auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
                auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
                auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());
                
                // column major
                libxsmm_gemm((const char*)&transA, 
                             (const char*)&transB, 
                             blksize, blksize, blksize,
                             (const double*)&one, 
                             (const double*)bb.data(), (const libxsmm_blasint*)&ldb,                                                                
                             (const double*)aa.data(), (const libxsmm_blasint*)&lda, 
                             (const double*)&one, 
                             (double*)cc.data(), (const libxsmm_blasint*)&ldc);
              });
            
            DeviceSpaceType::fence();
            const double t = timer.seconds();
            tmin = std::min(tmin, t);
            tavg += (iter >= 0)*t;
          }
          tavg /= iter_end;

          // adjust transpose
          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c(i,j,k));

          std::cout << std::setw(12) << "libxsmm"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
                    << " diff to ref = " << diff
                    << std::endl;
        }
        libxsmm_finalize();
      }
#endif
      // ///
      // /// Plain version (comparable to micro BLAS version)
      // ///
      // if (!std::is_same<AlgoTagType,Algo::Gemm::CompactMKL>::value) {
      //   Kokkos::View<ValueType***,Kokkos::LayoutRight,HostSpaceType> 
      //     a("a", N*VectorLength, BlkSize, BlkSize),
      //     b("b", N*VectorLength, BlkSize, BlkSize),
      //     c("c", N*VectorLength, BlkSize, BlkSize);

      //   {
      //     const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N*VectorLength);
          
      //     double tavg = 0, tmin = tmax;

      //     for (int iter=iter_begin;iter<iter_end;++iter) {
      //       // flush
      //       flush.run();

      //       // initialize matrices
      //       Kokkos::deep_copy(a, amat);
      //       Kokkos::deep_copy(b, bmat);
      //       Kokkos::deep_copy(c, 0);

      //       DeviceSpaceType::fence();
      //       timer.reset();

      //       Kokkos::parallel_for
      //         (policy, 
      //          KOKKOS_LAMBDA(const int k) {
      //           auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
      //           auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
      //           auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());
                
      //           KokkosKernels::Serial::
      //             Gemm<Trans::NoTranspose,Trans::NoTranspose,AlgoTagType>::
      //             invoke(1.0, aa, bb, 1.0, cc);
      //         });
            
      //       DeviceSpaceType::fence();
      //       const double t = timer.seconds();
      //       tmin = std::min(tmin, t);
      //       tavg += (iter >= 0)*t;
      //     }
      //     tavg /= iter_end;

      //     double diff = 0;
      //     for (int i=0;i<cref.dimension(0);++i)
      //       for (int j=0;j<cref.dimension(1);++j)
      //         for (int k=0;k<cref.dimension(2);++k)
      //           diff += std::abs(cref(i,j,k) - c(i,j,k));

      //     std::cout << std::setw(12) << "Plain"
      //               << " BlkSize = " << std::setw(3) << BlkSize
      //               << " time = " << std::scientific << tmin
      //               << " avg flop/s = " << (flop/tavg)
      //               << " max flop/s = " << (flop/tmin)
      //               << " diff to ref = " << diff
      //               << std::endl;
      //   }
      // }

      ///
      /// Serial SIMD with appropriate data layout
      ///
      {
        Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
          a("a", N, BlkSize, BlkSize),
          b("b", N, BlkSize, BlkSize),
          c("c", N, BlkSize, BlkSize);
        
        {
          const Kokkos::RangePolicy<DeviceSpaceType,ScheduleType> policy(0, N);
          
          double tavg = 0, tmin = tmax;

          for (int iter=iter_begin;iter<iter_end;++iter) {
            // flush
            if (!hot)
              flush.run();

            // initialize matrices
            if (!hot && iter == iter_begin) {
              Kokkos::deep_copy(a, amat_simd);
              Kokkos::deep_copy(b, bmat_simd);
            }
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

          double diff = 0;
          for (int i=0;i<cref.dimension(0);++i)
            for (int j=0;j<cref.dimension(1);++j)
              for (int k=0;k<cref.dimension(2);++k)
                diff += std::abs(cref(i,j,k) - c(i/VectorLength,j,k)[i%VectorLength]);

          std::cout << std::setw(12) << "Serial SIMD"
                    << " BlkSize = " << std::setw(3) << BlkSize
                    << " time = " << std::scientific << tmin
                    << " avg flop/s = " << (flop/tavg)
                    << " max flop/s = " << (flop/tmin)
                    << " diff to ref = " << diff
                    << std::endl;
        }
      }

    //   ///
    //   /// Team SIMD with appropriate data layout
    //   ///
    //   if (!std::is_same<AlgoTagType,Algo::Gemm::CompactMKL>::value) {
    //     Kokkos::View<VectorType***,Kokkos::LayoutRight,HostSpaceType> 
    //       a("a", N, BlkSize, BlkSize),
    //       b("b", N, BlkSize, BlkSize),
    //       c("c", N, BlkSize, BlkSize);

    //     {
    //       double tavg = 0, tmin = tmax;

    //       typedef Kokkos::TeamPolicy<DeviceSpaceType,ScheduleType> policy_type;
    //       typedef typename policy_type::member_type member_type;
    //       const policy_type policy(N, Kokkos::AUTO, VectorLength);
          
    //       for (int iter=iter_begin;iter<iter_end;++iter) {
    //         // flush
    //         flush.run();

    //         // initialize matrices
    //         Kokkos::deep_copy(a, amat_simd);
    //         Kokkos::deep_copy(b, bmat_simd);
    //         Kokkos::deep_copy(c, 0);

    //         DeviceSpaceType::fence();
    //         timer.reset();

    //         Kokkos::parallel_for
    //           (policy, 
    //            KOKKOS_LAMBDA(const member_type &member) {
    //             const int k = member.league_rank();

    //             auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
    //             auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
    //             auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());
                
    //             KokkosKernels::Team::
    //               Gemm<member_type,Trans::NoTranspose,Trans::NoTranspose,AlgoTagType>::
    //               invoke(member, 1.0, aa, bb, 1.0, cc);
    //           });
            
    //         DeviceSpaceType::fence();
    //         const double t = timer.seconds();
    //         tmin = std::min(tmin, t);
    //         tavg += (iter >= 0)*t;
    //       }
    //       tavg /= iter_end;

    //       double diff = 0;
    //       for (int i=0;i<cref.dimension(0);++i)
    //         for (int j=0;j<cref.dimension(1);++j)
    //           for (int k=0;k<cref.dimension(2);++k)
    //             diff += std::abs(cref(i,j,k) - c(i/VectorLength,j,k)[i%VectorLength]);

    //       std::cout << std::setw(12) << "Team SIMD"
    //                 << " BlkSize = " << std::setw(3) << BlkSize
    //                 << " time = " << std::scientific << tmin
    //                 << " avg flop/s = " << (flop/tavg)
    //                 << " max flop/s = " << (flop/tmin)
    //                 << " diff to ref = " << diff
    //                 << std::endl << std::endl;
    //     }
    //   }
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
    if (token == std::string("-hot-cache")) hot = true;
  }

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
      
      // std::cout << "\n Testing SIMD-" << VectorLength << " and Algo::Gemm::Unblocked\n";
      // run<VectorTag<SIMD<double>,VectorLength>,Algo::Gemm::Unblocked>(N[i]/VectorLength);
      
      // std::cout << "\n Testing AVX-" << VectorLength << " and Algo::Gemm::Unblocked\n";
      // run<VectorTag<AVX<double>,VectorLength>,Algo::Gemm::Unblocked>(N[i]/VectorLength);
      
      // std::cout << "\n Testing SIMD-" << VectorLength << " and Algo::Gemm::Blocked\n";
      // run<VectorTag<SIMD<double>,VectorLength>,Algo::Gemm::Blocked>(N[i]/VectorLength);
      
      std::cout << "\n Testing AVX-" << VectorLength << " and Algo::Gemm::Blocked\n";
      run<VectorTag<AVX<double>,VectorLength>,Algo::Gemm::Blocked>(N[i]/VectorLength);

// #if defined(__KOKKOSKERNELS_INTEL_MKL_COMPACT_BATCHED__)
//       std::cout << "\n Testing AVX-" << VectorLength << " and Algo::Gemm::CompactMKL\n";
//       run<VectorTag<AVX<double>,VectorLength>,Algo::Gemm::CompactMKL>(N[i]/VectorLength);
// #endif
    }
  }

  Kokkos::finalize();

  return 0;
}
