//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>

#include "KokkosBatched_Util.hpp"
#if defined(__KOKKOSBATCHED_LIBXSMM__)
#include "libxsmm.h"
#endif

#if defined(__KOKKOSBATCHED_INTEL_MKL__)
#include "mkl.h"
#endif

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"

// #undef __KOKKOSBATCHED_INTEL_MKL_BATCHED__

namespace KokkosBatched {
namespace PerfTest {

#undef FLOP_MUL
#undef FLOP_ADD

#if defined(KokkosBatched_Test_Gemm_Host_Complex)
#define FLOP_MUL 6.0
#define FLOP_ADD 2.0
typedef Kokkos::complex<double> value_type;
#endif

#if defined(KokkosBatched_Test_Gemm_Host_Real)
#define FLOP_MUL 1.0
#define FLOP_ADD 1.0
typedef double value_type;
#endif

double FlopCount(int mm, int nn, int kk) {
  double m = (double)mm;
  double n = (double)nn;
  double k = (double)kk;
  return (FLOP_MUL * (m * n * k) + FLOP_ADD * (m * n * k));
}

template <int BlkSize, typename HostSpaceType, typename AlgoTagType>
void Gemm(const int NN) {
  typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;

  constexpr int VectorLength = DefaultVectorLength<value_type, typename HostSpaceType::memory_space>::value;
  const int N                = NN / VectorLength;

  {
    std::string value_type_name;
    if (std::is_same<value_type, double>::value) value_type_name = "double";
    if (std::is_same<value_type, Kokkos::complex<double> >::value) value_type_name = "Kokkos::complex<double>";
#if defined(__AVX512F__)
    std::cout << "AVX512 is defined: datatype " << value_type_name << " a vector length " << VectorLength << "\n";
#elif defined(__AVX__) || defined(__AVX2__)
    std::cout << "AVX or AVX2 is defined: datatype " << value_type_name << " a vector length " << VectorLength << "\n";
#else
    std::cout << "SIMD (compiler vectorization) is defined: datatype " << value_type_name << " a vector length "
              << VectorLength << "\n";
#endif
  }

  const double flop = (N * VectorLength) * FlopCount(BlkSize, BlkSize, BlkSize);
  const double tmax = 1.0e15;

  const int iter_begin = -10, iter_end = 100;
  Kokkos::Timer timer;

  Kokkos::View<value_type ***, Kokkos::LayoutRight, HostSpaceType> cref;
  Kokkos::View<value_type ***, Kokkos::LayoutRight, HostSpaceType> amat("amat", N * VectorLength, BlkSize, BlkSize),
      bmat("bmat", N * VectorLength, BlkSize, BlkSize);

  Kokkos::Random_XorShift64_Pool<HostSpaceType> random(13718);
  Kokkos::fill_random(amat, random, value_type(1.0));
  Kokkos::fill_random(bmat, random, value_type(1.0));

  typedef Vector<SIMD<value_type>, VectorLength> VectorType;
  Kokkos::View<VectorType ***, Kokkos::LayoutRight, HostSpaceType> amat_simd("amat_simd", N, BlkSize, BlkSize),
      bmat_simd("bmat_simd", N, BlkSize, BlkSize);

  Kokkos::parallel_for(
      "KokkosBatched::PerfTest::GemmHost::Pack", Kokkos::RangePolicy<HostSpaceType>(0, N * VectorLength),
      KOKKOS_LAMBDA(const int k) {
        const int k0 = k / VectorLength, k1 = k % VectorLength;
        for (int i = 0; i < BlkSize; ++i)
          for (int j = 0; j < BlkSize; ++j) {
            amat_simd(k0, i, j)[k1] = amat(k, i, j);
            bmat_simd(k0, i, j)[k1] = bmat(k, i, j);
          }
      });

  // for KNL (1MB per tile)
  constexpr size_t LLC_CAPACITY = 34 * 1024 * 1024;
  Flush<LLC_CAPACITY> flush;

  ///
  /// Reference version using MKL DGEMM
  ///
#if defined(__KOKKOSBATCHED_INTEL_MKL__)
  {
    Kokkos::View<value_type ***, Kokkos::LayoutRight, HostSpaceType> a("a", N * VectorLength, BlkSize, BlkSize),
        b("b", N * VectorLength, BlkSize, BlkSize), c("c", N * VectorLength, BlkSize, BlkSize);

    {
      const Kokkos::RangePolicy<HostSpaceType, ScheduleType> policy(0, N * VectorLength);

      double tavg = 0, tmin = tmax;
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);
        Kokkos::deep_copy(c, 0);

        HostSpaceType().fence();
        timer.reset();

        Kokkos::parallel_for(
            "KokkosBatched::PerfTest::GemmHost::CblasOpenMP", policy, KOKKOS_LAMBDA(const int k) {
              auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
              auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
              auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());

              const double one = 1.0;
              if (std::is_same<value_type, double>::value) {
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, BlkSize, BlkSize, BlkSize, one,
                            (double *)aa.data(), aa.stride_0(), (double *)bb.data(), bb.stride_0(), one,
                            (double *)cc.data(), cc.stride_0());
              } else if (std::is_same<value_type, Kokkos::complex<double> >::value) {
                cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, BlkSize, BlkSize, BlkSize, (void *)&one,
                            (void *)aa.data(), aa.stride_0(), (void *)bb.data(), bb.stride_0(), (void *)&one,
                            (void *)cc.data(), cc.stride_0());
              }
            });

        HostSpaceType().fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      std::cout << std::setw(12) << "MKL DGEMM"
                << " BlkSize = " << std::setw(3) << BlkSize << " time = " << std::scientific << tmin
                << " avg flop/s = " << (flop / tavg) << " max flop/s = " << (flop / tmin) << std::endl;

      cref = c;
    }
  }

#if defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__)
  {
    typedef Kokkos::View<value_type ***, Kokkos::LayoutRight, HostSpaceType> ViewType;
    ViewType a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, BlkSize),
        c("c", N * VectorLength, BlkSize, BlkSize);

    value_type *aa[N * VectorLength], *bb[N * VectorLength], *cc[N * VectorLength];

    for (int k = 0; k < N * VectorLength; ++k) {
      aa[k] = &a(k, 0, 0);
      bb[k] = &b(k, 0, 0);
      cc[k] = &c(k, 0, 0);
    }

    {
      double tavg = 0, tmin = tmax;

      MKL_INT blksize[1] = {BlkSize};
      MKL_INT lda[1]     = {a.stride_1()};
      MKL_INT ldb[1]     = {b.stride_1()};
      MKL_INT ldc[1]     = {c.stride_1()};

      CBLAS_TRANSPOSE transA[1] = {CblasNoTrans};
      CBLAS_TRANSPOSE transB[1] = {CblasNoTrans};

      double one[1]           = {1.0};
      MKL_INT size_per_grp[1] = {N * VectorLength};

      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);
        Kokkos::deep_copy(c, 0);

        HostSpaceType().fence();
        timer.reset();

        if (std::is_same<value_type, double>::value) {
          cblas_dgemm_batch(CblasRowMajor, transA, transB, blksize, blksize, blksize, one, (const double **)aa, lda,
                            (const double **)bb, ldb, one, (double **)cc, ldc, 1, size_per_grp);
        } else if (std::is_same<value_type, Kokkos::complex<double> >::value) {
          cblas_zgemm_batch(CblasRowMajor, transA, transB, blksize, blksize, blksize, one, (const void **)aa, lda,
                            (const void **)bb, ldb, one, (void **)cc, ldc, 1, size_per_grp);
        }

        HostSpaceType().fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      double diff = 0;
      for (int i = 0, iend = cref.extent(0); i < iend; ++i)
        for (int j = 0, jend = cref.extent(1); j < jend; ++j)
          for (int k = 0, kend = cref.extent(2); k < kend; ++k) diff += abs(cref(i, j, k) - c(i, j, k));

      std::cout << std::setw(12) << "MKL Batch"
                << " BlkSize = " << std::setw(3) << BlkSize << " time = " << std::scientific << tmin
                << " avg flop/s = " << (flop / tavg) << " max flop/s = " << (flop / tmin) << " diff to ref = " << diff
                << std::endl;
    }
  }
#endif
#if defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
  {
    Kokkos::View<VectorType ***, Kokkos::LayoutRight, HostSpaceType> a("a", N, BlkSize, BlkSize),
        b("b", N, BlkSize, BlkSize), c("c", N, BlkSize, BlkSize);

    {
      double tavg = 0, tmin = tmax;

      double done(1.0);
      std::complex<double> zone(1.0);

      MKL_COMPACT_PACK format;
      if (std::is_same<value_type, double>::value) {
        if (VectorLength == 4)
          format = MKL_COMPACT_AVX;
        else if (VectorLength == 8)
          format = MKL_COMPACT_AVX512;
      } else if (std::is_same<value_type, Kokkos::complex<double> >::value) {
        if (VectorLength == 2)
          format = MKL_COMPACT_AVX;
        else if (VectorLength == 4)
          format = MKL_COMPACT_AVX512;
      }

      if (format == MKL_COMPACT_AVX512 || format == MKL_COMPACT_AVX) {
        for (int iter = iter_begin; iter < iter_end; ++iter) {
          // flush
          flush.run();

          // initialize matrices
          Kokkos::deep_copy(a, amat_simd);
          Kokkos::deep_copy(b, bmat_simd);
          Kokkos::deep_copy(c, 0);

          HostSpaceType().fence();
          timer.reset();

          if (std::is_same<value_type, double>::value) {
            mkl_dgemm_compact(MKL_ROW_MAJOR, MKL_NOTRANS, MKL_NOTRANS, BlkSize, BlkSize, BlkSize, done,
                              (const double *)a.data(), (MKL_INT)a.stride_1(), (const double *)b.data(),
                              (MKL_INT)b.stride_1(), done, (double *)c.data(), (MKL_INT)c.stride_1(), format,
                              N * VectorLength);
          } else if (std::is_same<value_type, Kokkos::complex<double> >::value) {
            mkl_zgemm_compact(MKL_ROW_MAJOR, MKL_NOTRANS, MKL_NOTRANS, BlkSize, BlkSize, BlkSize,
                              (MKL_Complex16 *)&zone, (const double *)a.data(), (MKL_INT)a.stride_1(),
                              (const double *)b.data(), (MKL_INT)b.stride_1(), (MKL_Complex16 *)&zone,
                              (double *)c.data(), (MKL_INT)c.stride_1(), format, N * VectorLength);
          }

          HostSpaceType().fence();
          const double t = timer.seconds();
          tmin           = std::min(tmin, t);
          tavg += (iter >= 0) * t;
        }
        tavg /= iter_end;

        double diff = 0;
        for (int i = 0, iend = cref.extent(0); i < iend; ++i)
          for (int j = 0, jend = cref.extent(1); j < jend; ++j)
            for (int k = 0, kend = cref.extent(2); k < kend; ++k)
              diff += abs(cref(i, j, k) - c(i / VectorLength, j, k)[i % VectorLength]);

        std::cout << std::setw(12) << "MKL Cmpct"
                  << " BlkSize = " << std::setw(3) << BlkSize << " time = " << std::scientific << tmin
                  << " avg flop/s = " << (flop / tavg) << " max flop/s = " << (flop / tmin) << " diff to ref = " << diff
                  << std::endl;
      }
    }
  }
#endif
#endif

#if defined(__KOKKOSBATCHED_LIBXSMM__)
  {
    libxsmm_init();

    Kokkos::View<value_type ***, Kokkos::LayoutRight, HostSpaceType> a("a", N * VectorLength, BlkSize, BlkSize),
        b("b", N * VectorLength, BlkSize, BlkSize), c("c", N * VectorLength, BlkSize, BlkSize);

    libxsmm_blasint lda = a.stride_1(), ldb = b.stride_1(), ldc = c.stride_1();

    {
      const Kokkos::RangePolicy<HostSpaceType, ScheduleType> policy(0, N * VectorLength);

      double tavg = 0, tmin = tmax;

      // adjust column major order in xsmm
      char transA = 'N', transB = 'N';
      libxsmm_blasint blksize = BlkSize;
      double one              = 1.0;

      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);
        Kokkos::deep_copy(c, 0);

        HostSpaceType().fence();
        timer.reset();

        Kokkos::parallel_for(
            "KokkosBatched::PerfTest::GemmHost::libxswmmOpenMP", policy, KOKKOS_LAMBDA(const int k) {
              auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
              auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
              auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());

              // column major
              libxsmm_gemm((const char *)&transA, (const char *)&transB, blksize, blksize, blksize,
                           (const double *)&one, (const double *)bb.data(), (const libxsmm_blasint *)&ldb,
                           (const double *)aa.data(), (const libxsmm_blasint *)&lda, (const double *)&one,
                           (double *)cc.data(), (const libxsmm_blasint *)&ldc);
            });

        HostSpaceType().fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      // adjust transpose
      double diff = 0;
      for (int i = 0, iend = cref.extent(0); i < iend; ++i)
        for (int j = 0, jend = cref.extent(1); j < jend; ++j)
          for (int k = 0, kend = cref.extent(2); k < kend; ++k) diff += abs(cref(i, j, k) - c(i, j, k));

      std::cout << std::setw(12) << "libxsmm"
                << " BlkSize = " << std::setw(3) << BlkSize << " time = " << std::scientific << tmin
                << " avg flop/s = " << (flop / tavg) << " max flop/s = " << (flop / tmin) << " diff to ref = " << diff
                << std::endl;
    }
    libxsmm_finalize();
  }
#endif
  ///
  /// Do not test this. Test Compact vs MKL
  /// KK Scalar version (comparable to micro BLAS version)
  ///
  // if (!std::is_same<AlgoTagType,Algo::Gemm::CompactMKL>::value) {
  //   Kokkos::View<value_type***,Kokkos::LayoutRight,HostSpaceType>
  //     a("a", N*VectorLength, BlkSize, BlkSize),
  //     b("b", N*VectorLength, BlkSize, BlkSize),
  //     c("c", N*VectorLength, BlkSize, BlkSize);

  //   {
  //     const Kokkos::RangePolicy<HostSpaceType,ScheduleType> policy(0,
  //     N*VectorLength);

  //     double tavg = 0, tmin = tmax;

  //     for (int iter=iter_begin;iter<iter_end;++iter) {
  //       // flush
  //       flush.run();

  //       // initialize matrices
  //       Kokkos::deep_copy(a, amat);
  //       Kokkos::deep_copy(b, bmat);
  //       Kokkos::deep_copy(c, 0);

  //       HostSpaceType().fence();
  //       timer.reset();

  //       Kokkos::parallel_for
  //         (policy,
  //          KOKKOS_LAMBDA(const int k) {
  //           auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
  //           auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
  //           auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());

  //           SerialGemm<Trans::NoTranspose,Trans::NoTranspose,AlgoTagType>::
  //             invoke(1.0, aa, bb, 1.0, cc);
  //         });

  //       HostSpaceType().fence();
  //       const double t = timer.seconds();
  //       tmin = std::min(tmin, t);
  //       tavg += (iter >= 0)*t;
  //     }
  //     tavg /= iter_end;

  //     double diff = 0;
  //     for (int i=0,iend=cref.extent(2);i<iend;++i)
  //       for (int j=0,jend=cref.extent(2);j<jend;++j)
  //         for (int k=0,kend=cref.extent(2);k<kend;++k)
  //           diff += abs(cref(i,j,k) - c(i,j,k));

  //     std::cout << std::setw(12) << "KK Scalar"
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
    Kokkos::View<VectorType ***, Kokkos::LayoutRight, HostSpaceType> a("a", N, BlkSize, BlkSize),
        b("b", N, BlkSize, BlkSize), c("c", N, BlkSize, BlkSize);

    {
      const Kokkos::RangePolicy<HostSpaceType, ScheduleType> policy(0, N);

      double tavg = 0, tmin = tmax;

      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat_simd);
        Kokkos::deep_copy(b, bmat_simd);
        Kokkos::deep_copy(c, 0);

        HostSpaceType().fence();
        timer.reset();

        Kokkos::parallel_for(
            "KokkosBatched::PerfTest::GemmHost::SIMDSerialOpenMP", policy, KOKKOS_LAMBDA(const int k) {
              auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
              auto bb = Kokkos::subview(b, k, Kokkos::ALL(), Kokkos::ALL());
              auto cc = Kokkos::subview(c, k, Kokkos::ALL(), Kokkos::ALL());

              SerialGemm<Trans::NoTranspose, Trans::NoTranspose, AlgoTagType>::invoke(1.0, aa, bb, 1.0, cc);
            });

        HostSpaceType().fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      double diff = 0;
      for (int i = 0, iend = cref.extent(0); i < iend; ++i)
        for (int j = 0, jend = cref.extent(1); j < jend; ++j)
          for (int k = 0, kend = cref.extent(2); k < kend; ++k)
            diff += abs(cref(i, j, k) - c(i / VectorLength, j, k)[i % VectorLength]);

      std::cout << std::setw(12) << "KK Vector"
                << " BlkSize = " << std::setw(3) << BlkSize << " time = " << std::scientific << tmin
                << " avg flop/s = " << (flop / tavg) << " max flop/s = " << (flop / tmin) << " diff to ref = " << diff
                << std::endl;
    }
  }
  std::cout << std::endl;
}

}  // namespace PerfTest
}  // namespace KokkosBatched
