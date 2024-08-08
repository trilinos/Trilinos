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

// #define __KOKKOSBATCHED_INTEL_MKL__
// #define __KOKKOSBATCHED_INTEL_MKL_BATCHED__

#include <iomanip>

#include "KokkosBatched_Util.hpp"
#if defined(__KOKKOSBATCHED_INTEL_MKL__)
#include "mkl.h"
#endif

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"
#include "Kokkos_Random.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Gemv_Decl.hpp"
#include "KokkosBatched_Gemv_Serial_Impl.hpp"

#undef __KOKKOSBATCHED_INTEL_MKL_BATCHED__

namespace KokkosBatched {
namespace PerfTest {

#undef FLOP_MUL
#undef FLOP_ADD

#if defined(KokkosBatched_Test_Gemv_Host_Complex)
#define FLOP_MUL 6.0
#define FLOP_ADD 2.0
typedef Kokkos::complex<double> value_type;
#endif

#if defined(KokkosBatched_Test_Gemv_Host_Real)
#define FLOP_MUL 1.0
#define FLOP_ADD 1.0
typedef double value_type;
#endif

double FlopCount(int mm, int nn) {
  double m = (double)mm;
  double n = (double)nn;
  return (FLOP_MUL * (m * n) + FLOP_ADD * (m * n));
}

template <int BlkSize, int NumVecs, typename HostSpaceType, typename AlgoTagType>
void Gemv(const int NN) {
  typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;
  // typedef Kokkos::Schedule<Kokkos::Dynamic> ScheduleType;

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

  const double flop = (N * VectorLength) * FlopCount(BlkSize, BlkSize) * NumVecs;
  // const double tmax = 1.0e15;

  const int iter_begin = -10, iter_end = 100;
  Kokkos::Timer timer;

  Kokkos::View<value_type***, Kokkos::LayoutRight, HostSpaceType> yref;
  Kokkos::View<value_type***, Kokkos::LayoutRight, HostSpaceType> amat("amat", N * VectorLength, BlkSize, BlkSize);
  Kokkos::View<value_type***, Kokkos::LayoutRight, HostSpaceType> xvec("xvec", N * VectorLength, NumVecs, BlkSize);

  Kokkos::Random_XorShift64_Pool<HostSpaceType> random(13718);
  Kokkos::fill_random(xvec, random, value_type(1.0));
  Kokkos::fill_random(amat, random, value_type(1.0));

  // for KNL
  constexpr size_t LLC_CAPACITY = 34 * 1024 * 1024;
  Flush<LLC_CAPACITY> flush;

  ///
  /// Reference version using MKL DGEMM
  ///
#if defined(__KOKKOSBATCHED_INTEL_MKL__)
  {
    Kokkos::View<value_type***, Kokkos::LayoutRight, HostSpaceType> a("a", N * VectorLength, BlkSize, BlkSize),
        x("x", N * VectorLength, NumVecs, BlkSize), y("y", N * VectorLength, NumVecs, BlkSize);

    {
      const Kokkos::RangePolicy<HostSpaceType, ScheduleType> policy(0, N * VectorLength);

      double t = 0;
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(x, xvec);
        Kokkos::deep_copy(y, 0);

        HostSpaceType().fence();
        timer.reset();

        Kokkos::parallel_for(
            "KokkosBatched::PerfTest::GemvHost::CblasOpenMP", policy, KOKKOS_LAMBDA(const int k) {
              auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());
              for (int j = 0; j < NumVecs; ++j) {
                auto xx = Kokkos::subview(x, k, j, Kokkos::ALL());
                auto yy = Kokkos::subview(y, k, j, Kokkos::ALL());

                cblas_dgemv(CblasRowMajor, CblasNoTrans, BlkSize, BlkSize, 1.0, (double*)aa.data(), aa.stride_0(),
                            (double*)xx.data(), xx.stride_0(), 1.0, (double*)yy.data(), yy.stride_0());
              }
            });

        HostSpaceType().fence();
        t += (iter >= 0) * timer.seconds();
      }
      t /= iter_end;

      std::cout << std::setw(12) << "MKL DGEMV"
                << " BlkSize = " << std::setw(3) << BlkSize << " NumVecs = " << std::setw(3) << NumVecs
                << " time = " << std::scientific << t << " flop/s = " << (flop / t) << std::endl;

      yref = y;
    }
  }
#endif

  ///
  /// Plain version (comparable to micro BLAS version)
  ///
  {
    Kokkos::View<value_type***, Kokkos::LayoutRight, HostSpaceType> a("a", N * VectorLength, BlkSize, BlkSize),
        x("x", N * VectorLength, NumVecs, BlkSize), y("y", N * VectorLength, NumVecs, BlkSize);

    {
      const Kokkos::RangePolicy<HostSpaceType, ScheduleType> policy(0, N * VectorLength);

      double t = 0;
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(x, xvec);
        Kokkos::deep_copy(y, 0);

        HostSpaceType().fence();
        timer.reset();

        Kokkos::parallel_for(
            "KokkosBatched::PerfTest::GemvHost::SerialOpenMP", policy, KOKKOS_LAMBDA(const int k) {
              auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());

              for (int j = 0; j < NumVecs; ++j) {
                auto xx = Kokkos::subview(x, k, j, Kokkos::ALL());
                auto yy = Kokkos::subview(y, k, j, Kokkos::ALL());

                SerialGemv<Trans::NoTranspose, AlgoTagType>::invoke(1.0, aa, xx, 1.0, yy);
              }
            });

        HostSpaceType().fence();
        t += (iter >= 0) * timer.seconds();
      }
      t /= iter_end;

      double diff = 0;
      for (int i = 0, iend = yref.extent(0); i < iend; ++i)
        for (int j = 0, jend = yref.extent(1); j < jend; ++j)
          for (int k = 0, kend = yref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(yref(i, j, k) - y(i, j, k));

      std::cout << std::setw(12) << "Plain"
                << " BlkSize = " << std::setw(3) << BlkSize << " NumVecs = " << std::setw(3) << NumVecs
                << " time = " << std::scientific << t << " flop/s = " << (flop / t) << " diff to ref = " << diff
                << std::endl;
    }
  }

  typedef Vector<SIMD<value_type>, VectorLength> VectorType;
  Kokkos::View<VectorType***, Kokkos::LayoutRight, HostSpaceType> amat_simd("amat_simd", N, BlkSize, BlkSize),
      xvec_simd("xvec_simd", N, NumVecs, BlkSize);

  for (int k0 = 0; k0 < N; ++k0)
    for (int k1 = 0; k1 < VectorLength; ++k1)
      for (int i = 0; i < BlkSize; ++i) {
        for (int j = 0; j < NumVecs; ++j) xvec_simd(k0, j, i)[k1] = xvec(k0 * VectorLength + k1, j, i);
        for (int j = 0; j < BlkSize; ++j) amat_simd(k0, i, j)[k1] = amat(k0 * VectorLength + k1, i, j);
      }

  ///
  /// Serial SIMD with appropriate data layout
  ///
  {
    Kokkos::View<VectorType***, Kokkos::LayoutRight, HostSpaceType> a("a", N, BlkSize, BlkSize),
        x("x", N, NumVecs, BlkSize), y("y", N, NumVecs, BlkSize);

    {
      const Kokkos::RangePolicy<HostSpaceType, ScheduleType> policy(0, N);

      double t = 0;
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat_simd);
        Kokkos::deep_copy(x, xvec_simd);
        Kokkos::deep_copy(y, 0);

        HostSpaceType().fence();
        timer.reset();

        Kokkos::parallel_for(
            "KokkosBatched::PerfTest::GemvHost::SIMDSerialOpenMP", policy, KOKKOS_LAMBDA(const int k) {
              auto aa = Kokkos::subview(a, k, Kokkos::ALL(), Kokkos::ALL());

              for (int j = 0; j < NumVecs; ++j) {
                auto xx = Kokkos::subview(x, k, j, Kokkos::ALL());
                auto yy = Kokkos::subview(y, k, j, Kokkos::ALL());

                SerialGemv<Trans::NoTranspose, AlgoTagType>::invoke(1.0, aa, xx, 1.0, yy);
              }
            });

        HostSpaceType().fence();
        t += (iter >= 0) * timer.seconds();
      }
      t /= iter_end;

      double diff = 0;
      for (int i = 0, iend = yref.extent(0); i < iend; ++i)
        for (int j = 0, jend = yref.extent(1); j < jend; ++j)
          for (int k = 0, kend = yref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(yref(i, j, k) - y(i / VectorLength, j, k)[i % VectorLength]);

      std::cout << std::setw(12) << "Serial SIMD"
                << " BlkSize = " << std::setw(3) << BlkSize << " NumVecs = " << std::setw(3) << NumVecs
                << " time = " << std::scientific << t << " flop/s = " << (flop / t) << " diff to ref = " << diff
                << std::endl;
    }
  }
}

}  // namespace PerfTest
}  // namespace KokkosBatched
