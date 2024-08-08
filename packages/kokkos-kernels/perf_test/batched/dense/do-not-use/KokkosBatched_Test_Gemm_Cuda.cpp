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

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"

#if defined(KOKKOS_ENABLE_CUDA)

#include <iomanip>

#include "KokkosBatched_Util.hpp"
#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cublas_api.h"
#endif

#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Copy_Impl.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_Gemm_Team_Impl.hpp"

namespace KokkosBatched {
namespace PerfTest {

#undef FLOP_MUL
#undef FLOP_ADD
#define FLOP_MUL 1.0
#define FLOP_ADD 1.0
typedef double value_type;

double FlopCount(int mm, int nn, int kk) {
  double m = (double)mm;
  double n = (double)nn;
  double k = (double)kk;
  return (FLOP_MUL * (m * n * k) + FLOP_ADD * (m * n * k));
}

struct RangeTag {};
struct TeamTagV1 {};
struct TeamTagV2 {};
struct TeamTagV3 {};
struct TeamTagHandmade {};

template <typename ViewType, typename AlgoTagType, int VectorLength = 0>
struct Functor {
  ConstUnmanagedViewType<ViewType> _a, _b;
  UnmanagedViewType<ViewType> _c;

  KOKKOS_INLINE_FUNCTION
  Functor() = default;

  KOKKOS_INLINE_FUNCTION
  Functor(const ViewType &a, const ViewType &b, const ViewType &c) : _a(a), _b(b), _c(c) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const RangeTag &, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());
    auto cc = Kokkos::subview(_c, k, Kokkos::ALL(), Kokkos::ALL());

    SerialGemm<Trans::NoTranspose, Trans::NoTranspose, AlgoTagType>::invoke(1.0, aa, bb, 1.0, cc);
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV1 &, const MemberType &member) const {
    const int kbeg = (member.league_rank() * (member.team_size() * VectorLength) + member.team_rank() * VectorLength);
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < int(_c.extent(0))) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        auto bb = Kokkos::subview(_b, kk, Kokkos::ALL(), Kokkos::ALL());
        auto cc = Kokkos::subview(_c, kk, Kokkos::ALL(), Kokkos::ALL());

        SerialGemm<Trans::NoTranspose, Trans::NoTranspose, AlgoTagType>::invoke(1.0, aa, bb, 1.0, cc);
      }
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV2 &, const MemberType &member) const {
    const int kbeg = member.league_rank() * VectorLength;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < int(_c.extent(0))) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        auto bb = Kokkos::subview(_b, kk, Kokkos::ALL(), Kokkos::ALL());
        auto cc = Kokkos::subview(_c, kk, Kokkos::ALL(), Kokkos::ALL());

        TeamGemm<MemberType, Trans::NoTranspose, Trans::NoTranspose, AlgoTagType>::invoke(member, 1.0, aa, bb, 1.0, cc);
      }
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV3 &, const MemberType &member) const {
    const int lvl = 0;
    ScratchViewType<ViewType> sa(member.team_scratch(lvl), VectorLength, _a.extent(1), _a.extent(2));
    ScratchViewType<ViewType> sb(member.team_scratch(lvl), VectorLength, _b.extent(1), _b.extent(2));

    const int kbeg = member.league_rank() * VectorLength;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < int(_c.extent(0))) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        auto bb = Kokkos::subview(_b, kk, Kokkos::ALL(), Kokkos::ALL());
        auto cc = Kokkos::subview(_c, kk, Kokkos::ALL(), Kokkos::ALL());

        auto saa = Kokkos::subview(sa, k, Kokkos::ALL(), Kokkos::ALL());
        auto sbb = Kokkos::subview(sb, k, Kokkos::ALL(), Kokkos::ALL());

        TeamCopy<MemberType, Trans::NoTranspose>::invoke(member, aa, saa);
        TeamCopy<MemberType, Trans::NoTranspose>::invoke(member, bb, sbb);
        member.team_barrier();

        TeamGemm<MemberType, Trans::NoTranspose, Trans::NoTranspose, AlgoTagType>::invoke(member, 1.0, saa, sbb, 1.0,
                                                                                          cc);
      }
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagHandmade &, const MemberType &member) const {
    const int kbeg = member.league_rank() * VectorLength;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < int(_c.extent(0))) {
        const int m = _c.extent(1), n = _c.extent(2), q = _a.extent(2);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, m * n), [&](const int &ij) {
          const int i = ij % m, j = ij / m;
          typename ViewType::non_const_value_type cval = 0;
          for (int p = 0; p < q; ++p) cval += _a(kk, i, p) * _b(kk, p, j);
          _c(kk, i, j) += cval;
        });
      }
    });
  }
};

template <typename DeviceSpaceType, typename AlgoTagType>
void Gemm(const int NN, const int BlkSize) {
  typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;

  constexpr int VectorLength = DefaultVectorLength<value_type, typename DeviceSpaceType::memory_space>::value;
  const int N                = NN / VectorLength;

  {
    std::string value_type_name;
    if (std::is_same<value_type, double>::value) value_type_name = "double";
    if (std::is_same<value_type, Kokkos::complex<double> >::value) value_type_name = "Kokkos::complex<double>";

    std::cout << "SIMD is defined: datatype " << value_type_name << " a vector length " << VectorLength << "\n";
  }

  const double flop = (N * VectorLength) * FlopCount(BlkSize, BlkSize, BlkSize);
  const double tmax = 1.0e15;

  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
  typedef typename DeviceSpaceType::memory_space DeviceMemorySpaceType;

  const int iter_begin = -3, iter_end = 30;
  Kokkos::Timer timer;

  Kokkos::View<value_type ***, Kokkos::LayoutLeft, HostSpaceType> amat("amat", N * VectorLength, BlkSize, BlkSize),
      bmat("bmat", N * VectorLength, BlkSize, BlkSize), cref("cref", N * VectorLength, BlkSize, BlkSize);

  {
    Random<value_type> random;
    for (int k = 0; k < N * VectorLength; ++k)
      for (int i = 0; i < BlkSize; ++i)
        for (int j = 0; j < BlkSize; ++j) {
          amat(k, i, j) = random.value();
          bmat(k, i, j) = random.value();
        }
  }

  // P100 L2 cache 4MB per core
  constexpr size_t LLC_CAPACITY = 56 * 4 * 1024 * 1024;
  Flush<LLC_CAPACITY, DeviceSpaceType> flush;

#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
  if (1) {
    ///
    /// CUBLAS Strided version
    ///
    const Kokkos::LayoutStride stride(N * VectorLength, BlkSize * BlkSize, BlkSize, 1, BlkSize, BlkSize);

    Kokkos::View<value_type ***, Kokkos::LayoutStride, DeviceSpaceType> a("a", stride), b("b", stride), c("c", stride);

    double tavg = 0, tmin = tmax;

    cublasStatus_t stat;
    cublasHandle_t handle;

    stat = cublasCreate(&handle);
    if (stat != CUBLAS_STATUS_SUCCESS) Kokkos::abort("CUBLAS initialization failed\n");

    auto amat_device = Kokkos::create_mirror_view(DeviceMemorySpaceType(), amat);
    auto bmat_device = Kokkos::create_mirror_view(DeviceMemorySpaceType(), bmat);

    Kokkos::deep_copy(amat_device, amat);
    Kokkos::deep_copy(bmat_device, bmat);

    Kokkos::fence();

    const double one(1.0), zero(0.0);
    {
      tavg = 0;
      tmin = tmax;

      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat_device);
        Kokkos::deep_copy(b, bmat_device);
        Kokkos::deep_copy(c, 0);

        Kokkos::fence();
        timer.reset();

        stat = cublasDgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, BlkSize, BlkSize, BlkSize, &one,
                                         (const value_type *)a.data(), BlkSize, BlkSize * BlkSize,
                                         (const value_type *)b.data(), BlkSize, BlkSize * BlkSize, &zero,
                                         (value_type *)c.data(), BlkSize, BlkSize * BlkSize, N * VectorLength);

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
      Kokkos::deep_copy(csol, c);
      Kokkos::deep_copy(cref, csol);

      std::cout << std::setw(8) << "CUBLAS" << std::setw(8) << "Strided"
                << " BlkSize = " << std::setw(3) << BlkSize << " TeamSize = N/A"
                << " ScratchSize (KB) =   0"
                << " time = " << std::scientific << tmin << " avg flop/s = " << (flop / tavg)
                << " max flop/s = " << (flop / tmin) << std::endl;
    }
    cublasDestroy(handle);
  }
#endif

  if (1) {
    ///
    /// Range policy version
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, BlkSize),
        c("c", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Functor<view_type, AlgoTagType> functor_type;
      const Kokkos::RangePolicy<DeviceSpaceType, ScheduleType, RangeTag> policy(0, N * VectorLength);

      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);
        Kokkos::deep_copy(c, 0);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::GemmCuda::RangeTag", policy, functor_type(a, b, c));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
      Kokkos::deep_copy(csol, c);

      double diff = 0;
      for (int i = 0, iend = cref.extent(0); i < iend; ++i)
        for (int j = 0, jend = cref.extent(1); j < jend; ++j)
          for (int k = 0, kend = cref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(cref(i, j, k) - csol(i, j, k));

      std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Range"
                << " BlkSize = " << std::setw(3) << BlkSize << " TeamSize = N/A"
                << " ScratchSize (KB) =   0"
                << " time = " << std::scientific << tmin << " avg flop/s = " << (flop / tavg)
                << " max flop/s = " << (flop / tmin);
#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
      std::cout << " diff to ref = " << diff;
#endif
      std::cout << std::endl;
    }
  }

  if (1) {
    ///
    /// Team policy V1 - almost same scheduling with range policy;
    ///                  expect the same performance as range policy
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, BlkSize),
        c("c", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV1> policy_type;

      typedef Functor<view_type, AlgoTagType, VectorLength> functor_type;

      // 128 is rough estimates
      const int team_size = policy_type(N / 32, Kokkos::AUTO, VectorLength)
                                .team_size_recommended(functor_type(), Kokkos::ParallelForTag());

      const policy_type policy(N / team_size, team_size, VectorLength);
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);
        Kokkos::deep_copy(c, 0);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::GemmCuda::TeamPolicyV1", policy, functor_type(a, b, c));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
      Kokkos::deep_copy(csol, c);

      double diff = 0;
      for (int i = 0, iend = cref.extent(0); i < iend; ++i)
        for (int j = 0, jend = cref.extent(1); j < jend; ++j)
          for (int k = 0, kend = cref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(cref(i, j, k) - csol(i, j, k));

      std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team V1"
                << " BlkSize = " << std::setw(3) << BlkSize << " TeamSize = " << std::setw(3) << team_size
                << " ScratchSize (KB) =   0"
                << " time = " << std::scientific << tmin << " avg flop/s = " << (flop / tavg)
                << " max flop/s = " << (flop / tmin);
#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
      std::cout << " diff to ref = " << diff;
#endif
      std::cout << std::endl;
    }
  }

  if (1) {
    ///
    /// Team policy V2 - team parallel
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, BlkSize),
        c("c", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV2> policy_type;
      typedef Functor<view_type, AlgoTagType, VectorLength> functor_type;

      const int is_blocked_algo = (std::is_same<AlgoTagType, Algo::Gemm::Blocked>::value),
                mb = Algo::Gemm::Blocked::mb<DeviceMemorySpaceType>(), mp = BlkSize % mb > 0;

      const int mblk = is_blocked_algo ? (BlkSize / mb + mp) : BlkSize;

      const int max_team_size =
          policy_type(N, Kokkos::AUTO, VectorLength).team_size_max(functor_type(), Kokkos::ParallelForTag());
      const int team_size = std::min(std::max(mblk * mblk, 4), max_team_size);

      policy_type policy(N, team_size, VectorLength);
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);
        Kokkos::deep_copy(c, 0);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::GemmCuda::TeamPolicyV2", policy, functor_type(a, b, c));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
      Kokkos::deep_copy(csol, c);

      double diff = 0;
      for (int i = 0, iend = cref.extent(0); i < iend; ++i)
        for (int j = 0, jend = cref.extent(1); j < jend; ++j)
          for (int k = 0, kend = cref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(cref(i, j, k) - csol(i, j, k));

      std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team V2"
                << " BlkSize = " << std::setw(3) << BlkSize << " TeamSize = " << std::setw(3) << team_size
                << " ScratchSize (KB) =   0"
                << " time = " << std::scientific << tmin << " avg flop/s = " << (flop / tavg)
                << " max flop/s = " << (flop / tmin);
#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
      std::cout << " diff to ref = " << diff;
#endif
      std::cout << std::endl;
    }
  }

  if (1) {
    ///
    /// Team policy V3 - team parallel + scratch
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, BlkSize),
        c("c", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV3> policy_type;
      typedef Functor<view_type, AlgoTagType, VectorLength> functor_type;

      const int lvl = 0, per_team_scratch = 2 * ScratchViewType<view_type>::shmem_size(VectorLength, BlkSize, BlkSize);
      // std::cout << "per team scratch " << per_team_scratch << "\n";
      if (per_team_scratch / 1024 < 48) {
        const int is_blocked_algo = (std::is_same<AlgoTagType, Algo::Gemm::Blocked>::value),
                  mb = Algo::Gemm::Blocked::mb<DeviceMemorySpaceType>(), mp = BlkSize % mb > 0;

        const int mblk = is_blocked_algo ? (BlkSize / mb + mp) : BlkSize;

        const int max_team_size = policy_type(N, Kokkos::AUTO, VectorLength)
                                      .set_scratch_size(lvl, Kokkos::PerTeam(per_team_scratch))
                                      .team_size_max(functor_type(), Kokkos::ParallelForTag());
        const int team_size = std::min(std::max(mblk * mblk, 4), max_team_size);

        policy_type policy =
            policy_type(N, team_size, VectorLength).set_scratch_size(lvl, Kokkos::PerTeam(per_team_scratch));
        for (int iter = iter_begin; iter < iter_end; ++iter) {
          // flush
          flush.run();

          // initialize matrices
          Kokkos::deep_copy(a, amat);
          Kokkos::deep_copy(b, bmat);
          Kokkos::deep_copy(c, 0);

          Kokkos::fence();
          timer.reset();

          Kokkos::parallel_for("KokkosBatched::PerfTest::GemmCuda::TeamPolicyV3", policy, functor_type(a, b, c));

          Kokkos::fence();
          const double t = timer.seconds();
          tmin           = std::min(tmin, t);
          tavg += (iter >= 0) * t;
        }
        tavg /= iter_end;

        auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
        Kokkos::deep_copy(csol, c);

        double diff = 0;
        for (int i = 0, iend = cref.extent(0); i < iend; ++i)
          for (int j = 0, jend = cref.extent(1); j < jend; ++j)
            for (int k = 0, kend = cref.extent(2); k < kend; ++k)
              diff += Kokkos::ArithTraits<value_type>::abs(cref(i, j, k) - csol(i, j, k));

        std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team V3"
                  << " BlkSize = " << std::setw(3) << BlkSize << " TeamSize = " << std::setw(3) << team_size
                  << " ScratchSize (KB) = " << std::setw(3) << (per_team_scratch / 1024)
                  << " time = " << std::scientific << tmin << " avg flop/s = " << (flop / tavg)
                  << " max flop/s = " << (flop / tmin);
#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
        std::cout << " diff to ref = " << diff;
#endif
        std::cout << std::endl;
      } else {
        std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team V3"
                  << " Scratch per team is too big:" << std::setw(3) << (per_team_scratch / 1024) << std::endl;
      }
    }
  }

  if (1) {
    ///
    /// Team policy - handmade
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, BlkSize),
        c("c", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagHandmade> policy_type;
      typedef Functor<view_type, AlgoTagType, VectorLength> functor_type;

      const int max_team_size =
          policy_type(N, Kokkos::AUTO, VectorLength).team_size_max(functor_type(), Kokkos::ParallelForTag());

      const int team_size = std::min(max_team_size, BlkSize * BlkSize);

      const policy_type policy(N, team_size, VectorLength);
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);
        Kokkos::deep_copy(c, 0);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::GemmCuda::TeamPolicyHandmade", policy, functor_type(a, b, c));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto csol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
      Kokkos::deep_copy(csol, c);

      double diff = 0;
      for (int i = 0, iend = cref.extent(0); i < iend; ++i)
        for (int j = 0, jend = cref.extent(1); j < jend; ++j)
          for (int k = 0, kend = cref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(cref(i, j, k) - csol(i, j, k));

      std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team HM"
                << " BlkSize = " << std::setw(3) << BlkSize << " TeamSize = " << std::setw(3) << team_size
                << " ScratchSize (KB) =   0"
                << " time = " << std::scientific << tmin << " avg flop/s = " << (flop / tavg)
                << " max flop/s = " << (flop / tmin);
#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
      std::cout << " diff to ref = " << diff;
#endif
      std::cout << std::endl;
    }
  }

  std::cout << std::endl;
}
}  // namespace PerfTest
}  // namespace KokkosBatched

using namespace KokkosBatched;

template <typename AlgoTagType>
void run(const int N, const int B) {
  typedef Kokkos::DefaultExecutionSpace ExecSpace;

  Kokkos::print_configuration(std::cout);

  if (B != 0) {
    PerfTest::Gemm<ExecSpace, AlgoTagType>(N, B);
  } else {
    PerfTest::Gemm<ExecSpace, AlgoTagType>(N, 3);
    PerfTest::Gemm<ExecSpace, AlgoTagType>(N, 5);
    PerfTest::Gemm<ExecSpace, AlgoTagType>(N, 10);
    PerfTest::Gemm<ExecSpace, AlgoTagType>(N, 15);

    // PerfTest::Gemm<ExecSpace, AlgoTagType>(N,  4);
    // PerfTest::Gemm<ExecSpace, AlgoTagType>(N,  8);
    // PerfTest::Gemm<ExecSpace, AlgoTagType>(N, 16);
    // PerfTest::Gemm<ExecSpace, AlgoTagType>(N, 18);
  }
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);

  int N = 128 * 128, B = 0;

  for (int i = 1; i < argc; ++i) {
    const std::string &token = argv[i];
    if (token == std::string("-N")) N = std::atoi(argv[++i]);
    if (token == std::string("-B")) B = std::atoi(argv[++i]);
  }

  {
    std::cout << " N = " << N << std::endl;

    std::cout << "\n Testing LayoutLeft Algo::Gemm::Unblocked\n";
    run<Algo::Gemm::Unblocked>(N, B);

    std::cout << "\n Testing LayoutLeft Algo::Gemm::Blocked\n";
    run<Algo::Gemm::Blocked>(N, B);
  }

  Kokkos::finalize();

  return 0;
}

#else

int main(int argc, char *argv[]) {
  std::cout << "Kokkos::Cuda is not enabled\n";
  return -1;
}
#endif
