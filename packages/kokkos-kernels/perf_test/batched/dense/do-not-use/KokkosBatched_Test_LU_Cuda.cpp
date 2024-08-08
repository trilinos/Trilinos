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

#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cublas_api.h"
#endif

#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Copy_Impl.hpp"

#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_LU_Team_Impl.hpp"

namespace KokkosBatched {
namespace PerfTest {

#define FLOP_MUL 1.0
#define FLOP_ADD 1.0
typedef double value_type;

double FlopCount(int mm, int nn) {
  double m = (double)mm;
  double n = (double)nn;
  if (m > n)
    return (FLOP_MUL * (0.5 * m * n * n - (1.0 / 6.0) * n * n * n + 0.5 * m * n - 0.5 * n * n + (2.0 / 3.0) * n) +
            FLOP_ADD * (0.5 * m * n * n - (1.0 / 6.0) * n * n * n - 0.5 * m * n + (1.0 / 6.0) * n));
  else
    return (FLOP_MUL * (0.5 * n * m * m - (1.0 / 6.0) * m * m * m + 0.5 * n * m - 0.5 * m * m + (2.0 / 3.0) * m) +
            FLOP_ADD * (0.5 * n * m * m - (1.0 / 6.0) * m * m * m - 0.5 * n * m + (1.0 / 6.0) * m));
}

struct RangeTag {};
struct TeamTagV1 {};
struct TeamTagV2 {};
struct TeamTagV3 {};
struct TeamTagHandmade {};

template <typename ViewType, typename AlgoTagType, int VectorLength = 0>
struct Functor {
  UnmanagedViewType<ViewType> _a;

  KOKKOS_INLINE_FUNCTION
  Functor() = default;

  KOKKOS_INLINE_FUNCTION
  Functor(const ViewType &a) : _a(a) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const RangeTag &, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    SerialLU<AlgoTagType>::invoke(aa);
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV1 &, const MemberType &member) const {
    const int kbeg = (member.league_rank() * (member.team_size() * VectorLength) + member.team_rank() * VectorLength);
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < _a.extent_int(0)) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        SerialLU<AlgoTagType>::invoke(aa);
      }
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV2 &, const MemberType &member) const {
    const int kbeg = member.league_rank() * VectorLength;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < _a.extent_int(0)) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        TeamLU<MemberType, AlgoTagType>::invoke(member, aa);
      }
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV3 &, const MemberType &member) const {
    const int lvl = 0;
    ScratchViewType<ViewType> sa(member.team_scratch(lvl), VectorLength, _a.extent(1), _a.extent(2));

    const int kbeg = member.league_rank() * VectorLength;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < _a.extent_int(0)) {
        auto aa  = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        auto saa = Kokkos::subview(sa, k, Kokkos::ALL(), Kokkos::ALL());

        TeamCopy<MemberType, Trans::NoTranspose>::invoke(member, aa, saa);
        member.team_barrier();
        TeamLU<MemberType, AlgoTagType>::invoke(member, saa);
        member.team_barrier();
        TeamCopy<MemberType, Trans::NoTranspose>::invoke(member, saa, aa);
      }
    });
  }
};

template <typename DeviceSpaceType, typename AlgoTagType>
void LU(const int NN, const int BlkSize) {
  typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;

  constexpr int VectorLength = DefaultVectorLength<value_type, typename DeviceSpaceType::memory_space>::value;
  const int N                = NN / VectorLength;

  {
    std::string value_type_name;
    if (std::is_same<value_type, double>::value) value_type_name = "double";
    if (std::is_same<value_type, Kokkos::complex<double> >::value) value_type_name = "Kokkos::complex<double>";

    std::cout << "SIMD is defined: datatype " << value_type_name << " a vector length " << VectorLength << "\n";
  }

  const double flop = (N * VectorLength) * FlopCount(BlkSize, BlkSize);
  const double tmax = 1.0e15;

  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
  typedef typename DeviceSpaceType::memory_space DeviceMemorySpaceType;

  const int iter_begin = -3, iter_end = 50;
  Kokkos::Timer timer;

  Kokkos::View<value_type ***, Kokkos::LayoutLeft, HostSpaceType> amat("amat", N * VectorLength, BlkSize, BlkSize),
      aref("aref", N * VectorLength, BlkSize, BlkSize);

  {
    Random<value_type> random;
    for (int k = 0; k < N * VectorLength; ++k) {
      // use tridiagonal matrices; for now we just check elementwise l/u factors
      // do not allow pivots
      for (int i = 0; i < BlkSize; ++i) {
        amat(k, i, i) = random.value() + 10.0;
        if ((i + 1) < BlkSize) {
          amat(k, i, i + 1) = random.value() + 1.0;
          amat(k, i + 1, i) = random.value() + 1.0;
        }
      }

      // value_type d[BlkSize], v[BlkSize][BlkSize];
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

  constexpr size_t LLC_CAPACITY = 56 * 4 * 1024 * 1024;
  Flush<LLC_CAPACITY> flush;

#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
  if (1) {
    ///
    /// CUBLAS Batch version
    ///
    const Kokkos::LayoutStride stride(N * VectorLength, BlkSize * BlkSize, BlkSize, 1, BlkSize, BlkSize);

    Kokkos::View<value_type ***, Kokkos::LayoutStride, DeviceSpaceType> a("a", stride);
    Kokkos::View<int *, DeviceSpaceType> info("info", N * VectorLength);

    cublasStatus_t stat;
    cublasHandle_t handle;

    stat = cublasCreate(&handle);
    if (stat != CUBLAS_STATUS_SUCCESS) Kokkos::abort("CUBLAS initialization failed\n");

    auto amat_device = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), amat);
    Kokkos::deep_copy(amat_device, amat);

    Kokkos::fence();
    {
      double tavg = 0, tmin = tmax;
      value_type *aa[N * VectorLength];

      for (int k = 0; k < N * VectorLength; ++k) {
        aa[k] = a.data() + k * a.stride_0();
      }
      value_type **aa_device;
      if (cudaMalloc(&aa_device, N * VectorLength * sizeof(value_type *)) != cudaSuccess) {
        Kokkos::abort("CUDA memory allocation failed\n");
      }
      if (cudaMemcpy(aa_device, aa, sizeof(value_type *) * N * VectorLength, cudaMemcpyHostToDevice) != cudaSuccess) {
        Kokkos::abort("CUDA memcpy failed\n");
      }
      Kokkos::fence();
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrix
        Kokkos::deep_copy(a, amat_device);

        Kokkos::fence();
        timer.reset();

        stat = cublasDgetrfBatched(handle, BlkSize, (value_type **)aa_device, BlkSize, NULL, (int *)info.data(),
                                   N * VectorLength);
        if (stat != CUBLAS_STATUS_SUCCESS) {
          Kokkos::abort("CUBLAS LU Batched failed\n");
        }

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto asol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a);
      Kokkos::deep_copy(asol, a);
      Kokkos::deep_copy(aref, asol);

      if (cudaFree(aa_device) != cudaSuccess) {
        Kokkos::abort("CUDA memory free failed\n");
      }

      std::cout << std::setw(8) << "CUBLAS" << std::setw(8) << "Batch"
                << " BlkSize = " << std::setw(3) << BlkSize << " TeamSize = N/A"
                << " ScratchSize (KB) = N/A"
                << " time = " << std::scientific << tmin << " avg flop/s = " << (flop / tavg)
                << " max flop/s = " << (flop / tmin) << std::endl;
    }
  }
#endif

  if (1) {
    ///
    /// Range policy version
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Functor<view_type, AlgoTagType> functor_type;
      const Kokkos::RangePolicy<DeviceSpaceType, ScheduleType, RangeTag> policy(0, N * VectorLength);

      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrix
        Kokkos::deep_copy(a, amat);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::LUCuda::RangeTag", policy, functor_type(a));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto asol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a);
      Kokkos::deep_copy(asol, a);

      double diff = 0;
      for (int i = 0, iend = aref.extent(0); i < iend; ++i)
        for (int j = 0, jend = aref.extent(1); j < jend; ++j)
          for (int k = 0, kend = aref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(aref(i, j, k) - asol(i, j, k));

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
    /// Team V1
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV1> policy_type;
      typedef Functor<view_type, AlgoTagType, VectorLength> functor_type;

      const int team_size = policy_type(N / 32, Kokkos::AUTO, VectorLength)
                                .team_size_recommended(functor_type(), Kokkos::ParallelForTag());

      const policy_type policy(N / team_size, team_size, VectorLength);
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrix
        Kokkos::deep_copy(a, amat);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::LUCuda::TeamTagV1", policy, functor_type(a));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto asol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a);
      Kokkos::deep_copy(asol, a);

      double diff = 0;
      for (int i = 0, iend = aref.extent(0); i < iend; ++i)
        for (int j = 0, jend = aref.extent(1); j < jend; ++j)
          for (int k = 0, kend = aref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(aref(i, j, k) - asol(i, j, k));

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
    /// Team V2
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV2> policy_type;
      typedef Functor<view_type, AlgoTagType, VectorLength> functor_type;

      const int is_blocked_algo = (std::is_same<AlgoTagType, Algo::LU::Blocked>::value),
                mb              = Algo::LU::Blocked::mb<DeviceMemorySpaceType>();
      // mp = BlkSize%mb > 0;

      const int
          // mblk = is_blocked_algo ? (BlkSize/mb + mp) : BlkSize;
          mblk = is_blocked_algo ? (BlkSize - mb) : (BlkSize - 1);

      const int max_team_size =
          policy_type(N, Kokkos::AUTO, VectorLength).team_size_max(functor_type(), Kokkos::ParallelForTag());
      const int team_size = std::min(std::max(mblk * 2, 1), max_team_size);

      const policy_type policy(N, team_size, VectorLength);
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrix
        Kokkos::deep_copy(a, amat);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::LUCuda::TeamTagV2", policy, functor_type(a));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto asol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a);
      Kokkos::deep_copy(asol, a);

      double diff = 0;
      for (int i = 0, iend = aref.extent(0); i < iend; ++i)
        for (int j = 0, jend = aref.extent(1); j < jend; ++j)
          for (int k = 0, kend = aref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(aref(i, j, k) - asol(i, j, k));

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
    /// Team V3
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV3> policy_type;
      typedef Functor<view_type, AlgoTagType, VectorLength> functor_type;

      const int lvl = 0, per_team_scratch = ScratchViewType<view_type>::shmem_size(VectorLength, BlkSize, BlkSize);
      if (per_team_scratch / 1024 < 48) {
        const int is_blocked_algo = (std::is_same<AlgoTagType, Algo::LU::Blocked>::value),
                  mb              = Algo::LU::Blocked::mb<DeviceMemorySpaceType>();
        //                  mp = BlkSize%mb > 0;

        const int
            // mblk = is_blocked_algo ? (BlkSize/mb + mp) : BlkSize;
            mblk = is_blocked_algo ? (BlkSize - mb) : (BlkSize - 1);

        const int max_team_size = policy_type(N, Kokkos::AUTO, VectorLength)
                                      .set_scratch_size(lvl, Kokkos::PerTeam(per_team_scratch))
                                      .team_size_max(functor_type(), Kokkos::ParallelForTag());
        const int team_size = std::min(std::max(mblk * 2, 1), max_team_size);

        policy_type policy(N, team_size, VectorLength);
        for (int iter = iter_begin; iter < iter_end; ++iter) {
          // flush
          flush.run();

          // initialize matrix
          Kokkos::deep_copy(a, amat);

          Kokkos::fence();
          timer.reset();

          Kokkos::parallel_for("KokkosBatched::PerfTest::LUCuda::TeamTagV3",
                               policy.set_scratch_size(lvl, Kokkos::PerTeam(per_team_scratch)), functor_type(a));

          Kokkos::fence();
          const double t = timer.seconds();
          tmin           = std::min(tmin, t);
          tavg += (iter >= 0) * t;
        }
        tavg /= iter_end;

        auto asol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), a);
        Kokkos::deep_copy(asol, a);

        double diff = 0;
        for (int i = 0, iend = aref.extent(0); i < iend; ++i)
          for (int j = 0, jend = aref.extent(1); j < jend; ++j)
            for (int k = 0, kend = aref.extent(2); k < kend; ++k)
              diff += Kokkos::ArithTraits<value_type>::abs(aref(i, j, k) - asol(i, j, k));

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
                  << " Scratch per team is too big (KB): " << (per_team_scratch / 1024) << std::endl;
      }
    }
  }
  std::cout << "\n\n";
}
}  // namespace PerfTest
}  // namespace KokkosBatched

using namespace KokkosBatched;

template <typename AlgoTagType>
void run(const int N, const int B) {
  typedef Kokkos::DefaultExecutionSpace ExecSpace;

  Kokkos::print_configuration(std::cout, false);

  if (B != 0) {
    PerfTest::LU<ExecSpace, AlgoTagType>(N, B);
  } else {
    PerfTest::LU<ExecSpace, AlgoTagType>(N, 3);
    PerfTest::LU<ExecSpace, AlgoTagType>(N, 5);
    PerfTest::LU<ExecSpace, AlgoTagType>(N, 10);
    PerfTest::LU<ExecSpace, AlgoTagType>(N, 15);
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

    std::cout << "\n Testing Algo::LU::Unblocked\n";
    run<Algo::LU::Unblocked>(N, B);

    std::cout << "\n Testing LayoutLeft Algo::LU::Blocked\n";
    run<Algo::LU::Blocked>(N, B);
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
