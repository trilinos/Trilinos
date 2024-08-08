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

#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"
#include "KokkosBatched_Trsm_Team_Impl.hpp"

namespace KokkosBatched {
namespace PerfTest {

#undef FLOP_MUL
#undef FLOP_ADD
#define FLOP_MUL 1.0
#define FLOP_ADD 1.0
typedef double value_type;

double FlopCountLower(int mm, int nn) {
  double m = (double)mm;
  double n = (double)nn;
  return (FLOP_MUL * (0.5 * m * n * (n + 1.0)) + FLOP_ADD * (0.5 * m * n * (n - 1.0)));
}

double FlopCountUpper(int mm, int nn) {
  double m = (double)mm;
  double n = (double)nn;
  return (FLOP_MUL * (0.5 * m * n * (n + 1.0)) + FLOP_ADD * (0.5 * m * n * (n - 1.0)));
}

struct RangeTag {};
struct TeamTagV1 {};
struct TeamTagV2 {};
struct TeamTagV3 {};
struct TeamTagHandmade {};

template <int test, typename ViewType, typename AlgoTagType, int VectorLength = 0>
struct Functor {
  ConstUnmanagedViewType<ViewType> _a;
  UnmanagedViewType<ViewType> _b;

  KOKKOS_INLINE_FUNCTION
  Functor() = default;

  KOKKOS_INLINE_FUNCTION
  Functor(const ViewType &a, const ViewType &b) : _a(a), _b(b) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const RangeTag &, const int k) const {
    auto aa = Kokkos::subview(_a, k, Kokkos::ALL(), Kokkos::ALL());
    auto bb = Kokkos::subview(_b, k, Kokkos::ALL(), Kokkos::ALL());

    switch (test) {
      case 0:
        SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(1.0, aa, bb);
        break;
      case 1:
        SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(1.0, aa, bb);
        break;
      case 2:
        SerialTrsm<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(1.0, aa, bb);
        break;
      case 3:
        SerialTrsm<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(1.0, aa, bb);
        break;
      case 4:
        SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(1.0, aa, bb);
        break;
    }
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV1 &, const MemberType &member) const {
    const int kbeg = (member.league_rank() * (member.team_size() * VectorLength) + member.team_rank() * VectorLength);
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < int(_b.extent(0))) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        auto bb = Kokkos::subview(_b, kk, Kokkos::ALL(), Kokkos::ALL());

        switch (test) {
          case 0:
            SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(1.0, aa, bb);
            break;
          case 1:
            SerialTrsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(1.0, aa, bb);
            break;
          case 2:
            SerialTrsm<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(1.0, aa, bb);
            break;
          case 3:
            SerialTrsm<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(1.0, aa, bb);
            break;
          case 4:
            SerialTrsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(1.0, aa, bb);
            break;
        }
      }
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV2 &, const MemberType &member) const {
    const int kbeg = member.league_rank() * VectorLength;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < int(_b.extent(0))) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        auto bb = Kokkos::subview(_b, kk, Kokkos::ALL(), Kokkos::ALL());

        switch (test) {
          case 0:
            TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(
                member, 1.0, aa, bb);
            break;
          case 1:
            TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(
                member, 1.0, aa, bb);
            break;
          case 2:
            TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(
                member, 1.0, aa, bb);
            break;
          case 3:
            TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(
                member, 1.0, aa, bb);
            break;
          case 4:
            TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(
                member, 1.0, aa, bb);
            break;
        }
      }
    });
  }

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const TeamTagV3 &, const MemberType &member) const {
    const int lvl = 0;
    ScratchViewType<ViewType> sa(member.team_scratch(lvl), VectorLength, _a.extent(1), _a.extent(2));
    // ScratchViewType<ViewType> sb(member.team_scratch(lvl), VectorLength,
    // _b.extent(1), _b.extent(2));

    const int kbeg = member.league_rank() * VectorLength;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength), [&](const int &k) {
      const int kk = kbeg + k;
      if (kk < int(_b.extent(0))) {
        auto aa = Kokkos::subview(_a, kk, Kokkos::ALL(), Kokkos::ALL());
        auto bb = Kokkos::subview(_b, kk, Kokkos::ALL(), Kokkos::ALL());

        auto saa = Kokkos::subview(sa, k, Kokkos::ALL(), Kokkos::ALL());

        TeamCopy<MemberType, Trans::NoTranspose>::invoke(member, aa, saa);
        member.team_barrier();

        switch (test) {
          case 0:
            TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(
                member, 1.0, saa, bb);
            break;
          case 1:
            TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(
                member, 1.0, saa, bb);
            break;
          case 2:
            TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit, AlgoTagType>::invoke(
                member, 1.0, saa, bb);
            break;
          case 3:
            TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(
                member, 1.0, saa, bb);
            break;
          case 4:
            TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, AlgoTagType>::invoke(
                member, 1.0, saa, bb);
            break;
        }
      }
    });
  }
};

template <int test, typename DeviceSpaceType, typename AlgoTagType>
void Trsm(const int NN, const int BlkSize, const int NumCols) {
  typedef Kokkos::Schedule<Kokkos::Static> ScheduleType;

  constexpr int VectorLength = DefaultVectorLength<value_type, typename DeviceSpaceType::memory_space>::value;
  const int N                = NN / VectorLength;

  {
    std::string value_type_name;
    if (std::is_same<value_type, double>::value) value_type_name = "double";
    if (std::is_same<value_type, Kokkos::complex<double> >::value) value_type_name = "Kokkos::complex<double>";

    std::cout << "SIMD is defined: datatype " << value_type_name << " a vector length " << VectorLength << "\n";
  }

  switch (test) {
    case 0: std::cout << "TestID = Left,  Lower, NoTrans,    UnitDiag\n"; break;
    case 1: std::cout << "TestID = Left,  Lower, NoTrans, NonUnitDiag\n"; break;
    case 2: std::cout << "TestID = Right, Upper, NoTrans,    UnitDiag\n"; break;
    case 3: std::cout << "TestID = Right, Upper, NoTrans, NonUnitDiag\n"; break;
    case 4: std::cout << "TestID = Left,  Upper, NoTrans, NonUnitDiag\n"; break;
  }

  // when m == n, lower upper does not matter (unit and nonunit)
  double flop = 0;
  switch (test) {
    case 0:
    case 1: flop = FlopCountLower(BlkSize, NumCols); break;
    case 2:
    case 3:
    case 4: flop = FlopCountUpper(BlkSize, NumCols); break;
  }
  flop *= (N * VectorLength);
  const double tmax = 1.0e15;

  typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;
  typedef typename DeviceSpaceType::memory_space DeviceMemorySpaceType;

  const int iter_begin = -3, iter_end = 30;
  Kokkos::Timer timer;

  Kokkos::View<value_type ***, Kokkos::LayoutLeft, HostSpaceType> amat("amat", N * VectorLength, BlkSize, BlkSize),
      bmat("bmat", N * VectorLength, BlkSize, NumCols), bref("bmat", N * VectorLength, BlkSize, NumCols);

  {
    Random<value_type> random;
    for (int k = 0; k < N * VectorLength; ++k) {
      for (int i = 0; i < BlkSize; ++i)
        for (int j = 0; j < BlkSize; ++j) amat(k, i, j) = random.value() + 4.0 * (i == j);
      for (int i = 0; i < BlkSize; ++i)
        for (int j = 0; j < NumCols; ++j) bmat(k, i, j) = random.value();
    }
  }

  // P100 L2 cache 4MB per core
  constexpr size_t LLC_CAPACITY = 56 * 4 * 1024 * 1024;
  Flush<LLC_CAPACITY, DeviceSpaceType> flush;

#if defined(__KOKKOSKERNELS_NVIDIA_CUBLAS__)
  if (1) {
    ///
    /// CUBLAS Batch version
    ///
    const Kokkos::LayoutStride stride(N * VectorLength, BlkSize * BlkSize, BlkSize, 1, BlkSize, BlkSize);

    Kokkos::View<value_type ***, Kokkos::LayoutStride, DeviceSpaceType> a("a", stride), b("b", stride);

    cublasStatus_t stat;
    cublasHandle_t handle;

    stat = cublasCreate(&handle);
    if (stat != CUBLAS_STATUS_SUCCESS) Kokkos::abort("CUBLAS initialization failed\n");

    auto amat_device = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), amat);
    auto bmat_device = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), bmat);

    Kokkos::deep_copy(amat_device, amat);
    Kokkos::deep_copy(bmat_device, bmat);

    Kokkos::fence();

    const double one(1.0);  //, zero(0.0);
    {
      double tavg = 0, tmin = tmax;
      value_type *aa[N * VectorLength], *bb[N * VectorLength];
      for (int k = 0; k < N * VectorLength; ++k) {
        aa[k] = a.data() + k * a.stride_0();
        bb[k] = b.data() + k * b.stride_0();
      }
      value_type **aa_device, **bb_device;
      if (cudaMalloc(&aa_device, N * VectorLength * sizeof(value_type *)) != cudaSuccess ||
          cudaMalloc(&bb_device, N * VectorLength * sizeof(value_type *)) != cudaSuccess) {
        Kokkos::abort("CUDA memory allocation failed\n");
      }
      if (cudaMemcpy(aa_device, aa, sizeof(value_type *) * N * VectorLength, cudaMemcpyHostToDevice) != cudaSuccess ||
          cudaMemcpy(bb_device, bb, sizeof(value_type *) * N * VectorLength, cudaMemcpyHostToDevice) != cudaSuccess) {
        Kokkos::abort("CUDA memcpy failed\n");
      }
      Kokkos::fence();
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat_device);
        Kokkos::deep_copy(b, bmat_device);

        timer.reset();
        switch (test) {
          case 0: {
            // Left,  Lower, NoTrans,    UnitDiag
            stat = cublasDtrsmBatched(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT,
                                      BlkSize, NumCols, &one, (const value_type **)aa_device, BlkSize,
                                      (value_type **)bb_device, BlkSize, N * VectorLength);
            break;
          }
          case 1: {
            // Left,  Lower, NoTrans, NonUnitDiag
            stat = cublasDtrsmBatched(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,
                                      CUBLAS_DIAG_NON_UNIT, BlkSize, NumCols, &one, (const value_type **)aa_device,
                                      BlkSize, (value_type **)bb_device, BlkSize, N * VectorLength);
            break;
          }
          case 2: {
            // Right, Upper, NoTrans,    UnitDiag
            stat = cublasDtrsmBatched(handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT,
                                      BlkSize, NumCols, &one, (const value_type **)aa_device, BlkSize,
                                      (value_type **)bb_device, BlkSize, N * VectorLength);
            break;
          }
          case 3: {
            // Right, Upper, NoTrans, NonUnitDiag
            stat = cublasDtrsmBatched(handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N,
                                      CUBLAS_DIAG_NON_UNIT, BlkSize, NumCols, &one, (const value_type **)aa_device,
                                      BlkSize, (value_type **)bb_device, BlkSize, N * VectorLength);
            break;
          }
          case 4: {
            // Left,  Upper, NoTrans, NonUnitDiag
            stat = cublasDtrsmBatched(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N,
                                      CUBLAS_DIAG_NON_UNIT, BlkSize, NumCols, &one, (const value_type **)aa_device,
                                      BlkSize, (value_type **)bb_device, BlkSize, N * VectorLength);
            break;
          }
        }

        if (stat != CUBLAS_STATUS_SUCCESS) {
          Kokkos::abort("CUBLAS Trsm Batched failed\n");
        }
        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto bsol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b);
      Kokkos::deep_copy(bsol, b);
      Kokkos::deep_copy(bref, bsol);

      if (cudaFree(aa_device) != cudaSuccess || cudaFree(bb_device) != cudaSuccess) {
        Kokkos::abort("CUDA memory free failed\n");
      }

      std::cout << std::setw(8) << "CUBLAS" << std::setw(8) << "Batched"
                << " BlkSize = " << std::setw(3) << BlkSize << " NumCols = " << std::setw(3) << NumCols
                << " TeamSize = N/A"
                << " ScratchSize (KB) = N/A"
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
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, NumCols);

    double tavg = 0, tmin = tmax;
    {
      typedef Functor<test, view_type, AlgoTagType> functor_type;
      const Kokkos::RangePolicy<DeviceSpaceType, ScheduleType, RangeTag> policy(0, N * VectorLength);

      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::RangeTag", policy, functor_type(a, b));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto bsol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b);
      Kokkos::deep_copy(bsol, b);

      double diff = 0;
      for (int i = 0, iend = bref.extent(0); i < iend; ++i)
        for (int j = 0, jend = bref.extent(1); j < jend; ++j)
          for (int k = 0, kend = bref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(bref(i, j, k) - bsol(i, j, k));

      std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Range"
                << " BlkSize = " << std::setw(3) << BlkSize << " NumCols = " << std::setw(3) << NumCols
                << " TeamSize = N/A"
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
    /// Team policy V1 - almost same scheduling with range policy
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, NumCols);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV1> policy_type;
      typedef Functor<test, view_type, AlgoTagType, VectorLength> functor_type;

      const int team_size = policy_type(N / 32, Kokkos::AUTO, VectorLength)
                                .team_size_recommended(functor_type(), Kokkos::ParallelForTag());

      const policy_type policy(N / team_size, team_size, VectorLength);
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);

        Kokkos::fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::TeamTagV1", policy, functor_type(a, b));

        Kokkos::fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto bsol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b);
      Kokkos::deep_copy(bsol, b);

      double diff = 0;
      for (int i = 0, iend = bref.extent(0); i < iend; ++i)
        for (int j = 0, jend = bref.extent(1); j < jend; ++j)
          for (int k = 0, kend = bref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(bref(i, j, k) - bsol(i, j, k));

      std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team V1"
                << " BlkSize = " << std::setw(3) << BlkSize << " NumCols = " << std::setw(3) << NumCols
                << " TeamSize = " << std::setw(3) << team_size << " ScratchSize (KB) =   0"
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
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, NumCols);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV2> policy_type;
      typedef Functor<test, view_type, AlgoTagType, VectorLength> functor_type;

      const int is_blocked_algo = (std::is_same<AlgoTagType, Algo::Trsm::Blocked>::value),
                mb = Algo::Trsm::Blocked::mb<DeviceMemorySpaceType>(), mp = BlkSize % mb > 0;

      const int mblk = is_blocked_algo ? (BlkSize / mb + mp) : BlkSize;

      const int max_team_size =
          policy_type(N, Kokkos::AUTO, VectorLength).team_size_max(functor_type(), Kokkos::ParallelForTag());
      const int team_size = std::min(std::max(NumCols, (mblk - 1) * mblk), max_team_size);

      const policy_type policy(N, team_size, VectorLength);
      for (int iter = iter_begin; iter < iter_end; ++iter) {
        // flush
        flush.run();

        // initialize matrices
        Kokkos::deep_copy(a, amat);
        Kokkos::deep_copy(b, bmat);

        DeviceSpaceType().fence();
        timer.reset();

        Kokkos::parallel_for("KokkosBatched::PerfTest::TeamTagV2", policy, functor_type(a, b));

        DeviceSpaceType().fence();
        const double t = timer.seconds();
        tmin           = std::min(tmin, t);
        tavg += (iter >= 0) * t;
      }
      tavg /= iter_end;

      auto bsol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b);
      Kokkos::deep_copy(bsol, b);

      double diff = 0;
      for (int i = 0, iend = bref.extent(0); i < iend; ++i)
        for (int j = 0, jend = bref.extent(1); j < jend; ++j)
          for (int k = 0, kend = bref.extent(2); k < kend; ++k)
            diff += Kokkos::ArithTraits<value_type>::abs(bref(i, j, k) - bsol(i, j, k));

      std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team V2"
                << " BlkSize = " << std::setw(3) << BlkSize << " NumCols = " << std::setw(3) << NumCols
                << " TeamSize = " << std::setw(3) << team_size << " ScratchSize (KB) =   0"
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
    /// Team policy V3 - team parallel + sratch
    ///
    typedef Kokkos::View<value_type ***, DeviceSpaceType> view_type;
    view_type a("a", N * VectorLength, BlkSize, BlkSize), b("b", N * VectorLength, BlkSize, NumCols);

    double tavg = 0, tmin = tmax;
    {
      typedef Kokkos::TeamPolicy<DeviceSpaceType, ScheduleType, TeamTagV3> policy_type;
      typedef Functor<test, view_type, AlgoTagType, VectorLength> functor_type;

      const int lvl = 0, per_team_scratch = ScratchViewType<view_type>::shmem_size(VectorLength, BlkSize, BlkSize);

      if (per_team_scratch / 1024 < 48) {
        const int is_blocked_algo = (std::is_same<AlgoTagType, Algo::Trsm::Blocked>::value),
                  mb = Algo::Trsm::Blocked::mb<DeviceMemorySpaceType>(), mp = BlkSize % mb > 0;

        const int mblk = is_blocked_algo ? (BlkSize / mb + mp) : BlkSize;

        const int max_team_size = policy_type(N, Kokkos::AUTO, VectorLength)
                                      .set_scratch_size(lvl, Kokkos::PerTeam(per_team_scratch))
                                      .team_size_max(functor_type(), Kokkos::ParallelForTag());
        const int team_size = std::min(std::max(NumCols, (mblk - 1) * mblk), max_team_size);

        policy_type policy(N, team_size, VectorLength);
        for (int iter = iter_begin; iter < iter_end; ++iter) {
          // flush
          flush.run();

          // initialize matrices
          Kokkos::deep_copy(a, amat);
          Kokkos::deep_copy(b, bmat);

          DeviceSpaceType().fence();
          timer.reset();

          Kokkos::parallel_for("KokkosBatched::PerfTest::TeamTagV3", policy, functor_type(a, b));

          DeviceSpaceType().fence();
          const double t = timer.seconds();
          tmin           = std::min(tmin, t);
          tavg += (iter >= 0) * t;
        }
        tavg /= iter_end;

        auto bsol = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), b);
        Kokkos::deep_copy(bsol, b);

        double diff = 0;
        for (int i = 0, iend = bref.extent(0); i < iend; ++i)
          for (int j = 0, jend = bref.extent(1); j < jend; ++j)
            for (int k = 0, kend = bref.extent(2); k < kend; ++k)
              diff += Kokkos::ArithTraits<value_type>::abs(bref(i, j, k) - bsol(i, j, k));

        std::cout << std::setw(8) << "Kokkos" << std::setw(8) << "Team V3"
                  << " BlkSize = " << std::setw(3) << BlkSize << " NumCols = " << std::setw(3) << NumCols
                  << " TeamSize = " << std::setw(3) << team_size << " ScratchSize (KB) = " << std::setw(3)
                  << (per_team_scratch / 1024) << " time = " << std::scientific << tmin
                  << " avg flop/s = " << (flop / tavg) << " max flop/s = " << (flop / tmin);
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
void run(const int N, const int B, const int R) {
  typedef Kokkos::DefaultExecutionSpace ExecSpace;

  Kokkos::print_configuration(std::cout, false);

  if (B != 0 && R != 0) {
    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, B, R);
  } else {
    std::cout << "\n\n Used for Factorization \n\n";

    /// Left, Lower, NoTrans, UnitDiag (used in LU factorization and LU solve)

    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 3, 3);
    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 5, 5);
    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 10, 10);
    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 15, 15);

    /// Left, Lower, NoTrans, NonUnitDiag

    PerfTest::Trsm<1, ExecSpace, AlgoTagType>(N, 3, 3);
    PerfTest::Trsm<1, ExecSpace, AlgoTagType>(N, 5, 5);
    PerfTest::Trsm<1, ExecSpace, AlgoTagType>(N, 10, 10);
    PerfTest::Trsm<1, ExecSpace, AlgoTagType>(N, 15, 15);

    /// Right, Upper, NoTrans, UnitDiag

    PerfTest::Trsm<2, ExecSpace, AlgoTagType>(N, 3, 3);
    PerfTest::Trsm<2, ExecSpace, AlgoTagType>(N, 5, 5);
    PerfTest::Trsm<2, ExecSpace, AlgoTagType>(N, 10, 10);
    PerfTest::Trsm<2, ExecSpace, AlgoTagType>(N, 15, 15);

    /// Right, Upper, NoTrans, NonUnitDiag (used in LU factorization)

    PerfTest::Trsm<3, ExecSpace, AlgoTagType>(N, 3, 3);
    PerfTest::Trsm<3, ExecSpace, AlgoTagType>(N, 5, 5);
    PerfTest::Trsm<3, ExecSpace, AlgoTagType>(N, 10, 10);
    PerfTest::Trsm<3, ExecSpace, AlgoTagType>(N, 15, 15);

    std::cout << "\n\n Used for Solve \n\n";

    /// Left, Lower, NoTrans, UnitDiag (used in LU solve)

    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 3, 1);
    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 5, 1);
    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 10, 1);
    PerfTest::Trsm<0, ExecSpace, AlgoTagType>(N, 15, 1);

    /// Left, Upper, Notrans, NonUnitDiag (user in LU solve)

    PerfTest::Trsm<4, ExecSpace, AlgoTagType>(N, 3, 1);
    PerfTest::Trsm<4, ExecSpace, AlgoTagType>(N, 5, 1);
    PerfTest::Trsm<4, ExecSpace, AlgoTagType>(N, 10, 1);
    PerfTest::Trsm<4, ExecSpace, AlgoTagType>(N, 15, 1);
  }
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);

  int N = 128 * 128, B = 0, R = 0;

  for (int i = 1; i < argc; ++i) {
    const std::string &token = argv[i];
    if (token == std::string("-N")) N = std::atoi(argv[++i]);
    if (token == std::string("-B")) B = std::atoi(argv[++i]);
    if (token == std::string("-R")) R = std::atoi(argv[++i]);
  }

  if (R == 0 && B != 0) R = B;

  {
    std::cout << " N = " << N << std::endl;

    std::cout << "\n Testing LayoutLeft and Algo::Trsm::Unblocked\n";
    run<Algo::Trsm::Unblocked>(N, B, R);

    std::cout << "\n Testing LayoutLeft Algo::Trsm::Blocked\n";
    run<Algo::Trsm::Blocked>(N, B, R);
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
