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
/// Kokkos headers
#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"
#include "Kokkos_Random.hpp"

#if !(defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_CUDA_LAMBDA))
#define KOKKOSBATCHED_TEST_BLOCKTRIDIAGJACOBI
#endif

#if defined(KOKKOSBATCHED_TEST_BLOCKTRIDIAGJACOBI)

/// KokkosKernels headers
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Copy_Impl.hpp>
#include <KokkosBatched_SetIdentity_Decl.hpp>
#include <KokkosBatched_SetIdentity_Impl.hpp>
#include <KokkosBatched_AddRadial_Decl.hpp>
#include <KokkosBatched_AddRadial_Impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemm_Team_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Gemv_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Team_Impl.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Trsm_Serial_Impl.hpp>
#include <KokkosBatched_Trsm_Team_Impl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_Trsv_Serial_Impl.hpp>
#include <KokkosBatched_Trsv_Team_Impl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>

#define KOKKOSBATCHED_PROFILE 1
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
#include "cuda_profiler_api.h"
#endif

#define KOKKOSBATCHED_USE_128BIT_MEMORY_INST

typedef Kokkos::DefaultExecutionSpace exec_space;
typedef typename exec_space::memory_space memory_space;
typedef Kokkos::DefaultHostExecutionSpace host_space;

typedef double value_type;

/// 128*128*128/16*5 * (2*8) / 16
///
/// simd typedefs
///
using namespace KokkosBatched;

static constexpr int vector_length = DefaultVectorLength<value_type, memory_space>::value;
#if defined(KOKKOSBATCHED_USE_128BIT_MEMORY_INST)
static constexpr int internal_vector_length = DefaultInternalVectorLength<value_type, memory_space>::value;
#else
static constexpr int internal_vector_length = 1;
#endif

typedef Vector<SIMD<value_type>, vector_length> vector_type;
#if defined(KOKKOSBATCHED_USE_128BIT_MEMORY_INST)
typedef Vector<SIMD<value_type>, internal_vector_length> internal_vector_type;
#else
typedef value_type internal_vector_type;
#endif

template <typename ExecutionSpace>
struct InverseDiagonalsModeAndAlgo;

struct InverseDiagonalsModeAndAlgoHostImpl {
  typedef Mode::Serial mode_type;
  typedef Algo::Level3::Blocked algo_type;
};

#if defined(KOKKOS_ENABLE_SERIAL)
template <>
struct InverseDiagonalsModeAndAlgo<Kokkos::Serial> : InverseDiagonalsModeAndAlgoHostImpl {};
#endif

#if defined(KOKKOS_ENABLE_THREADS)
template <>
struct InverseDiagonalsModeAndAlgo<Kokkos::Threads> : InverseDiagonalsModeAndAlgoHostImpl {};
#endif

#if defined(KOKKOS_ENABLE_ONPENMP)
template <>
struct InverseDiagonalsModeAndAlgo<Kokkos::Threads> : InverseDiagonalsModeAndAlgoHostImpl {};
#endif

struct InverseDiagonalsModeAndAlgoDeviceImpl {
  typedef Mode::Team mode_type;
  typedef Algo::Level3::Unblocked algo_type;
};

#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct InverseDiagonalsModeAndAlgo<Kokkos::Cuda> : InverseDiagonalsModeAndAlgoDeviceImpl {};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <>
struct InverseDiagonalsModeAndAlgo<Kokkos::HIP> : InverseDiagonalsModeAndAlgoDeviceImpl {};
#endif

template <typename ExecutionSpace>
struct SolveModeAndAlgo;

struct SolveModeAndAlgoHostImpl {
  typedef Mode::Serial mode_type;
  typedef Algo::Level2::Blocked algo_type;
};

#if defined(KOKKOS_ENABLE_SERIAL)
template <>
struct SolveModeAndAlgo<Kokkos::Serial> : SolveModeAndAlgoHostImpl {};
#endif

#if defined(KOKKOS_ENABLE_THREADS)
template <>
struct SolveModeAndAlgo<Kokkos::Threads> : SolveModeAndAlgoHostImpl {};
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
template <>
struct SolveModeAndAlgo<Kokkos::OpenMP> : SolveModeAndAlgoHostImpl {};
#endif

struct SolveModeAndAlgoDeviceImpl {
  typedef Mode::Team mode_type;
  typedef Algo::Level2::Unblocked algo_type;
};

#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct SolveModeAndAlgo<Kokkos::Cuda> : SolveModeAndAlgoDeviceImpl {};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <>
struct SolveModeAndAlgo<Kokkos::HIP> : SolveModeAndAlgoDeviceImpl {};
#endif

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
    cudaProfilerStop();
#endif
    Kokkos::print_configuration(std::cout);

    // typedef Kokkos::ArithTraits<value_type> ats;
    Kokkos::Timer timer;

    ///
    /// input arguments parsing
    ///
    int N      = 128 * 128;  /// # of problems (batch size)
    int L      = 128;        /// length of block tridiags
    int Blk    = 5;          /// block dimension
    int Nvec   = 1;
    int S      = 0;  /// scratch size
    int niter  = 1;
    int nsweep = 10;
    for (int i = 1; i < argc; ++i) {
      const std::string &token = argv[i];
      if (token == std::string("-N")) N = std::atoi(argv[++i]);
      if (token == std::string("-L")) L = std::atoi(argv[++i]);
      if (token == std::string("-B")) Blk = std::atoi(argv[++i]);
      if (token == std::string("-Nvec")) Nvec = std::atoi(argv[++i]);
      if (token == std::string("-S")) S = std::atoi(argv[++i]);
      if (token == std::string("-Niter")) niter = std::atoi(argv[++i]);
      if (token == std::string("-Nsweep")) nsweep = std::atoi(argv[++i]);
    }

    printf(
        " :::: Testing (N = %d, L = %d, Blk = %d, vl = %d, vi = %d, niter = "
        "%d, nsweep = %d)\n",
        N, L, Blk, vector_length, internal_vector_length, niter, nsweep);

    ///
    /// problem container
    ///

    /// double 16
    Kokkos::View<vector_type *****, Kokkos::LayoutRight, exec_space> Av("A", N / vector_length, L, 4, Blk, Blk);

    /// double
    Kokkos::View<value_type ******, Kokkos::LayoutRight, exec_space> As(
        (value_type *)Av.data(), Av.extent(0), Av.extent(1), Av.extent(2), Av.extent(3), Av.extent(4), vector_length);

    /// double 2
    Kokkos::View<internal_vector_type ******, Kokkos::LayoutRight, exec_space> Ai(
        (internal_vector_type *)Av.data(), Av.extent(0), Av.extent(1), Av.extent(2), Av.extent(3), Av.extent(4),
        vector_length / internal_vector_length);
    /// double 16
    Kokkos::View<vector_type *****, Kokkos::LayoutRight, exec_space> xv("x", N / vector_length, Nvec, 2, L, Blk);

    /// double
    Kokkos::View<value_type ******, Kokkos::LayoutRight, exec_space> xs(
        (value_type *)xv.data(), xv.extent(0), xv.extent(1), xv.extent(2), xv.extent(3), xv.extent(4), vector_length);

    /// double 2
    Kokkos::View<internal_vector_type ******, Kokkos::LayoutRight, exec_space> xi(
        (internal_vector_type *)xv.data(), xv.extent(0), xv.extent(1), xv.extent(2), xv.extent(3), xv.extent(4),
        vector_length / internal_vector_length);

    /// double 16
    Kokkos::View<vector_type ****, Kokkos::LayoutRight, exec_space> bv("b", N / vector_length, Nvec, L, Blk);

    /// double
    Kokkos::View<value_type *****, Kokkos::LayoutRight, exec_space> bs(
        (value_type *)bv.data(), bv.extent(0), bv.extent(1), bv.extent(2), bv.extent(3), vector_length);

    /// double 2
    Kokkos::View<internal_vector_type *****, Kokkos::LayoutRight, exec_space> bi(
        (internal_vector_type *)bv.data(), bv.extent(0), bv.extent(1), bv.extent(2), bv.extent(3),
        vector_length / internal_vector_length);

    /// double copy of A
    Kokkos::View<value_type ******, Kokkos::LayoutRight, exec_space> Acopy(
        "Acopy", As.extent(0), As.extent(1), As.extent(2), As.extent(3), As.extent(4), As.extent(5));

    Kokkos::View<value_type *****, Kokkos::LayoutRight, exec_space> rs("rs", bs.extent(0), bs.extent(1), bs.extent(2),
                                                                       bs.extent(3), bs.extent(4));

#if defined(KOKKOSBATCHED_USE_128BIT_MEMORY_INST)
    auto AA = Ai;
    auto bb = bi;
    auto xx = xi;
#else
    auto AA = As;
    auto bb = bs;
    auto xx = xs;
#endif

    /// randomize input
    Kokkos::Random_XorShift64_Pool<exec_space> random(13245);
    Kokkos::fill_random(As, random, value_type(1.0));
    Kokkos::fill_random(bs, random, value_type(1.0));

    ///
    /// diagonal dominant
    ///
    if (1) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStart();
#endif
      using policy_type = Kokkos::TeamPolicy<exec_space>;
      using member_type = typename policy_type::member_type;
      policy_type policy(AA.extent(0) * L, Kokkos::AUTO(), AA.extent(5));
      Kokkos::parallel_for(
          "diagonal dominant", policy, KOKKOS_LAMBDA(const member_type &member) {
            const int i = member.league_rank() / L;
            const int k = member.league_rank() % L;
            Kokkos::parallel_for(Kokkos::TeamThreadRange(member, Blk), [&](const int &j) {
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, AA.extent(5)),
                                   [&](const int &v) { AA(i, k, 1, j, j, v) += internal_vector_type(9 * Blk); });
            });
          });
      Kokkos::fence();
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStop();
#endif
    }

    Kokkos::deep_copy(Acopy, As);

    ///
    /// compute the inverse of diagonals
    ///
    if (1) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStart();
#endif
      timer.reset();
      typedef internal_vector_type scratch_value_type;
      typedef Kokkos::View<scratch_value_type ***, Kokkos::LayoutRight, typename exec_space::scratch_memory_space,
                           Kokkos::MemoryUnmanaged>
          scratch_view_type;

      using policy_type          = Kokkos::TeamPolicy<exec_space>;
      using member_type          = typename policy_type::member_type;
      const int per_team_scratch = scratch_view_type::shmem_size(Blk, Blk, AA.extent(5));
      int team_size              = 0;
      if (Blk < 8) {
        team_size = 32 / AA.extent(5);
      } else if (Blk < 12) {
        team_size = 32 / AA.extent(5);
      } else {
        team_size = 64 / AA.extent(5);
      }

      policy_type policy(AA.extent(0) * L, team_size, AA.extent(5));
      Kokkos::parallel_for(
          "inverse diagonals", policy.set_scratch_size(0, Kokkos::PerTeam(S < per_team_scratch ? per_team_scratch : S)),
          KOKKOS_LAMBDA(const member_type &member) {
            typedef InverseDiagonalsModeAndAlgo<Kokkos::DefaultExecutionSpace> default_mode_and_algo_type;
            typedef default_mode_and_algo_type::mode_type mode_type;
            typedef default_mode_and_algo_type::algo_type algo_type;

            const int i = member.league_rank() / L;
            const int k = member.league_rank() % L;

            scratch_view_type WW(member.team_scratch(0), Blk, Blk, AA.extent(5));
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, AA.extent(5)), [&](const int &v) {
              auto A = Kokkos::subview(AA, i, k, 1, Kokkos::ALL(), Kokkos::ALL(), v);
              auto D = Kokkos::subview(AA, i, k, 3, Kokkos::ALL(), Kokkos::ALL(), v);
              auto W = Kokkos::subview(WW, Kokkos::ALL(), Kokkos::ALL(), v);

              Copy<member_type, Trans::NoTranspose, mode_type>::invoke(member, A, W);
              SetIdentity<member_type, mode_type>::invoke(member, D);
              member.team_barrier();
              LU<member_type, mode_type, algo_type>::invoke(member, W);
              Trsm<member_type, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, mode_type, algo_type>::invoke(
                  member, 1.0, W, D);
              Trsm<member_type, Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, mode_type,
                   algo_type>::invoke(member, 1.0, W, D);
            });
          });
      Kokkos::fence();
      const double t = timer.seconds();
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStop();
#endif
      printf("inverse time = %f , # of inverse per min = %f \n", t, 1.0 / t * 60);
    }

    ///
    /// solve the matrix 20 times
    ///
    if (1) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStart();
#endif
      timer.reset();
      typedef internal_vector_type scratch_value_type;
      typedef Kokkos::View<scratch_value_type **, Kokkos::LayoutRight, typename exec_space::scratch_memory_space,
                           Kokkos::MemoryUnmanaged>
          scratch_view_type;
      const int per_team_scratch = scratch_view_type::shmem_size(Blk, AA.extent(5));

      using policy_type = Kokkos::TeamPolicy<exec_space>;
      using member_type = typename policy_type::member_type;
      int team_size     = 0;
      if (Blk < 8) {
        team_size = 32 / AA.extent(5);
      } else if (Blk < 12) {
        team_size = 32 / AA.extent(5);
      } else {
        team_size = 32 / AA.extent(5);
      }
      policy_type policy(AA.extent(0) * L, team_size, AA.extent(5));

      for (int iter = 0; iter < niter; ++iter) {
        auto xxx = Kokkos::subview(xx, Kokkos::ALL(), Kokkos::ALL(), 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto yyy = Kokkos::subview(xx, Kokkos::ALL(), Kokkos::ALL(), 1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        for (int nis = 0; nis < nsweep; ++nis) {
          Kokkos::parallel_for(
              "solve", policy.set_scratch_size(0, Kokkos::PerTeam(S < per_team_scratch ? per_team_scratch : S)),
              KOKKOS_LAMBDA(const member_type &member) {
                typedef SolveModeAndAlgo<Kokkos::DefaultExecutionSpace> default_mode_and_algo_type;
                typedef default_mode_and_algo_type::mode_type mode_type;
                typedef default_mode_and_algo_type::algo_type algo_type;

                scratch_view_type WW(member.team_scratch(0), Blk, AA.extent(5));
                const int i = member.league_rank() / L;  //%AA.extent(0);
                const int k = member.league_rank() % L;
                Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, AA.extent(5)), [&](const int &v) {
                  auto A = Kokkos::subview(AA, i, k, 1, Kokkos::ALL(), Kokkos::ALL(), v);
                  auto D = Kokkos::subview(AA, i, k, 3, Kokkos::ALL(), Kokkos::ALL(), v);
                  auto B = Kokkos::subview(AA, i, k, 2, Kokkos::ALL(), Kokkos::ALL(), v);
                  auto C = Kokkos::subview(AA, i, k ? k - 1 : 0, 0, Kokkos::ALL(), Kokkos::ALL(), v);
                  auto u = Kokkos::subview(WW, Kokkos::ALL(), v);
                  for (int jvec = 0; jvec < Nvec; ++jvec) {
                    auto x0 = Kokkos::subview(xxx, i, jvec, k == 0 ? 0 : k - 1, Kokkos::ALL(), v);
                    auto x1 = Kokkos::subview(xxx, i, jvec, k, Kokkos::ALL(), v);
                    auto x2 = Kokkos::subview(xxx, i, jvec, k == L - 1 ? 0 : k + 1, Kokkos::ALL(), v);
                    auto y1 = Kokkos::subview(yyy, i, jvec, k, Kokkos::ALL(), v);
                    auto b  = Kokkos::subview(bb, i, jvec, k, Kokkos::ALL(), v);

                    if (L == 1) {
                      Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, 1.0, D, b, 0.0, x1);
                    } else {
                      Copy<member_type, Trans::NoTranspose, mode_type>::invoke(member, b, u);
                      if (k == 0) {
                        Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, -1.0, B, x2, 1.0,
                                                                                            u);
                      } else if (k == L - 1) {
                        Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, -1.0, C, x0, 1.0,
                                                                                            u);
                      } else {
                        Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, -1.0, B, x2, 1.0,
                                                                                            u);
                        Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, -1.0, C, x0, 1.0,
                                                                                            u);
                      }
                      Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, 1.0, D, u, 0.0, y1);
                    }
                  }
                });
              });
          auto tmp = xxx;
          xxx      = yyy;
          yyy      = tmp;
        }
        Kokkos::fence();
      }
      const double t = timer.seconds();
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStop();
#endif
      printf("solve time = %f , # of solves per min = %f\n", t, 1.0 / t * 60 * niter);
    }

    ///
    /// compute residual
    ///
    if (1) {
      typedef KokkosBatched::Algo::Level2::Unblocked algo_type;
      using policy_type = Kokkos::TeamPolicy<exec_space>;
      policy_type policy(Acopy.extent(0), Kokkos::AUTO(), Acopy.extent(5));
      Kokkos::parallel_for(
          "compute residual", policy, KOKKOS_LAMBDA(const typename policy_type::member_type &member) {
            const int i = member.league_rank();
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, Acopy.extent(5)), [&](const int &v) {
              auto A = Kokkos::subview(Acopy, i, Kokkos::ALL(), 1, Kokkos::ALL(), Kokkos::ALL(), v);
              auto B = Kokkos::subview(Acopy, i, Kokkos::ALL(), 2, Kokkos::ALL(), Kokkos::ALL(), v);
              auto C = Kokkos::subview(Acopy, i, Kokkos::ALL(), 0, Kokkos::ALL(), Kokkos::ALL(), v);

              for (int jvec = 0, jvecend = rs.extent(1); jvec < jvecend; ++jvec) {
                auto x = Kokkos::subview(xs, i, jvec, nsweep % 2, Kokkos::ALL(), Kokkos::ALL(), v);
                auto b = Kokkos::subview(bs, i, jvec, Kokkos::ALL(), Kokkos::ALL(), v);
                auto r = Kokkos::subview(rs, i, jvec, Kokkos::ALL(), Kokkos::ALL(), v);

                if (L == 1) {
                  auto A0 = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());
                  auto x0 = Kokkos::subview(x, 0, Kokkos::ALL());
                  auto b0 = Kokkos::subview(b, 0, Kokkos::ALL());
                  auto r0 = Kokkos::subview(r, 0, Kokkos::ALL());

                  TeamCopy<typename policy_type::member_type, Trans::NoTranspose>::invoke(member, b0, r0);
                  TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A0,
                                                                                                     x0, 1.0, r0);
                } else {
                  int k = 0;
                  {
                    /// first row
                    auto A1 = Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL());
                    auto B2 = Kokkos::subview(B, k, Kokkos::ALL(), Kokkos::ALL());

                    auto x1 = Kokkos::subview(x, k, Kokkos::ALL());
                    auto x2 = Kokkos::subview(x, k + 1, Kokkos::ALL());

                    auto bk = Kokkos::subview(b, k, Kokkos::ALL());
                    auto rk = Kokkos::subview(r, k, Kokkos::ALL());
                    TeamCopy<typename policy_type::member_type, Trans::NoTranspose>::invoke(member, bk, rk);
                    member.team_barrier();
                    TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A1,
                                                                                                       x1, 1.0, rk);
                    TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, B2,
                                                                                                       x2, 1.0, rk);
                    ++k;
                  }
                  for (; k < (L - 1); ++k) {
                    auto C0 = Kokkos::subview(C, k - 1, Kokkos::ALL(), Kokkos::ALL());
                    auto A1 = Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL());
                    auto B2 = Kokkos::subview(B, k, Kokkos::ALL(), Kokkos::ALL());

                    auto x0 = Kokkos::subview(x, k - 1, Kokkos::ALL());
                    auto x1 = Kokkos::subview(x, k, Kokkos::ALL());
                    auto x2 = Kokkos::subview(x, k + 1, Kokkos::ALL());

                    auto bk = Kokkos::subview(b, k, Kokkos::ALL());
                    auto rk = Kokkos::subview(r, k, Kokkos::ALL());
                    TeamCopy<typename policy_type::member_type, Trans::NoTranspose>::invoke(member, bk, rk);
                    member.team_barrier();
                    TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, C0,
                                                                                                       x0, 1.0, rk);
                    TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A1,
                                                                                                       x1, 1.0, rk);
                    TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, B2,
                                                                                                       x2, 1.0, rk);
                  }
                  {
                    // last row
                    auto C0 = Kokkos::subview(C, k - 1, Kokkos::ALL(), Kokkos::ALL());
                    auto A1 = Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL());

                    auto x0 = Kokkos::subview(x, k - 1, Kokkos::ALL());
                    auto x1 = Kokkos::subview(x, k, Kokkos::ALL());

                    auto bk = Kokkos::subview(b, k, Kokkos::ALL());
                    auto rk = Kokkos::subview(r, k, Kokkos::ALL());
                    TeamCopy<typename policy_type::member_type, Trans::NoTranspose>::invoke(member, bk, rk);
                    member.team_barrier();
                    TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, C0,
                                                                                                       x0, 1.0, rk);
                    TeamGemv<typename policy_type::member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A1,
                                                                                                       x1, 1.0, rk);
                  }
                }
              }
            });
          });
      Kokkos::fence();
      auto rs_host = Kokkos::create_mirror_view(rs);
      auto bs_host = Kokkos::create_mirror_view(bs);
      Kokkos::deep_copy(rs_host, rs);
      Kokkos::deep_copy(bs_host, bs);
      Kokkos::fence();
      {
        double norm2 = 0, diff2 = 0;
        for (int i0 = 0, i0end = rs.extent(0); i0 < i0end; ++i0)            // N/vector_length
          for (int i1 = 0, i1end = rs.extent(1); i1 < i1end; ++i1)          // Nvec
            for (int i2 = 0, i2end = rs.extent(2); i2 < i2end; ++i2)        // L
              for (int i3 = 0, i3end = rs.extent(3); i3 < i3end; ++i3)      // Blk
                for (int i4 = 0, i4end = rs.extent(4); i4 < i4end; ++i4) {  // vector_length
                  const auto val = bs_host(i0, i1, i2, i3, i4);
                  const auto res = rs_host(i0, i1, i2, i3, i4);
                  norm2 += val * val;
                  diff2 += res * res;
                }
        printf("rel error = %e\n", diff2 / norm2);

        // const int i0 = 0;
        // const int i1 = 0;
        // const int i4 = 0;
        // for (int i2=0;i2<rs.extent(2);++i2) // L
        //   for (int i3=0;i3<rs.extent(3);++i3) // Blk
        //     printf("row %d, block row %d residual %e\n",
        //            i2, i3, rs_host(i0,i1,i2,i3,i4));
      }
    }
  }
  Kokkos::finalize();

  return 0;
}
#else
int main() { return 0; }
#endif
