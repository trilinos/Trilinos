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

/// KokkosKernels headers
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Copy_Impl.hpp>
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

using exec_space_type   = Kokkos::DefaultExecutionSpace;
using memory_space_type = exec_space_type::memory_space;
using host_space_type   = Kokkos::DefaultHostExecutionSpace;

using value_type  = double;
using policy_type = Kokkos::TeamPolicy<exec_space_type>;
using member_type = typename policy_type::member_type;

/// 128*128*128/16*5 * (2*8) / 16
///
/// simd typedefs
///
using namespace KokkosBatched;

static constexpr int vector_length = DefaultVectorLength<value_type, memory_space_type>::value;
#if defined(KOKKOSBATCHED_USE_128BIT_MEMORY_INST)
static constexpr int internal_vector_length = DefaultInternalVectorLength<value_type, memory_space_type>::value;
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
struct FactorizeModeAndAlgo;

struct FactorizeModeAndAlgoHostImpl {
  typedef Mode::Serial mode_type;
  typedef Algo::Level3::Blocked algo_type;
};

#if defined(KOKKOS_ENABLE_SERIAL)
template <>
struct FactorizeModeAndAlgo<Kokkos::Serial> : FactorizeModeAndAlgoHostImpl {};
#endif

#if defined(KOKKOS_ENABLE_THREADS)
template <>
struct FactorizeModeAndAlgo<Kokkos::Threads> : FactorizeModeAndAlgoHostImpl {};
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
template <>
struct FactorizeModeAndAlgo<Kokkos::OpenMP> : FactorizeModeAndAlgoHostImpl {};
#endif

struct FactorizeModeAndAlgoDeviceImpl {
  typedef Mode::Team mode_type;
  typedef Algo::Level3::Unblocked algo_type;
};

#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct FactorizeModeAndAlgo<Kokkos::Cuda> : FactorizeModeAndAlgoDeviceImpl {};
#endif

#if defined(KOKKOS_ENABLE_HIP)
template <>
struct FactorizeModeAndAlgo<Kokkos::HIP> : FactorizeModeAndAlgoDeviceImpl {};
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

template <class VT>
struct SetTridiagToIdentity {
 private:
  VT __AA;

 public:
  SetTridiagToIdentity(VT AA) : __AA(AA) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member) const {
    const int i = member.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, __AA.extent(1)), [&](const int &j) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, __AA.extent(5)), [&](const int &v) {
        for (int k = 0, kend = __AA.extent(3); k < kend; ++k) __AA(i, j, 1, k, k, v) = 1;
      });
    });
  }
};

template <class VT, class LT>
struct Factorize {
 private:
  VT __AA;
  LT __L;

 public:
  Factorize(VT AA, LT L) : __AA(AA), __L(L) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type &member) const {
    typedef FactorizeModeAndAlgo<Kokkos::DefaultExecutionSpace> default_mode_and_algo_type;
    typedef default_mode_and_algo_type::mode_type mode_type;
    typedef default_mode_and_algo_type::algo_type algo_type;

    const int i = member.league_rank();

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, __AA.extent(5)), [&](const int &v) {
      auto AAA = Kokkos::subview(__AA, i, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), v);

      /// subview patterns
      auto A = Kokkos::subview(AAA, 0, 1, Kokkos::ALL(), Kokkos::ALL());
      auto B = Kokkos::subview(AAA, 0, 2, Kokkos::ALL(), Kokkos::ALL());
      auto C = Kokkos::subview(AAA, 0, 0, Kokkos::ALL(), Kokkos::ALL());
      auto D = Kokkos::subview(AAA, 0, 1, Kokkos::ALL(), Kokkos::ALL());

      if (__L == 1) {
        A.assign_data(&AAA(0, 1, 0, 0));
        LU<member_type, mode_type, algo_type>::invoke(member, A);
      } else {
        for (int k = 0; k < (__L - 1); ++k) {
          A.assign_data(&AAA(k, 1, 0, 0));
          B.assign_data(&AAA(k, 2, 0, 0));
          C.assign_data(&AAA(k, 0, 0, 0));
          D.assign_data(&AAA(k + 1, 1, 0, 0));

          LU<member_type, mode_type, algo_type>::invoke(member, A);
          Trsm<member_type, Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, mode_type, algo_type>::invoke(
              member, 1.0, A, B);
          Trsm<member_type, Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, mode_type, algo_type>::invoke(
              member, 1.0, A, C);
          Gemm<member_type, Trans::NoTranspose, Trans::NoTranspose, mode_type, algo_type>::invoke(member, -1.0, C, B,
                                                                                                  1.0, D);
        }
        LU<member_type, mode_type, algo_type>::invoke(member, D);
      }
    });
  }
};

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
    int N     = 128 * 128;  /// # of problems (batch size)
    int L     = 128;        /// length of block tridiags
    int Blk   = 5;          /// block dimension
    int Nvec  = 1;
    int S     = 0;  /// scratch size
    int niter = 1;
    for (int i = 1; i < argc; ++i) {
      const std::string &token = argv[i];
      if (token == std::string("-N")) N = std::atoi(argv[++i]);
      if (token == std::string("-L")) L = std::atoi(argv[++i]);
      if (token == std::string("-B")) Blk = std::atoi(argv[++i]);
      if (token == std::string("-Nvec")) Nvec = std::atoi(argv[++i]);
      if (token == std::string("-S")) S = std::atoi(argv[++i]);
      if (token == std::string("-Niter")) niter = std::atoi(argv[++i]);
    }

    printf(
        " :::: Testing (N = %d, L = %d, Blk = %d, vl = %d, vi = %d, niter = "
        "%d)\n",
        N, L, Blk, vector_length, internal_vector_length, niter);

    ///
    /// problem container
    ///

    /// double 16
    Kokkos::View<vector_type *****, Kokkos::LayoutRight, exec_space_type> Av("A", N / vector_length, L, 3, Blk, Blk);

    /// double
    Kokkos::View<value_type ******, Kokkos::LayoutRight, exec_space_type> As(
        (value_type *)Av.data(), Av.extent(0), Av.extent(1), Av.extent(2), Av.extent(3), Av.extent(4), vector_length);

    /// double 2
    Kokkos::View<internal_vector_type ******, Kokkos::LayoutRight, exec_space_type> Ai(
        (internal_vector_type *)Av.data(), Av.extent(0), Av.extent(1), Av.extent(2), Av.extent(3), Av.extent(4),
        vector_length / internal_vector_length);
    /// double 16
    Kokkos::View<vector_type ****, Kokkos::LayoutRight, exec_space_type> xv("x", N / vector_length, Nvec, L, Blk);

    /// double
    Kokkos::View<value_type *****, Kokkos::LayoutRight, exec_space_type> xs(
        (value_type *)xv.data(), xv.extent(0), xv.extent(1), xv.extent(2), xv.extent(3), vector_length);

    /// double 2
    Kokkos::View<internal_vector_type *****, Kokkos::LayoutRight, exec_space_type> xi(
        (internal_vector_type *)xv.data(), xv.extent(0), xv.extent(1), xv.extent(2), xv.extent(3),
        vector_length / internal_vector_length);

    /// double 16
    Kokkos::View<vector_type ****, Kokkos::LayoutRight, exec_space_type> bv("b", N / vector_length, Nvec, L, Blk);

    /// double
    Kokkos::View<value_type *****, Kokkos::LayoutRight, exec_space_type> bs(
        (value_type *)bv.data(), bv.extent(0), bv.extent(1), bv.extent(2), bv.extent(3), vector_length);

    /// double 2
    Kokkos::View<internal_vector_type *****, Kokkos::LayoutRight, exec_space_type> bi(
        (internal_vector_type *)bv.data(), bv.extent(0), bv.extent(1), bv.extent(2), bv.extent(3),
        vector_length / internal_vector_length);

    /// double copy of A
    Kokkos::View<value_type ******, Kokkos::LayoutRight, exec_space_type> Acopy(
        "Acopy", As.extent(0), As.extent(1), As.extent(2), As.extent(3), As.extent(4), As.extent(5));

    Kokkos::View<value_type *****, Kokkos::LayoutRight, exec_space_type> rs("rs", bs.extent(0), bs.extent(1),
                                                                            bs.extent(2), bs.extent(3), bs.extent(4));

#if defined(KOKKOSBATCHED_USE_128BIT_MEMORY_INST)
    auto AA = Ai;
    auto bb = bi;
    auto xx = xi;
#else
    auto AA = As;
    auto bb = bs;
    auto xx = xs;
#endif

    ///
    /// set identity
    ///
    if (0) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStart();
#endif
      timer.reset();
      policy_type policy(AA.extent(0), Kokkos::AUTO(), AA.extent(5));
      Kokkos::parallel_for("setTridiagToIdentity", policy, SetTridiagToIdentity<decltype(AA)>(AA));
      Kokkos::fence();
      const double t = timer.seconds();
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStop();
#endif
      printf("identity time = %f\n", t);
    }

    /// randomize input
    {
      const value_type one(1);
      Kokkos::Random_XorShift64_Pool<exec_space_type> random(13245);
      Kokkos::fill_random(As, random, one);
      Kokkos::fill_random(bs, random, one);

      Kokkos::deep_copy(Acopy, As);
    }

    ///
    /// factorize the matrix
    ///
    if (1) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStart();
#endif
      timer.reset();
      int team_size = 0;
      if (Blk < 8) {
        team_size = 32 / AA.extent(5);
      } else if (Blk < 12) {
        team_size = 64 / AA.extent(5);
      } else {
        team_size = 128 / AA.extent(5);
      }

      policy_type policy(AA.extent(0), team_size, AA.extent(5));
      Kokkos::parallel_for("factorize", policy.set_scratch_size(0, Kokkos::PerTeam(S)),
                           Factorize<decltype(AA), decltype(L)>(AA, L));
      Kokkos::fence();
      const double t = timer.seconds();
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStop();
#endif
      printf("factorize time = %f , # of factorization per min = %f \n", t, 1.0 / t * 60);
    }

    ///
    /// solve the matrix 20 times
    ///
    if (1) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSBATCHED_PROFILE)
      cudaProfilerStart();
#endif
      timer.reset();
      int team_size = 0;
      if (Blk < 8) {
        team_size = 32 / AA.extent(5);
      } else if (Blk < 12) {
        team_size = 64 / AA.extent(5);
      } else {
        team_size = 128 / AA.extent(5);
      }

      policy_type policy(AA.extent(0), team_size, AA.extent(5));
      for (int iter = 0; iter < niter; ++iter) {
        Kokkos::parallel_for(
            "solve", policy.set_scratch_size(0, Kokkos::PerTeam(S)), KOKKOS_LAMBDA(const member_type &member) {
              typedef SolveModeAndAlgo<Kokkos::DefaultExecutionSpace> default_mode_and_algo_type;
              typedef default_mode_and_algo_type::mode_type mode_type;
              typedef default_mode_and_algo_type::algo_type algo_type;

              const int i = member.league_rank();
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, AA.extent(5)), [&](const int &v) {
                auto A = Kokkos::subview(AA, i, Kokkos::ALL(), 1, Kokkos::ALL(), Kokkos::ALL(), v);
                auto B = Kokkos::subview(AA, i, Kokkos::ALL(), 2, Kokkos::ALL(), Kokkos::ALL(), v);
                auto C = Kokkos::subview(AA, i, Kokkos::ALL(), 0, Kokkos::ALL(), Kokkos::ALL(), v);

                for (int jvec = 0; jvec < Nvec; ++jvec) {
                  auto x = Kokkos::subview(xx, i, jvec, Kokkos::ALL(), Kokkos::ALL(), v);
                  auto b = Kokkos::subview(bb, i, jvec, Kokkos::ALL(), Kokkos::ALL(), v);

                  auto xt = Kokkos::subview(x, 0, Kokkos::ALL());
                  auto xb = Kokkos::subview(x, 0, Kokkos::ALL());

                  ///
                  /// forward substitution
                  ///
                  {
                    // const bool is_same_x_and_b = (x.data() == b.data());
                    auto LT = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());
                    auto LB = Kokkos::subview(C, 0, Kokkos::ALL(), Kokkos::ALL());

                    auto bk = Kokkos::subview(b, 0, Kokkos::ALL());
                    {
                      {  // if (!is_same_x_and_b) {
                        Copy<member_type, Trans::NoTranspose, mode_type>::invoke(member, bk, xb);
                        member.team_barrier();
                      }
                    }
                    const int kend = L - 1;
                    for (int k = 0; k < kend; ++k) {
                      LT.assign_data(&A(k, 0, 0));
                      LB.assign_data(&C(k, 0, 0));

                      xt.assign_data(&x(k, 0));
                      xb.assign_data(&x(k + 1, 0));

                      {  // if (!is_same_x_and_b) {
                        bk.assign_data(&b(k + 1, 0));
                        Copy<member_type, Trans::NoTranspose, mode_type>::invoke(member, bk, xb);
                      }

                      Trsv<member_type, Uplo::Lower, Trans::NoTranspose, Diag::Unit, mode_type, algo_type>::invoke(
                          member, 1.0, LT, xt);

                      Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, -1.0, LB, xt, 1.0,
                                                                                          xb);
                    }
                    {
                      LT.assign_data(&A(kend, 0, 0));
                      xt.assign_data(&x(kend, 0));
                      Trsv<member_type, Uplo::Lower, Trans::NoTranspose, Diag::Unit, mode_type, algo_type>::invoke(
                          member, 1.0, LT, xt);
                    }
                  }  /// end forward substitution

                  ///
                  /// backward substitution
                  ///
                  {
                    auto UT = Kokkos::subview(B, 0, Kokkos::ALL(), Kokkos::ALL());
                    auto UB = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());

                    const int kbegin = L - 1;
                    for (int k = kbegin; k > 0; --k) {
                      UT.assign_data(&B(k - 1, 0, 0));
                      UB.assign_data(&A(k, 0, 0));

                      xt.assign_data(&x(k - 1, 0));
                      xb.assign_data(&x(k, 0));

                      Trsv<member_type, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, mode_type, algo_type>::invoke(
                          member, 1.0, UB, xb);

                      Gemv<member_type, Trans::NoTranspose, mode_type, algo_type>::invoke(member, -1.0, UT, xb, 1.0,
                                                                                          xt);
                    }
                    {
                      UB.assign_data(&A(0, 0, 0));
                      xb.assign_data(&x(0, 0));
                      Trsv<member_type, Uplo::Upper, Trans::NoTranspose, Diag::NonUnit, mode_type, algo_type>::invoke(
                          member, 1.0, UB, xb);
                    }
                  }  // end backward substitution
                }
              });
            });
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
      policy_type policy(Acopy.extent(0), Kokkos::AUTO(), Acopy.extent(5));
      Kokkos::parallel_for(
          "compute residual", policy, KOKKOS_LAMBDA(const member_type &member) {
            const int i = member.league_rank();
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, Acopy.extent(5)), [&](const int &v) {
              auto A = Kokkos::subview(Acopy, i, Kokkos::ALL(), 1, Kokkos::ALL(), Kokkos::ALL(), v);
              auto B = Kokkos::subview(Acopy, i, Kokkos::ALL(), 2, Kokkos::ALL(), Kokkos::ALL(), v);
              auto C = Kokkos::subview(Acopy, i, Kokkos::ALL(), 0, Kokkos::ALL(), Kokkos::ALL(), v);

              for (int jvec = 0, jvecend = rs.extent(1); jvec < jvecend; ++jvec) {
                auto x = Kokkos::subview(xs, i, jvec, Kokkos::ALL(), Kokkos::ALL(), v);
                auto b = Kokkos::subview(bs, i, jvec, Kokkos::ALL(), Kokkos::ALL(), v);
                auto r = Kokkos::subview(rs, i, jvec, Kokkos::ALL(), Kokkos::ALL(), v);

                if (L == 1) {
                  auto A0 = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());
                  auto x0 = Kokkos::subview(x, 0, Kokkos::ALL());
                  auto b0 = Kokkos::subview(b, 0, Kokkos::ALL());
                  auto r0 = Kokkos::subview(r, 0, Kokkos::ALL());

                  TeamCopy<member_type, Trans::NoTranspose>::invoke(member, b0, r0);
                  TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A0, x0, 1.0, r0);
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
                    TeamCopy<member_type, Trans::NoTranspose>::invoke(member, bk, rk);
                    member.team_barrier();
                    TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A1, x1, 1.0, rk);
                    TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, B2, x2, 1.0, rk);
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
                    TeamCopy<member_type, Trans::NoTranspose>::invoke(member, bk, rk);
                    member.team_barrier();
                    TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, C0, x0, 1.0, rk);
                    TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A1, x1, 1.0, rk);
                    TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, B2, x2, 1.0, rk);
                  }
                  {
                    // last row
                    auto C0 = Kokkos::subview(C, k - 1, Kokkos::ALL(), Kokkos::ALL());
                    auto A1 = Kokkos::subview(A, k, Kokkos::ALL(), Kokkos::ALL());

                    auto x0 = Kokkos::subview(x, k - 1, Kokkos::ALL());
                    auto x1 = Kokkos::subview(x, k, Kokkos::ALL());

                    auto bk = Kokkos::subview(b, k, Kokkos::ALL());
                    auto rk = Kokkos::subview(r, k, Kokkos::ALL());
                    TeamCopy<member_type, Trans::NoTranspose>::invoke(member, bk, rk);
                    member.team_barrier();
                    TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, C0, x0, 1.0, rk);
                    TeamGemv<member_type, Trans::NoTranspose, algo_type>::invoke(member, -1.0, A1, x1, 1.0, rk);
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
