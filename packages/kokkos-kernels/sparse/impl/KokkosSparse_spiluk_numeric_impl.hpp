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

#ifndef KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_HPP_
#define KOKKOSSPARSE_IMPL_SPILUK_NUMERIC_HPP_

/// \file KokkosSparse_spiluk_numeric_impl.hpp
/// \brief Implementation(s) of the numeric phase of sparse ILU(k).

#include <KokkosKernels_config.h>
#include <KokkosKernels_Error.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_spiluk_handle.hpp>
#include "KokkosBatched_SetIdentity_Decl.hpp"
#include "KokkosBatched_SetIdentity_Impl.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"
#include "KokkosBatched_Axpy.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBlas1_set.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_Trmm_Decl.hpp"
#include "KokkosBatched_Trmm_Serial_Impl.hpp"

// #define NUMERIC_OUTPUT_INFO

namespace KokkosSparse {
namespace Impl {
namespace Experimental {

template <class IlukHandle>
struct IlukWrap {
  //
  // Useful types
  //
  using execution_space   = typename IlukHandle::execution_space;
  using memory_space      = typename IlukHandle::memory_space;
  using lno_t             = typename IlukHandle::nnz_lno_t;
  using size_type         = typename IlukHandle::size_type;
  using scalar_t          = typename IlukHandle::nnz_scalar_t;
  using WorkViewType      = typename IlukHandle::work_view_t;
  using LevelHostViewType = typename IlukHandle::nnz_lno_view_host_t;
  using LevelViewType     = typename IlukHandle::nnz_lno_view_t;
  using karith            = typename Kokkos::ArithTraits<scalar_t>;
  using team_policy       = typename IlukHandle::TeamPolicy;
  using member_type       = typename team_policy::member_type;
  using range_policy      = typename IlukHandle::RangePolicy;

  static team_policy get_team_policy(const size_type nrows, const int team_size) {
    team_policy rv;
    if (team_size == -1) {
      rv = team_policy(nrows, Kokkos::AUTO);
    } else {
      rv = team_policy(nrows, team_size);
    }

    return rv;
  }

  static team_policy get_team_policy(execution_space exe_space, const size_type nrows, const int team_size) {
    team_policy rv;
    if (team_size == -1) {
      rv = team_policy(exe_space, nrows, Kokkos::AUTO);
    } else {
      rv = team_policy(exe_space, nrows, team_size);
    }

    return rv;
  }

  /**
   * Common base class for SPILUK functors. Default version does not support
   * blocks
   */
  template <class ARowMapType, class AEntriesType, class AValuesType, class LRowMapType, class LEntriesType,
            class LValuesType, class URowMapType, class UEntriesType, class UValuesType, bool BlockEnabled>
  struct Common {
    ARowMapType A_row_map;
    AEntriesType A_entries;
    AValuesType A_values;
    LRowMapType L_row_map;
    LEntriesType L_entries;
    LValuesType L_values;
    URowMapType U_row_map;
    UEntriesType U_entries;
    UValuesType U_values;
    LevelViewType level_idx;
    WorkViewType iw;
    lno_t lev_start;

    using reftype = scalar_t &;
    using valtype = scalar_t;

    static constexpr size_type BUFF_SIZE = 1;

    struct SBlock {
      template <typename T>
      KOKKOS_INLINE_FUNCTION SBlock(T, size_type, size_type) {}

      KOKKOS_INLINE_FUNCTION
      scalar_t *data() { return nullptr; }
    };

    Common(const ARowMapType &A_row_map_, const AEntriesType &A_entries_, const AValuesType &A_values_,
           const LRowMapType &L_row_map_, const LEntriesType &L_entries_, LValuesType &L_values_,
           const URowMapType &U_row_map_, const UEntriesType &U_entries_, UValuesType &U_values_,
           const LevelViewType &level_idx_, WorkViewType &iw_, const lno_t &lev_start_, const size_type &block_size_)
        : A_row_map(A_row_map_),
          A_entries(A_entries_),
          A_values(A_values_),
          L_row_map(L_row_map_),
          L_entries(L_entries_),
          L_values(L_values_),
          U_row_map(U_row_map_),
          U_entries(U_entries_),
          U_values(U_values_),
          level_idx(level_idx_),
          iw(iw_),
          lev_start(lev_start_) {
      KK_REQUIRE_MSG(block_size_ == 0, "Tried to use blocks with the unblocked Common?");
    }

    KOKKOS_INLINE_FUNCTION
    size_type get_block_size() const { return 0; }

    // lset
    KOKKOS_INLINE_FUNCTION
    void lset(const size_type nnz, const scalar_t &value) const { L_values(nnz) = value; }

    // uset
    KOKKOS_INLINE_FUNCTION
    void uset(const size_type nnz, const scalar_t &value) const { U_values(nnz) = value; }

    // lset_id
    KOKKOS_INLINE_FUNCTION
    void lset_id(const member_type &team, const size_type nnz) const {
      // Not sure a Kokkos::single is really needed here since the
      // race is harmless
      Kokkos::single(Kokkos::PerTeam(team), [&]() { L_values(nnz) = scalar_t(1.0); });
    }

    // divide. lhs /= rhs
    KOKKOS_INLINE_FUNCTION
    void divide(const member_type &team, scalar_t &lhs, const scalar_t &rhs, scalar_t *) const {
      Kokkos::single(Kokkos::PerTeam(team), [&]() { lhs /= rhs; });
      team.team_barrier();
    }

    // divide_left. lhs /= rhs
    KOKKOS_INLINE_FUNCTION
    void divide_left(scalar_t &lhs, const scalar_t &rhs, scalar_t *) const { lhs /= rhs; }

    // multiply_subtract. C -= A * B
    KOKKOS_INLINE_FUNCTION
    void multiply_subtract(const scalar_t &A, const scalar_t &B, scalar_t &C) const { C -= A * B; }

    // lget
    KOKKOS_INLINE_FUNCTION
    scalar_t &lget(const size_type nnz) const { return L_values(nnz); }

    // lcopy
    KOKKOS_INLINE_FUNCTION
    scalar_t lcopy(const size_type nnz, scalar_t *) const { return L_values(nnz); }

    // ucopy
    KOKKOS_INLINE_FUNCTION
    scalar_t ucopy(const size_type nnz, scalar_t *) const { return U_values(nnz); }

    // uget
    KOKKOS_INLINE_FUNCTION
    scalar_t &uget(const size_type nnz) const { return U_values(nnz); }

    // aget
    KOKKOS_INLINE_FUNCTION
    scalar_t aget(const size_type nnz) const { return A_values(nnz); }

    // uequal
    KOKKOS_INLINE_FUNCTION
    bool uequal(const size_type nnz, const scalar_t &value) const { return U_values(nnz) == value; }

    // print
    KOKKOS_INLINE_FUNCTION
    void print(const scalar_t &item) const { std::cout << item << std::endl; }

    // report
    KOKKOS_INLINE_FUNCTION
    void report() const { std::cout << "JGF using unblocked version" << std::endl; }
  };

  // Partial specialization for block support
  template <class ARowMapType, class AEntriesType, class AValuesType, class LRowMapType, class LEntriesType,
            class LValuesType, class URowMapType, class UEntriesType, class UValuesType>
  struct Common<ARowMapType, AEntriesType, AValuesType, LRowMapType, LEntriesType, LValuesType, URowMapType,
                UEntriesType, UValuesType, true> {
    // BSR data is in LayoutRight!
    using Layout      = Kokkos::LayoutRight;
    using value_type  = typename LValuesType::value_type;
    using cvalue_type = typename LValuesType::const_value_type;

    using Block = Kokkos::View<value_type **, Layout, typename LValuesType::device_type,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

    // const block
    using CBlock = Kokkos::View<cvalue_type **, Layout, typename UValuesType::device_type,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

    // scratch block
    using SBlock = Kokkos::View<value_type **, Layout, typename execution_space::scratch_memory_space,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess> >;

    using reftype = Block;
    using valtype = Block;

    static constexpr size_type BUFF_SIZE = 128;

    ARowMapType A_row_map;
    AEntriesType A_entries;
    AValuesType A_values;
    LRowMapType L_row_map;
    LEntriesType L_entries;
    LValuesType L_values;
    URowMapType U_row_map;
    UEntriesType U_entries;
    UValuesType U_values;
    LevelViewType level_idx;
    WorkViewType iw;
    lno_t lev_start;
    size_type block_size;
    size_type block_items;

    Common(const ARowMapType &A_row_map_, const AEntriesType &A_entries_, const AValuesType &A_values_,
           const LRowMapType &L_row_map_, const LEntriesType &L_entries_, LValuesType &L_values_,
           const URowMapType &U_row_map_, const UEntriesType &U_entries_, UValuesType &U_values_,
           const LevelViewType &level_idx_, WorkViewType &iw_, const lno_t &lev_start_, const size_type &block_size_)
        : A_row_map(A_row_map_),
          A_entries(A_entries_),
          A_values(A_values_),
          L_row_map(L_row_map_),
          L_entries(L_entries_),
          L_values(L_values_),
          U_row_map(U_row_map_),
          U_entries(U_entries_),
          U_values(U_values_),
          level_idx(level_idx_),
          iw(iw_),
          lev_start(lev_start_),
          block_size(block_size_),
          block_items(block_size * block_size) {
      KK_REQUIRE_MSG(block_size > 0, "Tried to use block_size=0 with the blocked Common?");
      KK_REQUIRE_MSG(block_size <= 11, "Max supported block size is 11");
    }

    KOKKOS_INLINE_FUNCTION
    size_type get_block_size() const { return block_size; }

    // lset
    KOKKOS_INLINE_FUNCTION
    void lset(const size_type block, const scalar_t &value) const { KokkosBlas::SerialSet::invoke(value, lget(block)); }

    KOKKOS_INLINE_FUNCTION
    void lset(const size_type block, const CBlock &rhs) const {
      auto lblock = lget(block);
      assign(lblock, rhs);
    }

    // uset
    KOKKOS_INLINE_FUNCTION
    void uset(const size_type block, const scalar_t &value) const { KokkosBlas::SerialSet::invoke(value, uget(block)); }

    KOKKOS_INLINE_FUNCTION
    void uset(const size_type block, const CBlock &rhs) const {
      auto ublock = uget(block);
      assign(ublock, rhs);
    }

    // lset_id
    KOKKOS_INLINE_FUNCTION
    void lset_id(const member_type &team, const size_type block) const {
      KokkosBatched::TeamSetIdentity<member_type>::invoke(team, lget(block));
    }

    // assign
    template <typename ViewT>
    KOKKOS_INLINE_FUNCTION void assign(const ViewT &lhs, const CBlock &rhs) const {
      for (size_type i = 0; i < block_size; ++i) {
        for (size_type j = 0; j < block_size; ++j) {
          lhs(i, j) = rhs(i, j);
        }
      }
    }

    // divide. lhs /= rhs (lhs = lhs * rhs^-1)
    KOKKOS_INLINE_FUNCTION
    void divide(const member_type &team, const Block &lhs, const CBlock &rhs, scalar_t *buff) const {
      // Need a temp block to do LU of rhs
      Block LU(buff, block_size, block_size);
      assign(LU, rhs);
      KokkosBatched::TeamLU<member_type, KokkosBatched::Algo::LU::Blocked>::invoke(team, LU);

      // rhs = LU
      // rhs^-1 = U^-1 * L^-1
      // lhs = (lhs * U^-1) * L^-1, so do U trsm first
      KokkosBatched::TeamTrsm<member_type, KokkosBatched::Side::Right, KokkosBatched::Uplo::Upper,
                              KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                              KokkosBatched::Algo::Trsm::Blocked>::invoke(team, 1.0, LU, lhs);

      KokkosBatched::TeamTrsm<member_type, KokkosBatched::Side::Right, KokkosBatched::Uplo::Lower,
                              KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit,
                              KokkosBatched::Algo::Trsm::Blocked>::invoke(team, 1.0, LU, lhs);
    }

    // divide_left. lhs /= rhs (lhs = rhs^-1 * lhs)
    KOKKOS_INLINE_FUNCTION
    void divide_left(const Block &lhs, const CBlock &rhs, scalar_t *buff) const {
      Block LU(buff, block_size, block_size);
      assign(LU, rhs);
      KokkosBatched::SerialLU<KokkosBatched::Algo::LU::Blocked>::invoke(LU);

      // rhs = LU
      // rhs^-1 = U^-1 * L^-1
      // lhs = U^-1 * (L^-1 * lhs), so do L trsm first
      KokkosBatched::SerialTrsm<KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit,
                                KokkosBatched::Algo::Trsm::Blocked>::invoke(1.0, LU, lhs);

      KokkosBatched::SerialTrsm<KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                KokkosBatched::Algo::Trsm::Blocked>::invoke(1.0, LU, lhs);
    }

    // multiply_subtract. C -= A * B
    KOKKOS_INLINE_FUNCTION
    void multiply_subtract(const CBlock &A, const CBlock &B, const Block &C) const {
      // Use gemm. alpha is hardcoded to -1, beta hardcoded to 1
      KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose,
                                KokkosBatched::Algo::Gemm::Blocked>::invoke<scalar_t, CBlock, CBlock, Block>(-1.0, A, B,
                                                                                                             1.0, C);
    }

    // lget
    KOKKOS_INLINE_FUNCTION
    Block lget(const size_type block) const {
      return Block(L_values.data() + (block * block_items), block_size, block_size);
    }

    // lcopy
    KOKKOS_INLINE_FUNCTION
    Block lcopy(const size_type block, scalar_t *buff) const {
      Block result(buff, block_size, block_size);
      auto lblock = lget(block);
      assign(result, lblock);
      return result;
    }

    // ucopy
    KOKKOS_INLINE_FUNCTION
    Block ucopy(const size_type block, scalar_t *buff) const {
      Block result(buff, block_size, block_size);
      auto ublock = uget(block);
      assign(result, ublock);
      return result;
    }

    // uget
    KOKKOS_INLINE_FUNCTION
    Block uget(const size_type block) const {
      return Block(U_values.data() + (block * block_items), block_size, block_size);
    }

    // aget
    KOKKOS_INLINE_FUNCTION
    CBlock aget(const size_type block) const {
      return CBlock(A_values.data() + (block * block_items), block_size, block_size);
    }

    // uequal
    KOKKOS_INLINE_FUNCTION
    bool uequal(const size_type block, const scalar_t &value) const {
      auto u_block = uget(block);
      for (size_type i = 0; i < block_size; ++i) {
        for (size_type j = 0; j < block_size; ++j) {
          if (u_block(i, j) != value) {
            return false;
          }
        }
      }
      return true;
    }

    // print
    KOKKOS_INLINE_FUNCTION
    void print(const CBlock &item) const {
      for (size_type i = 0; i < block_size; ++i) {
        std::cout << "      ";
        for (size_type j = 0; j < block_size; ++j) {
          std::cout << item(i, j) << " ";
        }
        std::cout << std::endl;
      }
    }

    // report
    KOKKOS_INLINE_FUNCTION
    void report() const { std::cout << "JGF using blocked version with block_size=" << block_size << std::endl; }
  };

  template <class ARowMapType, class AEntriesType, class AValuesType, class LRowMapType, class LEntriesType,
            class LValuesType, class URowMapType, class UEntriesType, class UValuesType, bool BlockEnabled>
  struct ILUKLvlSchedTP1NumericFunctor
      : public Common<ARowMapType, AEntriesType, AValuesType, LRowMapType, LEntriesType, LValuesType, URowMapType,
                      UEntriesType, UValuesType, BlockEnabled> {
    using Base = Common<ARowMapType, AEntriesType, AValuesType, LRowMapType, LEntriesType, LValuesType, URowMapType,
                        UEntriesType, UValuesType, BlockEnabled>;

    ILUKLvlSchedTP1NumericFunctor(const ARowMapType &A_row_map_, const AEntriesType &A_entries_,
                                  const AValuesType &A_values_, const LRowMapType &L_row_map_,
                                  const LEntriesType &L_entries_, LValuesType &L_values_, const URowMapType &U_row_map_,
                                  const UEntriesType &U_entries_, UValuesType &U_values_,
                                  const LevelViewType &level_idx_, WorkViewType &iw_, const lno_t &lev_start_,
                                  const size_type &block_size_ = 0)
        : Base(A_row_map_, A_entries_, A_values_, L_row_map_, L_entries_, L_values_, U_row_map_, U_entries_, U_values_,
               level_idx_, iw_, lev_start_, block_size_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type &team) const {
      // Thread-local buffers. Use for Serial (non-team) work
      scalar_t buff1[Base::BUFF_SIZE];
      scalar_t buff2[Base::BUFF_SIZE];
      scalar_t buff3[Base::BUFF_SIZE];

      // Team-shared buffer. Use for team work.
      const auto bs = Base::get_block_size();
      typename Base::SBlock shared_buff(team.team_shmem(), bs, bs);

      const auto my_team = team.league_rank();
      const auto rowid   = Base::level_idx(my_team + Base::lev_start);  // map to rowid

      // Set active entries in L to zero, store active cols in iw
      // Set L diagonal for this row to identity
      size_type k1 = Base::L_row_map(rowid);
      size_type k2 = Base::L_row_map(rowid + 1) - 1;
      Base::lset_id(team, k2);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2), [&](const size_type k) {
        const auto col = Base::L_entries(k);
        Base::lset(k, 0.0);
        Base::iw(my_team, col) = k;
      });

      team.team_barrier();

      // Set active entries in U to zero, store active cols in iw
      k1 = Base::U_row_map(rowid);
      k2 = Base::U_row_map(rowid + 1);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2), [&](const size_type k) {
        const auto col = Base::U_entries(k);
        Base::uset(k, 0.0);
        Base::iw(my_team, col) = k;
      });

      team.team_barrier();

      // Unpack the rowid-th row of A, copy into L,U
      k1 = Base::A_row_map(rowid);
      k2 = Base::A_row_map(rowid + 1);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2), [&](const size_type k) {
        const auto col  = Base::A_entries(k);
        const auto ipos = Base::iw(my_team, col);
        if (col < rowid) {
          Base::lset(ipos, Base::aget(k));
        } else {
          Base::uset(ipos, Base::aget(k));
        }
      });

      team.team_barrier();

      // Eliminate prev rows
      k1 = Base::L_row_map(rowid);
      k2 = Base::L_row_map(rowid + 1) - 1;
      for (auto k = k1; k < k2; k++) {
        const auto prev_row = Base::L_entries(k);
        const auto udiag    = Base::uget(Base::U_row_map(prev_row));
        typename Base::valtype fact;
        if (BlockEnabled) {
          fact = Base::lcopy(k, &buff1[0]);  // fact = copy(Lval(k))
          Base::divide(team, Base::lget(k), udiag,
                       shared_buff.data());  // Lval(k) *= udiag^-1
        } else {
          Base::divide(team, Base::lget(k), udiag, nullptr);
          fact = Base::lget(k);  // fact = Lval(k) / udiag
        }
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, Base::U_row_map(prev_row) + 1, Base::U_row_map(prev_row + 1)),
            [&](const size_type kk) {
              const auto col  = Base::U_entries(kk);
              const auto ipos = Base::iw(my_team, col);
              if (ipos != -1) {
                typename Base::reftype C = col < rowid ? Base::lget(ipos) : Base::uget(ipos);
                if (BlockEnabled) {
                  auto ucopy = Base::ucopy(kk, &buff2[0]);
                  Base::divide_left(ucopy, udiag,
                                    &buff3[0]);             // ucopy = udiag^-1 * Uval(kk)
                  Base::multiply_subtract(fact, ucopy, C);  // C -= Lval(k) * (udiag^-1 * Uval(kk))
                } else {
                  Base::multiply_subtract(fact, Base::uget(kk),
                                          C);  // C -= (Lval(k) / udiag) * Uval(kk)
                }
              }
            });  // end for kk

        team.team_barrier();
      }  // end for k

      // Seems hacky. Change zero diagonals to 1M?
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        const auto ipos = Base::iw(my_team, rowid);
        if (Base::uequal(ipos, 0.0)) {
          Base::uset(ipos, 1e6);
        }
      });

      team.team_barrier();
      // Reset
      k1 = Base::L_row_map(rowid);
      k2 = Base::L_row_map(rowid + 1) - 1;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2), [&](const size_type k) {
        const auto col         = Base::L_entries(k);
        Base::iw(my_team, col) = -1;
      });

      k1 = Base::U_row_map(rowid);
      k2 = Base::U_row_map(rowid + 1);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, k1, k2), [&](const size_type k) {
        const auto col         = Base::U_entries(k);
        Base::iw(my_team, col) = -1;
      });
    }
  };

#define FunctorTypeMacro(Functor, BlockEnabled)                                                                      \
  Functor<ARowMapType, AEntriesType, AValuesType, LRowMapType, LEntriesType, LValuesType, URowMapType, UEntriesType, \
          UValuesType, BlockEnabled>

#define KernelLaunchMacro(arow, aent, aval, lrow, lent, lval, urow, uent, uval, polc, name, lidx, iwv, lstrt, ftf, \
                          ftb, be, bs)                                                                             \
  if (be) {                                                                                                        \
    ftb functor(arow, aent, aval, lrow, lent, lval, urow, uent, uval, lidx, iwv, lstrt, bs);                       \
    const int scratch_size = ftb::SBlock::shmem_size(bs, bs);                                                      \
    polc                   = polc.set_scratch_size(0, Kokkos::PerTeam(scratch_size));                              \
    Kokkos::parallel_for(name, polc, functor);                                                                     \
  } else {                                                                                                         \
    ftf functor(arow, aent, aval, lrow, lent, lval, urow, uent, uval, lidx, iwv, lstrt);                           \
    Kokkos::parallel_for(name, polc, functor);                                                                     \
  }

  template <class ARowMapType, class AEntriesType, class AValuesType, class LRowMapType, class LEntriesType,
            class LValuesType, class URowMapType, class UEntriesType, class UValuesType>
  static void iluk_numeric(IlukHandle &thandle, const ARowMapType &A_row_map, const AEntriesType &A_entries,
                           const AValuesType &A_values, const LRowMapType &L_row_map, const LEntriesType &L_entries,
                           LValuesType &L_values, const URowMapType &U_row_map, const UEntriesType &U_entries,
                           UValuesType &U_values) {
    using TPF = FunctorTypeMacro(ILUKLvlSchedTP1NumericFunctor, false);
    using TPB = FunctorTypeMacro(ILUKLvlSchedTP1NumericFunctor, true);

    size_type nlevels        = thandle.get_num_levels();
    int team_size            = thandle.get_team_size();
    const auto block_size    = thandle.get_block_size();
    const auto block_enabled = thandle.is_block_enabled();

    LevelHostViewType level_ptr_h = thandle.get_host_level_ptr();
    LevelViewType level_idx       = thandle.get_level_idx();

    LevelHostViewType level_nchunks_h, level_nrowsperchunk_h;
    WorkViewType iw;

    level_nchunks_h       = thandle.get_level_nchunks();
    level_nrowsperchunk_h = thandle.get_level_nrowsperchunk();
    iw                    = thandle.get_iw();

    // Main loop must be performed sequential. Question: Try out Cuda's graph
    // stuff to reduce kernel launch overhead
    for (size_type lvl = 0; lvl < nlevels; ++lvl) {
      lno_t lev_start = level_ptr_h(lvl);
      lno_t lev_end   = level_ptr_h(lvl + 1);

      if ((lev_end - lev_start) != 0) {
        lno_t lvl_rowid_start = 0;
        lno_t lvl_nrows_chunk;
        for (int chunkid = 0; chunkid < level_nchunks_h(lvl); chunkid++) {
          if ((lvl_rowid_start + level_nrowsperchunk_h(lvl)) > (lev_end - lev_start))
            lvl_nrows_chunk = (lev_end - lev_start) - lvl_rowid_start;
          else
            lvl_nrows_chunk = level_nrowsperchunk_h(lvl);

          team_policy tpolicy = get_team_policy(lvl_nrows_chunk, team_size);
          KernelLaunchMacro(A_row_map, A_entries, A_values, L_row_map, L_entries, L_values, U_row_map, U_entries,
                            U_values, tpolicy, "parfor_tp1", level_idx, iw, lev_start + lvl_rowid_start, TPF, TPB,
                            block_enabled, block_size);
          Kokkos::fence();
          lvl_rowid_start += lvl_nrows_chunk;
        }
      }  // end if
    }    // end for lvl

// Output check
#ifdef NUMERIC_OUTPUT_INFO
    std::cout << "  iluk_numeric result: " << std::endl;

    std::cout << "  nnzL: " << thandle.get_nnzL() << std::endl;
    std::cout << "  L_row_map = ";
    for (size_type i = 0; i < thandle.get_nrows() + 1; ++i) {
      std::cout << L_row_map(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  L_entries = ";
    for (size_type i = 0; i < thandle.get_nnzL(); ++i) {
      std::cout << L_entries(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  L_values = ";
    for (size_type i = 0; i < thandle.get_nnzL(); ++i) {
      std::cout << L_values(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  nnzU: " << thandle.get_nnzU() << std::endl;
    std::cout << "  U_row_map = ";
    for (size_type i = 0; i < thandle.get_nrows() + 1; ++i) {
      std::cout << U_row_map(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  U_entries = ";
    for (size_type i = 0; i < thandle.get_nnzU(); ++i) {
      std::cout << U_entries(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "  U_values = ";
    for (size_type i = 0; i < thandle.get_nnzU(); ++i) {
      std::cout << U_values(i) << " ";
    }
    std::cout << std::endl;
#endif

  }  // end iluk_numeric

  template <class ExecutionSpace, class ARowMapType, class AEntriesType, class AValuesType, class LRowMapType,
            class LEntriesType, class LValuesType, class URowMapType, class UEntriesType, class UValuesType>
  static void iluk_numeric_streams(const std::vector<ExecutionSpace> &execspace_v,
                                   const std::vector<IlukHandle *> &thandle_v,
                                   const std::vector<ARowMapType> &A_row_map_v,
                                   const std::vector<AEntriesType> &A_entries_v,
                                   const std::vector<AValuesType> &A_values_v,
                                   const std::vector<LRowMapType> &L_row_map_v,
                                   const std::vector<LEntriesType> &L_entries_v, std::vector<LValuesType> &L_values_v,
                                   const std::vector<URowMapType> &U_row_map_v,
                                   const std::vector<UEntriesType> &U_entries_v, std::vector<UValuesType> &U_values_v) {
    using TPF = FunctorTypeMacro(ILUKLvlSchedTP1NumericFunctor, false);
    using TPB = FunctorTypeMacro(ILUKLvlSchedTP1NumericFunctor, true);

    // Create vectors for handles' data in streams
    int nstreams = execspace_v.size();
    std::vector<size_type> nlevels_v(nstreams);
    std::vector<LevelHostViewType> lvl_ptr_h_v(nstreams);
    std::vector<LevelViewType> lvl_idx_v(nstreams);  // device views
    std::vector<lno_t> lvl_start_v(nstreams);
    std::vector<lno_t> lvl_end_v(nstreams);
    std::vector<WorkViewType> iw_v(nstreams);  // device views
    std::vector<bool> stream_have_level_v(nstreams);
    std::vector<bool> is_block_enabled_v(nstreams);
    std::vector<size_type> block_size_v(nstreams);

    // Retrieve data from handles and find max. number of levels among streams
    size_type nlevels_max = 0;
    for (int i = 0; i < nstreams; i++) {
      nlevels_v[i]           = thandle_v[i]->get_num_levels();
      lvl_ptr_h_v[i]         = thandle_v[i]->get_host_level_ptr();
      lvl_idx_v[i]           = thandle_v[i]->get_level_idx();
      iw_v[i]                = thandle_v[i]->get_iw();
      is_block_enabled_v[i]  = thandle_v[i]->is_block_enabled();
      block_size_v[i]        = thandle_v[i]->get_block_size();
      stream_have_level_v[i] = true;
      if (nlevels_max < nlevels_v[i]) nlevels_max = nlevels_v[i];
    }

    std::vector<LevelHostViewType> lvl_nchunks_h_v(nstreams);
    std::vector<LevelHostViewType> lvl_nrowsperchunk_h_v(nstreams);
    std::vector<lno_t> lvl_rowid_start_v(nstreams);
    std::vector<int> team_size_v(nstreams);

    for (int i = 0; i < nstreams; i++) {
      lvl_nchunks_h_v[i]       = thandle_v[i]->get_level_nchunks();
      lvl_nrowsperchunk_h_v[i] = thandle_v[i]->get_level_nrowsperchunk();
      team_size_v[i]           = thandle_v[i]->get_team_size();
    }

    // Main loop must be performed sequential
    for (size_type lvl = 0; lvl < nlevels_max; lvl++) {
      // Initial work across streams at each level
      lno_t lvl_nchunks_max = 0;
      for (int i = 0; i < nstreams; i++) {
        // Only do this if this stream has this level
        if (lvl < nlevels_v[i]) {
          lvl_start_v[i] = lvl_ptr_h_v[i](lvl);
          lvl_end_v[i]   = lvl_ptr_h_v[i](lvl + 1);
          if ((lvl_end_v[i] - lvl_start_v[i]) != 0) {
            stream_have_level_v[i] = true;
            lvl_rowid_start_v[i]   = 0;
            if (lvl_nchunks_max < lvl_nchunks_h_v[i](lvl)) lvl_nchunks_max = lvl_nchunks_h_v[i](lvl);
          } else
            stream_have_level_v[i] = false;
        } else
          stream_have_level_v[i] = false;
      }

      // Main work of the level across streams -- looping through chunnks
      for (int chunkid = 0; chunkid < lvl_nchunks_max; chunkid++) {
        // 1. Launch work on all streams (for each chunk)
        for (int i = 0; i < nstreams; i++) {
          // Launch only if stream i-th has this level
          if (stream_have_level_v[i]) {
            // Launch only if stream i-th has this chunk
            if (chunkid < lvl_nchunks_h_v[i](lvl)) {
              // 1.a. Specify number of rows (i.e. number of teams) to launch
              lno_t lvl_nrows_chunk = 0;
              if ((lvl_rowid_start_v[i] + lvl_nrowsperchunk_h_v[i](lvl)) > (lvl_end_v[i] - lvl_start_v[i]))
                lvl_nrows_chunk = (lvl_end_v[i] - lvl_start_v[i]) - lvl_rowid_start_v[i];
              else
                lvl_nrows_chunk = lvl_nrowsperchunk_h_v[i](lvl);

              // 1.b. Create functor for stream i-th and launch
              team_policy tpolicy = get_team_policy(execspace_v[i], lvl_nrows_chunk, team_size_v[i]);
              KernelLaunchMacro(A_row_map_v[i], A_entries_v[i], A_values_v[i], L_row_map_v[i], L_entries_v[i],
                                L_values_v[i], U_row_map_v[i], U_entries_v[i], U_values_v[i], tpolicy, "parfor_tp1",
                                lvl_idx_v[i], iw_v[i], lvl_start_v[i] + lvl_rowid_start_v[i], TPF, TPB,
                                is_block_enabled_v[i], block_size_v[i]);
              // 1.c. Ready to move to next chunk
              lvl_rowid_start_v[i] += lvl_nrows_chunk;
            }  // end if (chunkid < lvl_nchunks_h_v[i](lvl))
          }    // end if (stream_have_level_v[i])
        }      // end for streams
      }        // end for chunkid
    }          // end for lvl
  }            // end iluk_numeric_streams

};  // IlukWrap

}  // namespace Experimental
}  // namespace Impl
}  // namespace KokkosSparse

#undef FunctorTypeMacro
#undef KernelLaunchMacro

#endif
