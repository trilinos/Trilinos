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

#ifndef KOKKOSSPARSE_IMPL_SPTRSV_SOLVE_HPP_
#define KOKKOSSPARSE_IMPL_SPTRSV_SOLVE_HPP_

/// \file KokkosSparse_impl_sptrsv.hpp
/// \brief Implementation(s) of sparse triangular solve.

#include <KokkosKernels_config.h>
#include <Kokkos_ArithTraits.hpp>
#include <KokkosSparse_sptrsv_handle.hpp>
#include <KokkosSparse_spmv.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
// Enable supernodal sptrsv
#include "KokkosSparse_spmv.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBlas2_team_gemv_spec.hpp"
#endif
#include "KokkosBlas3_trsm.hpp"
#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_Trsm_Team_Impl.hpp"
#include "KokkosBlas1_team_axpby.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosBlas1_set.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBlas2_serial_gemv_impl.hpp"

#define KOKKOSKERNELS_SPTRSV_TRILVLSCHED

// #define KOKKOSPSTRSV_SOLVE_IMPL_PROFILE 1
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
#include "cuda_profiler_api.h"
#endif

#if defined(KOKKOS_ENABLE_CUDA) && 10000 < CUDA_VERSION && defined(KOKKOSKERNELS_ENABLE_EXP_CUDAGRAPH)
#define KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
#endif

namespace KokkosSparse {
namespace Impl {
namespace Experimental {

template <class TriSolveHandle>
struct SptrsvWrap {
  //
  // Useful types
  //
  using execution_space = typename TriSolveHandle::execution_space;
  using memory_space    = typename TriSolveHandle::memory_space;
  using temp_mem_space  = typename TriSolveHandle::HandleTempMemorySpace;
  using lno_t           = typename TriSolveHandle::nnz_lno_t;
  using size_type       = typename TriSolveHandle::size_type;
  using scalar_t        = typename TriSolveHandle::scalar_t;
  using row_map_t       = typename TriSolveHandle::nnz_row_view_t;
  using entries_t       = typename TriSolveHandle::nnz_lno_view_t;
  using values_t        = typename TriSolveHandle::nnz_scalar_view_t;
  using work_view_t     = Kokkos::View<scalar_t *, Kokkos::Device<execution_space, temp_mem_space>>;
  using work_view_int_t = Kokkos::View<int *, Kokkos::Device<execution_space, temp_mem_space>>;
  using karith          = typename Kokkos::ArithTraits<scalar_t>;
  using team_policy     = typename TriSolveHandle::TeamPolicy;
  using member_type     = typename team_policy::member_type;
  using range_policy    = typename TriSolveHandle::RangePolicy;
  using range_type      = Kokkos::pair<int, int>;

  // Tag structs
  struct UnsortedTag {};  //  This doesn't appear to be supported
  struct LargerCutoffTag {};
  struct UnsortedLargerCutoffTag {};

  template <class ViewType>
  static void print_view1d_solve(const ViewType dv, size_t range = 0) {
    auto v = Kokkos::create_mirror_view(dv);
    Kokkos::deep_copy(v, dv);
    std::cout << "Output for view " << v.label() << std::endl;
    range = range == 0 ? dv.extent(0) : range;
    for (size_t i = 0; i < range; ++i) {
      std::cout << "v(" << i << ") = " << v(i) << " , ";
    }
    std::cout << std::endl;
  }

  // Needed for cudagraphs
  struct EmptyFunctor {
    KOKKOS_INLINE_FUNCTION
    void operator()(const int) const {}
  };

  /**
   * Common base class for sptrsv functors that need to work for both
   * point and block matrices. Default version does not support
   * blocks
   */
  template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, bool BlockEnabled>
  struct Common {
    RowMapType row_map;
    EntriesType entries;
    ValuesType values;
    LHSType lhs;
    RHSType rhs;
    entries_t nodes_grouped_by_level;

    using reftype = scalar_t &;

    struct SBlock {
      template <typename T>
      KOKKOS_INLINE_FUNCTION SBlock(T, size_type, size_type) {}

      KOKKOS_INLINE_FUNCTION
      scalar_t *data() { return nullptr; }

      static int shmem_size(size_type, size_type) { return 0; }
    };

    Common(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_,
           const RHSType &rhs_, const entries_t &nodes_grouped_by_level_, const size_type block_size_ = 0)
        : row_map(row_map_),
          entries(entries_),
          values(values_),
          lhs(lhs_),
          rhs(rhs_),
          nodes_grouped_by_level(nodes_grouped_by_level_) {
      KK_REQUIRE_MSG(block_size_ == 0, "Tried to use blocks with the unblocked Common?");
    }

    KOKKOS_INLINE_FUNCTION
    size_type get_block_size() const { return 0; }

    // lget
    KOKKOS_INLINE_FUNCTION
    scalar_t &lget(const size_type row) const { return lhs(row); }

    // rget
    KOKKOS_INLINE_FUNCTION
    scalar_t rget(const size_type row) const { return rhs(row); }

    // vget
    KOKKOS_INLINE_FUNCTION
    scalar_t vget(const size_type nnz) const { return values(nnz); }

    // lhs = (lhs + rhs) / diag (team)
    KOKKOS_INLINE_FUNCTION
    static void add_and_divide(const member_type &team, scalar_t &lhs_val, const scalar_t &rhs_val,
                               const scalar_t &diag_val) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() { lhs_val = (lhs_val + rhs_val) / diag_val; });
    }

    // lhs = (lhs + rhs) / diag (serial)
    KOKKOS_INLINE_FUNCTION
    static void add_and_divide(scalar_t &lhs_val, const scalar_t &rhs_val, const scalar_t &diag_val) {
      lhs_val = (lhs_val + rhs_val) / diag_val;
    }
  };

  // Partial specialization for block support
  template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType>
  struct Common<RowMapType, EntriesType, ValuesType, LHSType, RHSType, true> {
    // BSR data is in LayoutRight!
    using Layout = Kokkos::LayoutRight;

    using Block = Kokkos::View<scalar_t **, Layout, typename ValuesType::device_type,
                               Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

    // const block
    using CBlock = Kokkos::View<const scalar_t **, Layout, typename ValuesType::device_type,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

    // scratch block
    using SBlock = Kokkos::View<scalar_t **, Layout, typename execution_space::scratch_memory_space,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

    using Vector = Kokkos::View<scalar_t *, Layout, typename ValuesType::device_type,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

    using CVector = Kokkos::View<const scalar_t *, Layout, typename ValuesType::device_type,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;

    static constexpr size_type MAX_VEC_SIZE = 16;
    static constexpr size_type BUFF_SIZE    = 256;

    using reftype = Vector;

    RowMapType row_map;
    EntriesType entries;
    ValuesType values;
    LHSType lhs;
    RHSType rhs;
    entries_t nodes_grouped_by_level;
    size_type block_size;
    size_type block_items;

    Common(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_,
           const RHSType &rhs_, const entries_t &nodes_grouped_by_level_, const size_type block_size_)
        : row_map(row_map_),
          entries(entries_),
          values(values_),
          lhs(lhs_),
          rhs(rhs_),
          nodes_grouped_by_level(nodes_grouped_by_level_),
          block_size(block_size_),
          block_items(block_size * block_size) {
      KK_REQUIRE_MSG(block_size > 0, "Tried to use block_size=0 with the blocked Common?");
    }

    KOKKOS_INLINE_FUNCTION
    size_type get_block_size() const { return block_size; }

    // assign
    template <typename View1, typename View2>
    KOKKOS_INLINE_FUNCTION static void assign(const View1 &lhs_, const View2 &rhs_) {
      for (size_t i = 0; i < lhs_.size(); ++i) {
        lhs_.data()[i] = rhs_.data()[i];
      }
    }

    template <typename View1, typename View2>
    KOKKOS_INLINE_FUNCTION static void assign(const member_type &team, const View1 &lhs_, const View2 &rhs_) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, lhs_.size()),
                           [&](const size_type i) { lhs_.data()[i] = rhs_.data()[i]; });
    }

    // add. y += x
    KOKKOS_INLINE_FUNCTION
    static void add(const member_type &team, const CVector &x, const Vector &y) {
      KokkosBlas::Experimental::axpy(team, 1.0, x, y);
    }

    // serial add. y += x
    KOKKOS_INLINE_FUNCTION
    static void add(const CVector &x, const Vector &y) { KokkosBlas::serial_axpy(1.0, x, y); }

    // divide. b /= A (b = b * A^-1)
    KOKKOS_INLINE_FUNCTION
    static void divide(const member_type &team, const Vector &b, const CBlock &A) {
      // Team-shared buffer. Use for team work.
      const auto block_size_ = b.size();
      SBlock shared_buff(team.team_shmem(), block_size_, block_size_);

      // Need a temp block to do LU of A
      Block LU(shared_buff.data(), block_size_, block_size_);
      assign(team, LU, A);
      team.team_barrier();
      KokkosBatched::TeamLU<member_type, KokkosBatched::Algo::LU::Blocked>::invoke(team, LU);

      // A = LU
      // A^-1 = U^-1 * L^-1
      // b = (b * U^-1) * L^-1, so do U trsv first
      team.team_barrier();
      KokkosBatched::TeamTrsv<member_type, KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                              KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, LU,
                                                                                                        b);

      team.team_barrier();
      KokkosBatched::TeamTrsv<member_type, KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose,
                              KokkosBatched::Diag::Unit, KokkosBatched::Algo::Trsv::Blocked>::invoke(team, 1.0, LU, b);
    }

    // serial divide. b /= A (b = b * A^-1)
    KOKKOS_INLINE_FUNCTION
    static void divide(const Vector &b, const CBlock &A) {
      // Thread-local buffers. Use for Serial (non-team) work
      scalar_t buff[BUFF_SIZE];

      // Need a temp block to do LU of A
      const auto block_size_ = b.size();
      KK_KERNEL_REQUIRE_MSG(block_size_ <= MAX_VEC_SIZE,
                            "Max supported block size for range-policy is 16. Use team-policy alg if you need more.");

      Block LU(&buff[0], block_size_, block_size_);
      assign(LU, A);
      KokkosBatched::SerialLU<KokkosBatched::Algo::LU::Blocked>::invoke(LU);

      // A = LU
      // A^-1 = U^-1 * L^-1
      // b = (b * U^-1) * L^-1, so do U trsv first
      KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Upper, KokkosBatched::Trans::NoTranspose,
                                KokkosBatched::Diag::NonUnit, KokkosBatched::Algo::Trsv::Blocked>::invoke(1.0, LU, b);

      KokkosBatched::SerialTrsv<KokkosBatched::Uplo::Lower, KokkosBatched::Trans::NoTranspose,
                                KokkosBatched::Diag::Unit, KokkosBatched::Algo::Trsv::Blocked>::invoke(1.0, LU, b);
    }

    // lget
    KOKKOS_INLINE_FUNCTION
    Vector lget(const size_type row) const { return Vector(lhs.data() + (row * block_size), block_size); }

    // rget
    KOKKOS_INLINE_FUNCTION
    CVector rget(const size_type row) const { return CVector(rhs.data() + (row * block_size), block_size); }

    // vget
    KOKKOS_INLINE_FUNCTION
    CBlock vget(const size_type block) const {
      return CBlock(values.data() + (block * block_items), block_size, block_size);
    }

    // lhs = (lhs + rhs) / diag
    KOKKOS_INLINE_FUNCTION
    static void add_and_divide(const member_type &team, const Vector &lhs_val, const CVector &rhs_val,
                               const CBlock &diag_val) {
      add(team, rhs_val, lhs_val);
      team.team_barrier();
      divide(team, lhs_val, diag_val);
    }

    KOKKOS_INLINE_FUNCTION
    static void add_and_divide(const Vector &lhs_val, const CVector &rhs_val, const CBlock &diag_val) {
      add(rhs_val, lhs_val);
      divide(lhs_val, diag_val);
    }
  };

  /**
   * Intermediate class that contains implementation that identical
   * for blocked / non-blocked
   */
  template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, bool BlockEnabled>
  struct Intermediate : public Common<RowMapType, EntriesType, ValuesType, LHSType, RHSType, BlockEnabled> {
    using Base = Common<RowMapType, EntriesType, ValuesType, LHSType, RHSType, BlockEnabled>;

    Intermediate(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_, LHSType &lhs_,
                 const RHSType &rhs_, const entries_t &nodes_grouped_by_level_, const size_type block_size_ = 0)
        : Base(row_map_, entries_, values_, lhs_, rhs_, nodes_grouped_by_level_, block_size_) {}

    struct ReduceFunctorBasic {
      const Base *m_obj;

      KOKKOS_INLINE_FUNCTION
      ReduceFunctorBasic(const Base *obj, const lno_t = 0) : m_obj(obj) {}

      KOKKOS_INLINE_FUNCTION
      static void multiply_subtract(const scalar_t &val, const scalar_t &lhs_col_val, scalar_t &accum) {
        accum -= val * lhs_col_val;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(size_type i, scalar_t &accum) const {
        const auto colid = m_obj->entries(i);
        multiply_subtract(m_obj->vget(i), m_obj->lget(colid), accum);
      }
    };

    struct ReduceFunctorBlock : public ReduceFunctorBasic {
      using P = ReduceFunctorBasic;

      const size_type block_size;
      const size_type b;

      KOKKOS_INLINE_FUNCTION
      ReduceFunctorBlock(const Base *obj, const size_type block_size_, const size_type b_, const lno_t = 0)
          : P(obj), block_size(block_size_), b(b_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(size_type i, scalar_t &accum) const {
        const auto idx   = i / block_size;
        const auto colid = P::m_obj->entries(idx);
        P::multiply_subtract(P::m_obj->vget(idx)(b, i % block_size), P::m_obj->lget(colid)(b), accum);
      }
    };

    /**
     * If we want to support Unsorted, we'll need a Functor that returns the ptr
     * of the diag item (colid == rowid). Possibly via multi-reduce? The UnsortedTag
     * is defined above but no policies actually use it.
     */

    template <bool IsSerial, bool IsSorted, bool IsLower, bool UseThreadVec = false>
    KOKKOS_INLINE_FUNCTION void solve_impl(const member_type *team, const int my_rank, const long node_count) const {
      static_assert(!(IsSerial && UseThreadVec), "Requested thread vector range in serial?");
      static_assert(IsSorted, "Unsorted is not yet supported.");

      const auto rowid   = Base::nodes_grouped_by_level(my_rank + node_count);
      const auto soffset = Base::row_map(rowid);
      const auto eoffset = Base::row_map(rowid + 1);
      const auto rhs_val = Base::rget(rowid);

      // Set up range to auto-skip diag if is sorted
      const auto itr_b = soffset + (IsSorted ? (IsLower ? 0 : 1) : 0);
      const auto itr_e = eoffset - (IsSorted ? (IsLower ? 1 : 0) : 0);

      // We don't need the reducer to find the diag item if sorted
      typename Base::reftype lhs_val = Base::lget(rowid);

      const auto block_size_ = BlockEnabled ? Base::get_block_size() : 1;
      (void)block_size_;  // Some settings do not use this var

      if constexpr (IsSerial) {
        KK_KERNEL_ASSERT_MSG(my_rank == 0, "Non zero rank in serial");
        KK_KERNEL_ASSERT_MSG(team == nullptr, "Team provided in serial?");
        if constexpr (BlockEnabled) {
          for (size_type b = 0; b < block_size_; ++b) {
            ReduceFunctorBlock rf(this, block_size_, b, rowid);
            for (size_type i = itr_b * block_size_; i < itr_e * block_size_; ++i) {
              rf(i, lhs_val(b));
            }
          }
        } else {
          ReduceFunctorBasic rf(this, rowid);
          for (size_type i = itr_b; i < itr_e; ++i) {
            rf(i, lhs_val);
          }
        }
      } else {
        KK_KERNEL_ASSERT_MSG(team != nullptr, "Cannot do team operations without team");
        if constexpr (!UseThreadVec) {
          if constexpr (BlockEnabled) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(*team, block_size_), [&](size_type b) {
              ReduceFunctorBlock rf(this, block_size_, b, rowid);
              Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(*team, itr_b * block_size_, itr_e * block_size_), rf,
                                      lhs_val(b));
            });
          } else {
            ReduceFunctorBasic rf(this, rowid);
            Kokkos::parallel_reduce(Kokkos::TeamThreadRange(*team, itr_b, itr_e), rf, lhs_val);
          }
          team->team_barrier();
        } else {
          if constexpr (BlockEnabled) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(*team, block_size_), [&](size_type b) {
              ReduceFunctorBlock rf(this, block_size_, b, rowid);
              for (size_type i = itr_b * block_size_; i < itr_e * block_size_; ++i) {
                rf(i, lhs_val(b));
              }
            });
          } else {
            ReduceFunctorBasic rf(this, rowid);
            Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(*team, itr_b, itr_e), rf, lhs_val);
          }
        }
      }

      // If sorted, we already know the diag. Otherwise, get it from the reducer
      const lno_t diag = IsLower ? eoffset - 1 : soffset;

      // At end, handle the diag element. We need to be careful to avoid race
      // conditions here.
      if constexpr (IsSerial) {
        // Serial case is easy, there's only 1 thread so just do the
        // add_and_divide
        KK_KERNEL_ASSERT_MSG(diag != -1, "Serial should always know diag");
        Base::add_and_divide(lhs_val, rhs_val, Base::vget(diag));
      } else {
        // Parallel sorted case is complex. All threads know what the diag is.
        // If we have a team sharing the work, we need to ensure only one
        // thread performs the add_and_divide (except in BlockEnabled, then
        // we can use team operations).
        KK_KERNEL_ASSERT_MSG(diag != -1, "Sorted should always know diag");
        if constexpr (!UseThreadVec) {
          Base::add_and_divide(*team, lhs_val, rhs_val, Base::vget(diag));
        } else {
          Base::add_and_divide(lhs_val, rhs_val, Base::vget(diag));
        }
      }
    }
  };

  //
  // TriLvlSched functors
  //

  template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, bool IsLower,
            bool BlockEnabled>
  struct TriLvlSchedTP1SolverFunctor
      : public Intermediate<RowMapType, EntriesType, ValuesType, LHSType, RHSType, BlockEnabled> {
    using Base = Intermediate<RowMapType, EntriesType, ValuesType, LHSType, RHSType, BlockEnabled>;

    long node_count;  // like "block" offset into ngbl, my_league is the "local"
                      // offset

    TriLvlSchedTP1SolverFunctor(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_,
                                LHSType &lhs_, const RHSType &rhs_, const entries_t &nodes_grouped_by_level_,
                                const long &node_count_, const size_type block_size_ = 0)
        : Base(row_map_, entries_, values_, lhs_, rhs_, nodes_grouped_by_level_, block_size_),
          node_count(node_count_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type &team) const {
      Base::template solve_impl<false, true, IsLower>(&team, team.league_rank(), node_count);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const UnsortedTag &, const member_type &team) const {
      Base::template solve_impl<false, false, IsLower>(&team, team.league_rank(), node_count);
    }
  };

  template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, bool IsLower,
            bool BlockEnabled>
  struct TriLvlSchedRPSolverFunctor
      : public Intermediate<RowMapType, EntriesType, ValuesType, LHSType, RHSType, BlockEnabled> {
    using Base = Intermediate<RowMapType, EntriesType, ValuesType, LHSType, RHSType, BlockEnabled>;

    TriLvlSchedRPSolverFunctor(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_,
                               LHSType &lhs_, const RHSType &rhs_, const entries_t &nodes_grouped_by_level_,
                               const size_type block_size_ = 0)
        : Base(row_map_, entries_, values_, lhs_, rhs_, nodes_grouped_by_level_, block_size_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const lno_t i) const { Base::template solve_impl<true, true, IsLower>(nullptr, 0, i); }

    KOKKOS_INLINE_FUNCTION
    void operator()(const UnsortedTag &, const lno_t i) const {
      Base::template solve_impl<true, false, IsLower>(nullptr, 0, i);
    }
  };

  template <class RowMapType, class EntriesType, class ValuesType, class LHSType, class RHSType, bool IsLower>
  struct TriLvlSchedTP1SingleBlockFunctor
      : public Intermediate<RowMapType, EntriesType, ValuesType, LHSType, RHSType, false> {
    using Base = Intermediate<RowMapType, EntriesType, ValuesType, LHSType, RHSType, false>;

    entries_t nodes_per_level;

    long node_count;  // like "block" offset into ngbl, my_league is the "local"
                      // offset
    long lvl_start;
    long lvl_end;
    const int dense_nrows;
    const int cutoff;
    // team_size: each team can be assigned a row, if there are enough rows...

    TriLvlSchedTP1SingleBlockFunctor(const RowMapType &row_map_, const EntriesType &entries_, const ValuesType &values_,
                                     LHSType &lhs_, const RHSType &rhs_, const entries_t &nodes_grouped_by_level_,
                                     entries_t &nodes_per_level_, long node_count_, long lvl_start_, long lvl_end_,
                                     const int dense_nrows_ = 0, const int cutoff_ = 0)
        : Base(row_map_, entries_, values_, lhs_, rhs_, nodes_grouped_by_level_),
          nodes_per_level(nodes_per_level_),
          node_count(node_count_),
          lvl_start(lvl_start_),
          lvl_end(lvl_end_),
          dense_nrows(dense_nrows_),
          cutoff(cutoff_) {}

    // SingleBlock: Only one block (or league) executing; team_rank used to map
    // thread to row

    template <bool IsSorted, bool LargerCutoff>
    KOKKOS_INLINE_FUNCTION void common_impl(const member_type &team) const {
      auto mut_node_count = node_count;

      for (auto lvl = lvl_start; lvl < lvl_end; ++lvl) {
        const auto nodes_this_lvl = nodes_per_level(lvl);
        const auto my_team_rank   = team.team_rank();
        const auto loop_cutoff    = LargerCutoff ? cutoff : my_team_rank + 1;
        // If cutoff > team_size, then a thread will be responsible for multiple
        // rows - this may be a helpful scenario depending on occupancy etc.
        for (int my_rank = my_team_rank; my_rank < loop_cutoff; my_rank += team.team_size()) {
          if (my_rank < nodes_this_lvl) {
            Base::template solve_impl<false, IsSorted, IsLower, true>(&team, my_rank, mut_node_count);
          }
        }
        mut_node_count += nodes_this_lvl;
        team.team_barrier();
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type &team) const { common_impl<true, false>(team); }

    KOKKOS_INLINE_FUNCTION
    void operator()(const UnsortedTag &, const member_type &team) const { common_impl<false, false>(team); }

    KOKKOS_INLINE_FUNCTION
    void operator()(const LargerCutoffTag &, const member_type &team) const { common_impl<true, true>(team); }

    KOKKOS_INLINE_FUNCTION
    void operator()(const UnsortedLargerCutoffTag &, const member_type &team) const { common_impl<false, true>(team); }
  };

  //
  // Supernodal functors
  //

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
  // -----------------------------------------------------------
  // Helper functors for Lower-triangular solve with SpMV
  template <class LHSType>
  struct SparseTriSupernodalSpMVFunctor {
    int flag;
    long node_count;
    entries_t nodes_grouped_by_level;

    const int *supercols;
    const int *workoffset;

    LHSType X;
    work_view_t work;

    // constructor
    SparseTriSupernodalSpMVFunctor(int flag_, long node_count_, const entries_t &nodes_grouped_by_level_,
                                   const int *supercols_, const int *workoffset_, LHSType &X_, work_view_t work_)
        : flag(flag_),
          node_count(node_count_),
          nodes_grouped_by_level(nodes_grouped_by_level_),
          supercols(supercols_),
          workoffset(workoffset_),
          X(X_),
          work(work_) {}

    // operator
    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type &team) const {
      const int league_rank = team.league_rank();  // batch id
      const int team_size   = team.team_size();
      const int team_rank   = team.team_rank();
      const scalar_t zero(0.0);

      auto s = nodes_grouped_by_level(node_count + league_rank);

      // copy vector elements for the diagonal to input vector (work)
      // and zero out the corresponding elements in output (X)
      int w1 = workoffset[s];
      int j1 = supercols[s];
      // number of columns in the s-th supernode column
      int nscol = supercols[s + 1] - j1;

      if (flag == -2) {
        // copy X to work
        for (int j = team_rank; j < nscol; j += team_size) {
          work(w1 + j) = X(j1 + j);
        }
      } else if (flag == -1) {
        // copy work to X
        for (int j = team_rank; j < nscol; j += team_size) {
          X(j1 + j) = work(w1 + j);
        }
      } else if (flag == 1) {
        for (int j = team_rank; j < nscol; j += team_size) {
          work(w1 + j) = X(j1 + j);
          X(j1 + j)    = zero;
        }
      } else {
        // reinitialize work to zero
        for (int j = team_rank; j < nscol; j += team_size) {
          work(w1 + j) = zero;
        }
      }
      team.team_barrier();
    }
  };

  // -----------------------------------------------------------
  // Functor for Lower-triangular solve
  template <class ColptrView, class RowindType, class ValuesType, class LHSType>
  struct LowerTriSupernodalFunctor {
    const bool unit_diagonal;
    const bool invert_diagonal;
    const bool invert_offdiagonal;
    const int *supercols;
    ColptrView colptr;
    RowindType rowind;
    ValuesType values;

    int level;
    work_view_int_t kernel_type;
    work_view_int_t diag_kernel_type;

    LHSType X;

    work_view_t work;  // needed with gemv for update&scatter
    work_view_int_t work_offset;

    entries_t nodes_grouped_by_level;

    long node_count;

    // constructor
    LowerTriSupernodalFunctor(  // supernode info
        const bool unit_diagonal_, const bool invert_diagonal_, const bool invert_offdiagonal_, const int *supercols_,
        // L in CSC
        const ColptrView &colptr_, const RowindType &rowind_, const ValuesType &values_,
        // options to pick kernel type
        int level_, work_view_int_t &kernel_type_, work_view_int_t &diag_kernel_type_,
        // right-hand-side (input), solution (output)
        LHSType &X_,
        // workspace
        work_view_t work_, work_view_int_t &work_offset_,
        //
        const entries_t &nodes_grouped_by_level_, long node_count_)
        : unit_diagonal(unit_diagonal_),
          invert_diagonal(invert_diagonal_),
          invert_offdiagonal(invert_offdiagonal_),
          supercols(supercols_),
          colptr(colptr_),
          rowind(rowind_),
          values(values_),
          level(level_),
          kernel_type(kernel_type_),
          diag_kernel_type(diag_kernel_type_),
          X(X_),
          work(work_),
          work_offset(work_offset_),
          nodes_grouped_by_level(nodes_grouped_by_level_),
          node_count(node_count_) {}

    // operator
    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type &team) const {
      /* ----------------------------------------------------------------------
       */
      /* get inputs */
      /* ----------------------------------------------------------------------
       */
      const int league_rank = team.league_rank();  // batch id
      const int team_size   = team.team_size();
      const int team_rank   = team.team_rank();
      const scalar_t zero(0.0);
      const scalar_t one(1.0);

      auto s = nodes_grouped_by_level(node_count + league_rank);

      // supernodal column size
      const int j1 = supercols[s];
      const int j2 = supercols[s + 1];
      // > number of columns in the s-th supernode column
      const int nscol = j2 - j1;
      // "total" number of rows in all the supernodes (diagonal+off-diagonal)
      const int i1    = colptr(j1);
      const int nsrow = colptr(j1 + 1) - i1;

      // create a view for the s-th supernocal column
      // NOTE: we currently supports only default_layout = LayoutLeft
      scalar_t *dataL = const_cast<scalar_t *>(values.data());
      Kokkos::View<scalar_t **, default_layout, temp_mem_space, Kokkos::MemoryUnmanaged> viewL(&dataL[i1], nsrow,
                                                                                               nscol);

      // extract part of the solution, corresponding to the diagonal block
      auto Xj = Kokkos::subview(X, range_type(j1, j2));

      // workspace
      const int workoffset = work_offset(s);
      auto Z               = Kokkos::subview(work, range_type(workoffset + nscol, workoffset + nsrow));

      if (diag_kernel_type(level) != 3) {  // not a device-level TRSM-solve
        if (invert_offdiagonal) {
          // combined TRSM solve with diagonal + GEMV update with off-diagonal
          auto Y   = Kokkos::subview(work, range_type(workoffset,
                                                      workoffset + nsrow));  // needed for gemv instead of trmv/trsv
          auto Ljj = Kokkos::subview(viewL, range_type(0, nsrow), Kokkos::ALL());
          KokkosBlas::TeamGemv<member_type, KokkosBlas::Trans::NoTranspose, KokkosBlas::Algo::Gemv::Unblocked>::invoke(
              team, one, Ljj, Xj, zero, Y);
          team.team_barrier();
          for (int ii = team_rank; ii < nscol; ii += team_size) {
            Xj(ii) = Y(ii);
          }
          team.team_barrier();
        } else {
          /* TRSM with diagonal block */
          // extract diagonal and off-diagonal blocks of L
          auto Ljj = Kokkos::subview(viewL, range_type(0, nscol), Kokkos::ALL());
          if (invert_diagonal) {
            // workspace
            auto Y = Kokkos::subview(work, range_type(workoffset,
                                                      workoffset + nscol));  // needed for gemv instead of trmv/trsv
            for (int ii = team_rank; ii < nscol; ii += team_size) {
              Y(ii) = Xj(ii);
            }
            team.team_barrier();
            // calling team-level "Unblocked" gemv on small-size diagonal in
            // KokkosBatched
            KokkosBlas::TeamGemv<member_type, KokkosBlas::Trans::NoTranspose,
                                 KokkosBlas::Algo::Gemv::Unblocked>::invoke(team, one, Ljj, Y, zero, Xj);
          } else {
            // NOTE: we currently supports only default_layout = LayoutLeft
            Kokkos::View<scalar_t **, default_layout, temp_mem_space, Kokkos::MemoryUnmanaged> Xjj(Xj.data(), nscol, 1);
            if (unit_diagonal) {
              KokkosBatched::TeamTrsm<member_type, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                      KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::Unit,
                                      KokkosBatched::Algo::Trsm::Unblocked>::invoke(team, one, Ljj, Xjj);
            } else {
              KokkosBatched::TeamTrsm<member_type, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                      KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                      KokkosBatched::Algo::Trsm::Unblocked>::invoke(team, one, Ljj, Xjj);
            }
          }
          team.team_barrier();

          /* GEMM to update with off diagonal blocks */
          auto Lij = Kokkos::subview(viewL, range_type(nscol, nsrow), Kokkos::ALL());
          KokkosBlas::TeamGemv<member_type, KokkosBatched::Trans::NoTranspose,
                               KokkosBlas::Algo::Gemv::Unblocked>::invoke(team, one, Lij, Xj, zero, Z);
          team.team_barrier();
        }
      }

      /* scatter vectors back into X */
      int i2     = i1 + nscol;     // offset into rowind
      int nsrow2 = nsrow - nscol;  // "total" number of rows in all the off-diagonal supernodes
      Kokkos::View<scalar_t *, temp_mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>> Xatomic(
          X.data(), X.extent(0));
      for (int ii = team_rank; ii < nsrow2; ii += team_size) {
        int i = rowind(i2 + ii);
        Xatomic(i) -= Z(ii);
      }
      team.team_barrier();
    }
  };

  // -----------------------------------------------------------
  // Functor for Upper-triangular solve in CSR
  template <class ColptrType, class RowindType, class ValuesType, class LHSType>
  struct UpperTriSupernodalFunctor {
    // NOTE: we currently supports only default_layout = LayoutLeft
    using SupernodeView = typename Kokkos::View<scalar_t **, default_layout, temp_mem_space, Kokkos::MemoryUnmanaged>;

    bool invert_diagonal;
    const int *supercols;
    ColptrType colptr;
    RowindType rowind;
    ValuesType values;

    int level;
    work_view_int_t kernel_type;
    work_view_int_t diag_kernel_type;

    LHSType X;

    work_view_t work;  // needed with gemv for update&scatter
    work_view_int_t work_offset;

    entries_t nodes_grouped_by_level;

    long node_count;

    // constructor
    UpperTriSupernodalFunctor(  // supernode info
        bool invert_diagonal_, const int *supercols_,
        // U in CSR
        const ColptrType &colptr_, const RowindType &rowind_, const ValuesType &values_,
        // options to pick kernel type
        int level_, work_view_int_t &kernel_type_, work_view_int_t &diag_kernel_type_,
        // right-hand-side (input), solution (output)
        LHSType &X_,
        // workspace
        work_view_t &work_, work_view_int_t &work_offset_,
        //
        const entries_t &nodes_grouped_by_level_, long node_count_)
        : invert_diagonal(invert_diagonal_),
          supercols(supercols_),
          colptr(colptr_),
          rowind(rowind_),
          values(values_),
          level(level_),
          kernel_type(kernel_type_),
          diag_kernel_type(diag_kernel_type_),
          X(X_),
          work(work_),
          work_offset(work_offset_),
          nodes_grouped_by_level(nodes_grouped_by_level_),
          node_count(node_count_) {}

    // operator
    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type &team) const {
      /* ----------------------------------------------------------------------
       */
      /* get inputs */
      /* ----------------------------------------------------------------------
       */
      const int league_rank = team.league_rank();  // batch id
      const int team_size   = team.team_size();
      const int team_rank   = team.team_rank();
      const scalar_t zero(0.0);
      const scalar_t one(1.0);

      auto s = nodes_grouped_by_level(node_count + league_rank);

      // number of columns in the s-th supernode column
      int j1    = supercols[s];
      int j2    = supercols[s + 1];
      int nscol = j2 - j1;
      // "total" number of rows in all the supernodes (diagonal+off-diagonal)
      int i1    = colptr(j1);
      int nsrow = colptr(j1 + 1) - i1;

      // create a view of the s-th supernocal row of U
      scalar_t *dataU = const_cast<scalar_t *>(values.data());
      SupernodeView viewU(&dataU[i1], nsrow, nscol);

      // extract part of solution, corresponding to the diagonal block U(s, s)
      auto Xj       = Kokkos::subview(X, range_type(j1, j2));
      using Xj_type = decltype(Xj);

      // workspaces
      int workoffset = work_offset(s);

      // "total" number of rows in all the off-diagonal supernodes
      int nsrow2 = nsrow - nscol;
      /* gather vector into Z */
      int i2       = i1 + nscol;  // offset into rowind
      auto Z       = Kokkos::subview(work, range_type(workoffset + nscol,
                                                      workoffset + nsrow));  // needed with gemv for update&scatter
      using Z_type = decltype(Z);
      for (int ii = team_rank; ii < nsrow2; ii += team_size) {
        int i = rowind(i2 + ii);
        Z(ii) = X(i);
      }
      team.team_barrier();
      /* GEMM to update with off diagonal blocks, Xj = -Uij^T * Z */
      if (diag_kernel_type(level) != 3) {
        // not device-level GEMV-udpate
        auto Uij       = Kokkos::subview(viewU, range_type(nscol, nsrow), Kokkos::ALL());
        using Uij_type = decltype(Uij);
        KokkosBlas::TeamGemv<member_type, KokkosBatched::Trans::Transpose, KokkosBlas::Algo::Gemv::Unblocked>::
            template invoke<const scalar_t, Uij_type, Z_type, Xj_type>(team, -one, Uij, Z, one, Xj);
        team.team_barrier();

        /* TRSM with diagonal block */
        // extract diagonal and off-diagonal blocks of U
        auto Ujj       = Kokkos::subview(viewU, range_type(0, nscol), Kokkos::ALL());
        using Ujj_type = decltype(Ujj);

        if (invert_diagonal) {
          // workspace
          auto Y       = Kokkos::subview(work, range_type(workoffset,
                                                          workoffset + nscol));  // needed for gemv instead of trmv/trsv
          using Y_type = decltype(Y);
          for (int ii = team_rank; ii < nscol; ii += team_size) {
            Y(ii) = Xj(ii);
          }
          team.team_barrier();

          // caling team-level kernel in KokkosBatched on a small-size diagonal
          KokkosBlas::TeamGemv<member_type, KokkosBatched::Trans::Transpose, KokkosBlas::Algo::Gemv::Unblocked>::
              template invoke<const scalar_t, Ujj_type, Y_type, Xj_type>(team, one, Ujj, Y, zero, Xj);
        } else {
          // NOTE: we currently supports only default_layout = LayoutLeft
          Kokkos::View<scalar_t **, default_layout, temp_mem_space, Kokkos::MemoryUnmanaged> Xjj(Xj.data(), nscol, 1);
          KokkosBatched::TeamTrsm<member_type, KokkosBatched::Side::Left, KokkosBatched::Uplo::Lower,
                                  KokkosBatched::Trans::Transpose, KokkosBatched::Diag::NonUnit,
                                  KokkosBatched::Algo::Trsm::Unblocked>::invoke(team, one, Ujj, Xjj);
        }
        team.team_barrier();
      }
    }
  };

  // -----------------------------------------------------------
  // Functor for Upper-triangular solve in CSC
  template <class ColptrType, class RowindType, class ValuesType, class LHSType>
  struct UpperTriTranSupernodalFunctor {
    const bool invert_diagonal;
    const bool invert_offdiagonal;
    const int *supercols;
    ColptrType colptr;
    RowindType rowind;
    ValuesType values;

    int level;
    work_view_int_t kernel_type;
    work_view_int_t diag_kernel_type;

    LHSType X;

    work_view_t work;  // needed with gemv for update&scatter
    work_view_int_t work_offset;

    entries_t nodes_grouped_by_level;

    long node_count;

    // constructor
    UpperTriTranSupernodalFunctor(  // supernode info
        const bool invert_diagonal_, const bool invert_offdiagonal_, const int *supercols_,

        // U in CSC
        const ColptrType &colptr_, const RowindType &rowind_, const ValuesType &values_,
        // options to pick kernel type
        const int level_, const work_view_int_t &kernel_type_, const work_view_int_t &diag_kernel_type_,
        // right-hand-side (input), solution (output)
        const LHSType &X_,
        // workspace
        const work_view_t &work_, const work_view_int_t &work_offset_,
        //
        const entries_t &nodes_grouped_by_level_, const long node_count_)
        : invert_diagonal(invert_diagonal_),
          invert_offdiagonal(invert_offdiagonal_),
          supercols(supercols_),
          colptr(colptr_),
          rowind(rowind_),
          values(values_),
          level(level_),
          kernel_type(kernel_type_),
          diag_kernel_type(diag_kernel_type_),
          X(X_),
          work(work_),
          work_offset(work_offset_),
          nodes_grouped_by_level(nodes_grouped_by_level_),
          node_count(node_count_) {}

    // operator
    KOKKOS_INLINE_FUNCTION
    void operator()(const member_type &team) const {
      /* ----------------------------------------------------------------------
       */
      /* get inputs */
      /* ----------------------------------------------------------------------
       */
      const int league_rank = team.league_rank();  // batch id
      const int team_size   = team.team_size();
      const int team_rank   = team.team_rank();
      const scalar_t zero(0.0);
      const scalar_t one(1.0);

      auto s = nodes_grouped_by_level(node_count + league_rank);

      // number of columns in the s-th supernode column
      const int j1    = supercols[s];
      const int j2    = supercols[s + 1];
      const int nscol = j2 - j1;
      // "total" number of rows in all the supernodes (diagonal+off-diagonal)
      const int i1    = colptr(j1);
      const int nsrow = colptr(j1 + 1) - i1;
      // "total" number of rows in all the off-diagonal supernodes
      const int nsrow2 = nsrow - nscol;

      // create a view of the s-th supernocal column of U
      // NOTE: we currently supports only default_layout = LayoutLeft
      scalar_t *dataU = const_cast<scalar_t *>(values.data());
      Kokkos::View<scalar_t **, default_layout, temp_mem_space, Kokkos::MemoryUnmanaged> viewU(&dataU[i1], nsrow,
                                                                                               nscol);

      // extract part of solution, corresponding to the diagonal block U(s, s)
      auto Xj = Kokkos::subview(X, range_type(j1, j2));

      // workspaces
      int workoffset = work_offset(s);

      /* TRSM with diagonal block */
      if (diag_kernel_type(level) != 3) {
        // not device-level TRSM-solve
        if (invert_offdiagonal) {
          // extract diagonal + off-diagonal blocks of U
          auto Y   = Kokkos::subview(work, range_type(workoffset,
                                                      workoffset + nsrow));  // needed with gemv for update&scatter
          auto Uij = Kokkos::subview(viewU, range_type(0, nsrow), Kokkos::ALL());
          KokkosBlas::TeamGemv<member_type, KokkosBatched::Trans::NoTranspose,
                               KokkosBlas::Algo::Gemv::Unblocked>::invoke(team, one, Uij, Xj, zero, Y);
          team.team_barrier();
          // copy the diagonal back to output
          for (int ii = team_rank; ii < nscol; ii += team_size) {
            Xj(ii) = Y(ii);
          }
        } else {
          // extract diagonal block of U (stored on top)
          auto Ujj = Kokkos::subview(viewU, range_type(0, nscol), Kokkos::ALL());
          if (invert_diagonal) {
            auto Y = Kokkos::subview(work, range_type(workoffset,
                                                      workoffset + nscol));  // needed for gemv instead of trmv/trsv
            for (int ii = team_rank; ii < nscol; ii += team_size) {
              Y(ii) = Xj(ii);
            }
            team.team_barrier();
            KokkosBlas::TeamGemv<member_type, KokkosBatched::Trans::NoTranspose,
                                 KokkosBlas::Algo::Gemv::Unblocked>::invoke(team, one, Ujj, Y, zero, Xj);
          } else {
            // NOTE: we currently supports only default_layout = LayoutLeft
            Kokkos::View<scalar_t **, default_layout, temp_mem_space, Kokkos::MemoryUnmanaged> Xjj(Xj.data(), nscol, 1);
            KokkosBatched::TeamTrsm<member_type, KokkosBatched::Side::Left, KokkosBatched::Uplo::Upper,
                                    KokkosBatched::Trans::NoTranspose, KokkosBatched::Diag::NonUnit,
                                    KokkosBatched::Algo::Trsm::Unblocked>::invoke(team, one, Ujj, Xjj);
          }
        }
        team.team_barrier();
      }
      if (nsrow2 > 0) {
        /* GEMM to update off diagonal blocks, Z = Uij * Xj */
        auto Z = Kokkos::subview(work, range_type(workoffset + nscol, workoffset + nsrow));
        if (!invert_offdiagonal && diag_kernel_type(level) != 3) {
          // not device-level TRSM-solve
          auto Uij = Kokkos::subview(viewU, range_type(nscol, nsrow), Kokkos::ALL());
          KokkosBlas::TeamGemv<member_type, KokkosBatched::Trans::NoTranspose,
                               KokkosBlas::Algo::Gemv::Unblocked>::invoke(team, one, Uij, Xj, zero, Z);
          team.team_barrier();
        }

        /* scatter vector into Z */
        int i2 = i1 + nscol;  // offset into rowind
        Kokkos::View<scalar_t *, temp_mem_space, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>> Xatomic(
            X.data(), X.extent(0));
        for (int ii = team_rank; ii < nsrow2; ii += team_size) {
          int i = rowind(i2 + ii);
          Xatomic(i) -= Z(ii);
        }
        team.team_barrier();
      }
    }
  };
#endif

  //
  // End of functors, begin external API
  //

#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
  template <bool IsLower, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType>
  static void tri_solve_cg(TriSolveHandle &thandle, const RowMapType row_map, const EntriesType entries,
                           const ValuesType values, const RHSType &rhs, LHSType &lhs) {
    typename TriSolveHandle::SPTRSVcudaGraphWrapperType *lcl_cudagraph = thandle.get_sptrsvCudaGraph();

    auto nlevels = thandle.get_num_levels();

    auto stream1 = lcl_cudagraph->stream;
    Kokkos::Cuda cuda1(stream1);
    auto graph = lcl_cudagraph->cudagraph;

    Kokkos::parallel_for("Init", Kokkos::RangePolicy<execution_space>(0, 1), EmptyFunctor());
    Kokkos::Cuda().fence();
    cudaStreamSynchronize(stream1);

    auto hnodes_per_level       = thandle.get_host_nodes_per_level();
    auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();

    size_type node_count = 0;

    int team_size = thandle.get_team_size();
    team_size     = team_size == -1 ? 64 : team_size;

    // Start capturing stream
    if (thandle.cudagraphCreated == false) {
      Kokkos::fence();
      cudaStreamBeginCapture(stream1, cudaStreamCaptureModeGlobal);
      {
        for (int iter = 0; iter < nlevels; ++iter) {
          size_type lvl_nodes = hnodes_per_level(iter);

          auto policy = std::is_same<execution_space, Kokkos::Cuda>::value ? team_policy(lvl_nodes, team_size, cuda1)
                                                                           : team_policy(lvl_nodes, team_size);

          Kokkos::parallel_for(
              "parfor_l_team_cudagraph",
              Kokkos::Experimental::require(policy, Kokkos::Experimental::WorkItemProperty::HintLightWeight),
              TriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, IsLower>(
                  row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count));

          node_count += hnodes_per_level(iter);
        }
      }
      cudaStreamEndCapture(stream1, &graph);

      // Create graphExec
      cudaGraphInstantiate(&(lcl_cudagraph->cudagraphinstance), graph, NULL, NULL, 0);
      thandle.cudagraphCreated = true;
    }
    // Run graph
    Kokkos::fence();
    cudaGraphLaunch(lcl_cudagraph->cudagraphinstance, stream1);

    cudaStreamSynchronize(stream1);
    Kokkos::fence();
  }  // end tri_solve_cg

#endif

#define FunctorTypeMacro(Functor, IsLower, BlockEnabled) \
  Functor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, IsLower, BlockEnabled>

  template <bool BlockEnabled, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType>
  static void lower_tri_solve(execution_space &space, TriSolveHandle &thandle, const RowMapType row_map,
                              const EntriesType entries, const ValuesType values, const RHSType &rhs, LHSType &lhs) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
    cudaProfilerStop();
#endif
    const auto nlevels = thandle.get_num_levels();
    // Keep this a host View, create device version and copy to back to host
    // during scheduling This requires making sure the host view in the handle
    // is properly updated after the symbolic phase
    const auto nodes_per_level        = thandle.get_nodes_per_level();
    const auto hnodes_per_level       = thandle.get_host_nodes_per_level();
    const auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
    const auto block_size             = thandle.get_block_size();
    const auto block_enabled          = thandle.is_block_enabled();
    assert(block_enabled == BlockEnabled);

    // Set up functor types
    using LowerRPFunc = FunctorTypeMacro(TriLvlSchedRPSolverFunctor, true, BlockEnabled);
    using LowerTPFunc = FunctorTypeMacro(TriLvlSchedTP1SolverFunctor, true, BlockEnabled);

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    using namespace KokkosSparse::Experimental;
    using device_t            = Kokkos::Device<execution_space, temp_mem_space>;
    using integer_view_host_t = typename TriSolveHandle::integer_view_host_t;
    using row_map_host_view_t = Kokkos::View<size_type *, Kokkos::HostSpace>;

    row_map_host_view_t row_map_host;

    const scalar_t zero(0.0);
    const scalar_t one(1.0);

    auto nodes_grouped_by_level_host = thandle.get_host_nodes_grouped_by_level();

    if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE ||
        thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_ETREE ||
        thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG) {
      Kokkos::deep_copy(nodes_grouped_by_level_host, nodes_grouped_by_level);

      row_map_host =
          row_map_host_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "host rowmap"), row_map.extent(0));
      Kokkos::deep_copy(row_map_host, row_map);
    }

    // inversion options
    const bool invert_diagonal    = thandle.get_invert_diagonal();
    const bool invert_offdiagonal = thandle.get_invert_offdiagonal();
    const bool unit_diagonal      = thandle.is_unit_diagonal();

    // supernode sizes
    const int *supercols      = thandle.get_supercols();
    const int *supercols_host = thandle.get_supercols_host();

    // kernel types
    work_view_int_t kernel_type      = thandle.get_kernel_type();
    work_view_int_t diag_kernel_type = thandle.get_diag_kernel_type();

    integer_view_host_t kernel_type_host      = thandle.get_kernel_type_host();
    integer_view_host_t diag_kernel_type_host = thandle.get_diag_kernel_type_host();

    // workspaces
    work_view_int_t work_offset          = thandle.get_work_offset();
    integer_view_host_t work_offset_host = thandle.get_work_offset_host();
    auto work                            = thandle.get_workspace();
#endif

    size_type node_count = 0;

#ifdef profile_supernodal_etree
    Kokkos::Timer sptrsv_timer;
    sptrsv_timer.reset();
#endif

    for (size_type lvl = 0; lvl < nlevels; ++lvl) {
      size_type lvl_nodes = hnodes_per_level(lvl);

      if (lvl_nodes != 0) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
        cudaProfilerStart();
#endif
        if (thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_RP) {
          LowerRPFunc lrpp(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, block_size);

          Kokkos::parallel_for("parfor_fixed_lvl",
                               Kokkos::Experimental::require(range_policy(space, node_count, node_count + lvl_nodes),
                                                             Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                               lrpp);
        } else if (thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1) {
          LowerTPFunc ltpp(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count, block_size);
          int team_size = thandle.get_team_size();
          auto tp =
              team_size == -1 ? team_policy(space, lvl_nodes, Kokkos::AUTO) : team_policy(space, lvl_nodes, team_size);
          const int scratch_size = LowerTPFunc::SBlock::shmem_size(block_size, block_size);
          tp                     = tp.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for(
              "parfor_l_team",
              Kokkos::Experimental::require(tp, Kokkos::Experimental::WorkItemProperty::HintLightWeight), ltpp);
        }
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
        else if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE ||
                 thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_ETREE ||
                 thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG) {
          KK_REQUIRE_MSG(!block_enabled, "Block matrices not yet supported for supernodal");

#ifdef profile_supernodal_etree
          size_t flops = 0;
          Kokkos::Timer timer;
          timer.reset();
#endif

          // NOTE: we currently supports only default_layout = LayoutLeft
          using supernode_view_type = Kokkos::View<scalar_t **, default_layout, device_t, Kokkos::MemoryUnmanaged>;
          if (diag_kernel_type_host(lvl) == 3) {
            // using device-level kernels (functor is called to scatter the
            // results)
            scalar_t *dataL = const_cast<scalar_t *>(values.data());

            if (invert_diagonal && !invert_offdiagonal) {
              // copy diagonals to workspaces
              const int *work_offset_data = work_offset.data();
              SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(-2, node_count, nodes_grouped_by_level,
                                                                          supercols, work_offset_data, lhs, work);
              Kokkos::parallel_for(
                  "parfor_tri_supernode_spmv",
                  Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  sptrsv_init_functor);
            }

            for (size_type league_rank = 0; league_rank < lvl_nodes; league_rank++) {
              auto s = nodes_grouped_by_level_host(node_count + league_rank);

              // supernodal column size
              int j1 = supercols_host[s];
              int j2 = supercols_host[s + 1];
              // number of columns in the s-th supernode column
              int nscol = j2 - j1;
              // "total" number of rows in all the supernodes
              // (diagonal+off-diagonal)
              int i1    = row_map_host(j1);
              int nsrow = row_map_host(j1 + 1) - i1;
#ifdef profile_supernodal_etree
              flops += 2 * (nscol * nsrow);
#endif

              // workspace  (needed for gemv instead of trmv/trsv)
              int workoffset = work_offset_host(s);

              // create a view for the s-th supernocal block column
              supernode_view_type viewL(&dataL[i1], nsrow, nscol);

              // "triangular-solve" to compute Xj
              if (invert_offdiagonal) {
                auto Y  = Kokkos::subview(work, range_type(workoffset, workoffset + nsrow));
                auto Xj = Kokkos::subview(lhs, range_type(j1, j2));  // part of the solution, corresponding
                                                                     // to the diagonal block
                auto Ljj = Kokkos::subview(viewL, range_type(0, nsrow),
                                           Kokkos::ALL());  // s-th supernocal column of L
                KokkosBlas::gemv(space, "N", one, Ljj, Xj, zero, Y);
              } else {
                auto Xj = Kokkos::subview(lhs, range_type(j1, j2));  // part of the solution, corresponding
                                                                     // to the diagonal block
                auto Ljj = Kokkos::subview(viewL, range_type(0, nscol),
                                           Kokkos::ALL());  // diagonal block of s-th
                                                            // supernocal column of L
                if (invert_diagonal) {
                  auto Y = Kokkos::subview(work, range_type(workoffset, workoffset + nscol));
                  KokkosBlas::gemv(space, "N", one, Ljj, Y, zero, Xj);
                } else {
                  char unit_diag = (unit_diagonal ? 'U' : 'N');
                  // NOTE: we currently supports only default_layout =
                  // LayoutLeft
                  Kokkos::View<scalar_t **, default_layout, device_t, Kokkos::MemoryUnmanaged> Xjj(Xj.data(), nscol, 1);
                  KokkosBlas::trsm(space, "L", "L", "N", &unit_diag, one, Ljj, Xjj);
                  // TODO: space.fence();
                  Kokkos::fence();
                }
                // update off-diagonal blocks
                int nsrow2 = nsrow - nscol;  // "total" number of rows in all
                                             // the off-diagonal supernodes
                if (nsrow2 > 0) {
                  auto Z = Kokkos::subview(work, range_type(workoffset + nscol,
                                                            workoffset + nsrow));  // workspace, needed with
                                                                                   // gemv for update&scatter
                  auto Lij = Kokkos::subview(viewL, range_type(nscol, nsrow),
                                             Kokkos::ALL());  // off-diagonal blocks of s-th supernodal
                                                              // column of L
                  KokkosBlas::gemv(space, "N", one, Lij, Xj, zero, Z);
                }
              }
            }
            if (invert_offdiagonal) {
              // copy diagonals from workspaces
              const int *work_offset_data = work_offset.data();
              SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(-1, node_count, nodes_grouped_by_level,
                                                                          supercols, work_offset_data, lhs, work);
              Kokkos::parallel_for(
                  "parfor_tri_supernode_spmv",
                  Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  sptrsv_init_functor);
            }
          }

          // launching sparse-triangular solve functor
          LowerTriSupernodalFunctor<RowMapType, EntriesType, ValuesType, LHSType> sptrsv_functor(
              unit_diagonal, invert_diagonal, invert_offdiagonal, supercols, row_map, entries, values, lvl, kernel_type,
              diag_kernel_type, lhs, work, work_offset, nodes_grouped_by_level, node_count);
          Kokkos::parallel_for("parfor_lsolve_supernode",
                               Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                             Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                               sptrsv_functor);

#ifdef profile_supernodal_etree
          Kokkos::fence();
          double time_seconds = timer.seconds();
          std::cout << " > SUPERNODAL LowerTri: " << lvl << " " << time_seconds << " flop count: " << flops
                    << " kernel-type: " << kernel_type_host(lvl) << " # of supernodes: " << lvl_nodes << std::endl;
#endif
        } else if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                   thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
          KK_REQUIRE_MSG(!block_enabled, "Block matrices not yet supported for supernodal");
#ifdef profile_supernodal_etree
          Kokkos::Timer timer;
          timer.reset();
#endif

          // initialize input & output vectors

          // update with spmv (one or two SpMV)
          bool transpose_spmv = ((!thandle.transpose_spmv() && thandle.is_column_major()) ||
                                 (thandle.transpose_spmv() && !thandle.is_column_major()));
          const char *tran    = (transpose_spmv ? "T" : "N");
          if (!invert_offdiagonal) {
            // solve with diagonals
            auto digmat = thandle.get_diagblock(lvl);
            KokkosSparse::spmv(space, tran, one, digmat, lhs, one, work);
            // copy from work to lhs corresponding to diagonal blocks
            SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(-1, node_count, nodes_grouped_by_level,
                                                                        supercols, supercols, lhs, work);
            Kokkos::parallel_for("parfor_lsolve_supernode",
                                 Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                               Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                                 sptrsv_init_functor);
          } else {
            // copy lhs corresponding to diagonal blocks to work and zero out in
            // lhs
            SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(1, node_count, nodes_grouped_by_level,
                                                                        supercols, supercols, lhs, work);
            Kokkos::parallel_for("parfor_lsolve_supernode",
                                 Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                               Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                                 sptrsv_init_functor);
          }
          // update off-diagonals (potentiall combined with solve with
          // diagonals)
          auto submat = thandle.get_submatrix(lvl);
          KokkosSparse::spmv(space, tran, one, submat, work, one, lhs);

          // reinitialize workspace
          SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_finalize_functor(0, node_count, nodes_grouped_by_level,
                                                                          supercols, supercols, lhs, work);
          Kokkos::parallel_for("parfor_lsolve_supernode",
                               Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                             Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                               sptrsv_finalize_functor);

#ifdef profile_supernodal_etree
          Kokkos::fence();
          double time_seconds = timer.seconds();
          std::cout << " > SUPERNODAL LowerTri: " << lvl << " " << time_seconds
                    << " kernel-type: " << kernel_type_host(lvl) << " # of supernodes: " << lvl_nodes << std::endl;
#endif
        }
#endif
        node_count += lvl_nodes;

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
        cudaProfilerStop();
#endif
      }  // end if
    }    // end for lvl

#ifdef profile_supernodal_etree
    Kokkos::fence();
    double sptrsv_time_seconds = sptrsv_timer.seconds();
    std::cout << " + Execution space   : " << execution_space::name() << std::endl;
    std::cout << " + Memory space      : " << temp_mem_space::name() << std::endl;
    std::cout << " + SpTrsv(lower) time: " << sptrsv_time_seconds << std::endl << std::endl;
#endif
  }  // end lower_tri_solve

  template <bool BlockEnabled, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType>
  static void upper_tri_solve(execution_space &space, TriSolveHandle &thandle, const RowMapType row_map,
                              const EntriesType entries, const ValuesType values, const RHSType &rhs, LHSType &lhs) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
    cudaProfilerStop();
#endif
    using device_t = Kokkos::Device<execution_space, temp_mem_space>;

    auto nlevels = thandle.get_num_levels();
    // Keep this a host View, create device version and copy to back to host
    // during scheduling This requires making sure the host view in the handle
    // is properly updated after the symbolic phase
    auto nodes_per_level        = thandle.get_nodes_per_level();
    auto hnodes_per_level       = thandle.get_host_nodes_per_level();
    auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();
    const auto block_size       = thandle.get_block_size();
    const auto block_enabled    = thandle.is_block_enabled();
    assert(block_enabled == BlockEnabled);

    // Set up functor types
    using UpperRPFunc = FunctorTypeMacro(TriLvlSchedRPSolverFunctor, false, BlockEnabled);
    using UpperTPFunc = FunctorTypeMacro(TriLvlSchedTP1SolverFunctor, false, BlockEnabled);

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
    using namespace KokkosSparse::Experimental;
    using integer_view_host_t = typename TriSolveHandle::integer_view_host_t;
    using row_map_host_view_t = Kokkos::View<size_type *, Kokkos::HostSpace>;

    row_map_host_view_t row_map_host;

    const scalar_t zero(0.0);
    const scalar_t one(1.0);

    auto nodes_grouped_by_level_host = thandle.get_host_nodes_grouped_by_level();

    if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE ||
        thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_ETREE ||
        thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG) {
      Kokkos::deep_copy(nodes_grouped_by_level_host, nodes_grouped_by_level);

      row_map_host =
          row_map_host_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "host rowmap"), row_map.extent(0));
      Kokkos::deep_copy(row_map_host, row_map);
    }

    // supernode sizes
    const int *supercols      = thandle.get_supercols();
    const int *supercols_host = thandle.get_supercols_host();

    // inversion option
    const bool invert_diagonal    = thandle.get_invert_diagonal();
    const bool invert_offdiagonal = thandle.get_invert_offdiagonal();

    // kernel types
    work_view_int_t kernel_type      = thandle.get_kernel_type();
    work_view_int_t diag_kernel_type = thandle.get_diag_kernel_type();

    integer_view_host_t kernel_type_host      = thandle.get_kernel_type_host();
    integer_view_host_t diag_kernel_type_host = thandle.get_diag_kernel_type_host();

    // workspace
    work_view_int_t work_offset          = thandle.get_work_offset();
    integer_view_host_t work_offset_host = thandle.get_work_offset_host();
    auto work                            = thandle.get_workspace();
#endif

    size_type node_count = 0;

    // This must stay serial; would be nice to try out Cuda's graph stuff to
    // reduce kernel launch overhead
#ifdef profile_supernodal_etree
    Kokkos::Timer sptrsv_timer;
    sptrsv_timer.reset();
#endif
    for (size_type lvl = 0; lvl < nlevels; ++lvl) {
      size_type lvl_nodes = hnodes_per_level(lvl);

      if (lvl_nodes != 0) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
        cudaProfilerStart();
#endif

        if (thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_RP) {
          UpperRPFunc urpp(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, block_size);
          Kokkos::parallel_for("parfor_fixed_lvl",
                               Kokkos::Experimental::require(range_policy(space, node_count, node_count + lvl_nodes),
                                                             Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                               urpp);
        } else if (thandle.get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1) {
          UpperTPFunc utpp(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count, block_size);
          int team_size = thandle.get_team_size();
          auto tp =
              team_size == -1 ? team_policy(space, lvl_nodes, Kokkos::AUTO) : team_policy(space, lvl_nodes, team_size);
          const int scratch_size = UpperTPFunc::SBlock::shmem_size(block_size, block_size);
          tp                     = tp.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
          Kokkos::parallel_for(
              "parfor_u_team",
              Kokkos::Experimental::require(tp, Kokkos::Experimental::WorkItemProperty::HintLightWeight), utpp);
        }
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV)
        else if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_NAIVE ||
                 thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_ETREE ||
                 thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_DAG) {
          KK_REQUIRE_MSG(!block_enabled, "Block matrices not yet supported for supernodal");

#ifdef profile_supernodal_etree
          size_t flops = 0;
          Kokkos::Timer timer;
          timer.reset();
#endif

          if (thandle.is_column_major()) {  // U stored in CSC
            if (diag_kernel_type_host(lvl) == 3) {
              // using device-level kernels (functor is called to gather the
              // input into workspace)
              scalar_t *dataU = const_cast<scalar_t *>(values.data());

              if (invert_diagonal && !invert_offdiagonal) {
                // copy diagonals to workspaces
                const int *work_offset_data = work_offset.data();
                SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(-2, node_count, nodes_grouped_by_level,
                                                                            supercols, work_offset_data, lhs, work);
                Kokkos::parallel_for(
                    "parfor_tri_supernode_spmv",
                    Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                  Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                    sptrsv_init_functor);
              }
              for (size_type league_rank = 0; league_rank < lvl_nodes; league_rank++) {
                auto s = nodes_grouped_by_level_host(node_count + league_rank);

                // supernodal column size
                int j1    = supercols_host[s];
                int j2    = supercols_host[s + 1];
                int nscol = j2 - j1;  // number of columns in the s-th supernode column

                int i1    = row_map_host(j1);
                int i2    = row_map_host(j1 + 1);
                int nsrow = i2 - i1;         // "total" number of rows in all the
                                             // supernodes (diagonal+off-diagonal)
                int nsrow2 = nsrow - nscol;  // "total" number of rows in all
                                             // the off-diagonal supernodes
#ifdef profile_supernodal_etree
                flops += 2 * (nscol * nsrow);
#endif

                // workspace
                int workoffset = work_offset_host(s);

                // create a view for the s-th supernocal block column
                // NOTE: we currently supports only default_layout = LayoutLeft
                Kokkos::View<scalar_t **, default_layout, device_t, Kokkos::MemoryUnmanaged> viewU(&dataU[i1], nsrow,
                                                                                                   nscol);

                if (invert_offdiagonal) {
                  auto Uij = Kokkos::subview(viewU, range_type(0, nsrow), Kokkos::ALL());
                  auto Xj  = Kokkos::subview(lhs, range_type(j1, j2));
                  auto Z =
                      Kokkos::subview(work, range_type(workoffset,
                                                       workoffset + nsrow));  // needed with gemv for update&scatter
                  KokkosBlas::gemv(space, "N", one, Uij, Xj, zero, Z);
                } else {
                  // extract part of the solution, corresponding to the diagonal
                  // block
                  auto Xj = Kokkos::subview(lhs, range_type(j1, j2));

                  // "triangular-solve" to compute Xj
                  // extract the diagonal block of s-th supernocal column of U
                  auto Ujj = Kokkos::subview(viewU, range_type(0, nscol), Kokkos::ALL());
                  if (invert_diagonal) {
                    auto Y = Kokkos::subview(work, range_type(workoffset,
                                                              workoffset + nscol));  // needed for gemv
                                                                                     // instead of trmv/trsv
                    KokkosBlas::gemv(space, "N", one, Ujj, Y, zero, Xj);
                  } else {
                    // NOTE: we currently supports only default_layout =
                    // LayoutLeft
                    Kokkos::View<scalar_t **, default_layout, device_t, Kokkos::MemoryUnmanaged> Xjj(Xj.data(), nscol,
                                                                                                     1);
                    KokkosBlas::trsm(space, "L", "U", "N", "N", one, Ujj, Xjj);
                  }
                  // update off-diagonal blocks
                  if (nsrow2 > 0) {
                    // extract the off-diagonal blocks of s-th supernodal column
                    // of U
                    auto Uij = Kokkos::subview(viewU, range_type(nscol, nsrow), Kokkos::ALL());
                    auto Z   = Kokkos::subview(work, range_type(workoffset + nscol,
                                                                workoffset + nscol + nsrow2));  // needed with gemv for
                                                                                                // update&scatter
                    KokkosBlas::gemv(space, "N", one, Uij, Xj, zero, Z);
                  }
                }
              }
              if (invert_offdiagonal) {
                // copy diagonals from workspaces
                const int *work_offset_data = work_offset.data();
                SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(-1, node_count, nodes_grouped_by_level,
                                                                            supercols, work_offset_data, lhs, work);
                Kokkos::parallel_for(
                    "parfor_tri_supernode_spmv",
                    Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                  Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                    sptrsv_init_functor);
              }
            }

            // launching sparse-triangular solve functor
            UpperTriTranSupernodalFunctor<RowMapType, EntriesType, ValuesType, LHSType> sptrsv_functor(
                invert_diagonal, invert_offdiagonal, supercols, row_map, entries, values, lvl, kernel_type,
                diag_kernel_type, lhs, work, work_offset, nodes_grouped_by_level, node_count);

            Kokkos::parallel_for("parfor_usolve_tran_supernode",
                                 Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                               Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                                 sptrsv_functor);
          } else {  // U stored in CSR
            // launching sparse-triangular solve functor
            UpperTriSupernodalFunctor<RowMapType, EntriesType, ValuesType, LHSType> sptrsv_functor(
                invert_diagonal, supercols, row_map, entries, values, lvl, kernel_type, diag_kernel_type, lhs, work,
                work_offset, nodes_grouped_by_level, node_count);

            Kokkos::parallel_for("parfor_usolve_supernode",
                                 Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                               Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                                 sptrsv_functor);

            if (diag_kernel_type_host(lvl) == 3) {
              // using device-level kernels (functor is called to gather the
              // input into workspace)
              scalar_t *dataU = const_cast<scalar_t *>(values.data());

              for (size_type league_rank = 0; league_rank < lvl_nodes; league_rank++) {
                auto s = nodes_grouped_by_level_host(node_count + league_rank);

                // supernodal column size
                int j1    = supercols_host[s];
                int j2    = supercols_host[s + 1];
                int nscol = j2 - j1;  // number of columns in the s-th supernode column

                // "total" number of rows in all the supernodes
                // (diagonal+off-diagonal)
                int i1    = row_map_host(j1);
                int i2    = row_map_host(j1 + 1);
                int nsrow = i2 - i1;
                // "total" number of rows in all the off-diagonal supernodes
                int nsrow2 = nsrow - nscol;

                // workspace
                int workoffset = work_offset_host(s);

                // create a view for the s-th supernocal block column
                // NOTE: we currently supports only default_layout = LayoutLeft
                Kokkos::View<scalar_t **, default_layout, device_t, Kokkos::MemoryUnmanaged> viewU(&dataU[i1], nsrow,
                                                                                                   nscol);

                // extract part of the solution, corresponding to the diagonal
                // block
                auto Xj = Kokkos::subview(lhs, range_type(j1, j2));
                auto Y  = Kokkos::subview(work, range_type(workoffset,
                                                           workoffset + nscol));  // needed for gemv instead of trmv/trsv

                // update with off-diagonal blocks
                if (nsrow2 > 0) {
                  // extract the off-diagonal blocks of s-th supernodal column
                  // of
                  // U
                  auto Uij = Kokkos::subview(viewU, range_type(nscol, nsrow), Kokkos::ALL());
                  auto Z   = Kokkos::subview(
                      work, range_type(workoffset + nscol,
                                         workoffset + nscol + nsrow2));  // needed with gemv for update&scatter
                  KokkosBlas::gemv(space, "T", -one, Uij, Z, one, Xj);
                }

                // "triangular-solve" to compute Xj
                // extract the diagonal block of s-th supernocal column of U
                auto Ujj = Kokkos::subview(viewU, range_type(0, nscol), Kokkos::ALL());
                if (invert_diagonal) {
                  KokkosBlas::gemv(space, "T", one, Ujj, Xj, zero, Y);
                } else {
                  // NOTE: we currently supports only default_layout =
                  // LayoutLeft
                  Kokkos::View<scalar_t **, default_layout, device_t, Kokkos::MemoryUnmanaged> Xjj(Xj.data(), nscol, 1);
                  KokkosBlas::trsm(space, "L", "L", "T", "N", one, Ujj, Xjj);
                }
              }
              if (invert_diagonal) {
                // copy diagonals from workspaces
                const int *work_offset_data = work_offset.data();
                SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(-1, node_count, nodes_grouped_by_level,
                                                                            supercols, work_offset_data, lhs, work);
                Kokkos::parallel_for(
                    "parfor_tri_supernode_spmv",
                    Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                  Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                    sptrsv_init_functor);
              }
            }
          }
#ifdef profile_supernodal_etree
          Kokkos::fence();
          double time_seconds = timer.seconds();
          std::cout << " > SUPERNODAL UpperTri: " << lvl << " " << time_seconds << " flop count: " << flops
                    << " kernel-type: " << kernel_type_host(lvl) << " # of supernodes: " << lvl_nodes << std::endl;
#endif
        } else if (thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
                   thandle.get_algorithm() == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
          KK_REQUIRE_MSG(!block_enabled, "Block matrices not yet supported for supernodal");

#ifdef profile_supernodal_etree
          Kokkos::Timer timer;
          timer.reset();
#endif

          // initialize input & output vectors

          // update with one, or two, spmv
          bool transpose_spmv = ((!thandle.transpose_spmv() && thandle.is_column_major()) ||
                                 (thandle.transpose_spmv() && !thandle.is_column_major()));
          const char *tran    = (transpose_spmv ? "T" : "N");
          if (!transpose_spmv) {  // U stored in CSR
            if (!invert_offdiagonal) {
              // solve with diagonals
              auto digmat = thandle.get_diagblock(lvl);
              KokkosSparse::spmv(space, tran, one, digmat, lhs, one, work);
              // copy from work to lhs corresponding to diagonal blocks
              SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(-1, node_count, nodes_grouped_by_level,
                                                                          supercols, supercols, lhs, work);
              Kokkos::parallel_for(
                  "parfor_lsolve_supernode",
                  Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  sptrsv_init_functor);
            } else {
              // zero out lhs corresponding to diagonal blocks in lhs, and copy
              // to work
              SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(1, node_count, nodes_grouped_by_level,
                                                                          supercols, supercols, lhs, work);
              Kokkos::parallel_for(
                  "parfor_lsolve_supernode",
                  Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  sptrsv_init_functor);
            }
            // update with off-diagonals (potentiall combined with diagonal
            // solves)
            auto submat = thandle.get_submatrix(lvl);
            KokkosSparse::spmv(space, tran, one, submat, work, one, lhs);
          } else {
            if (!invert_offdiagonal) {
              // zero out lhs corresponding to diagonal blocks in lhs, and copy
              // to work
              SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_init_functor(1, node_count, nodes_grouped_by_level,
                                                                          supercols, supercols, lhs, work);
              Kokkos::parallel_for(
                  "parfor_lsolve_supernode",
                  Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  sptrsv_init_functor);

              // update with off-diagonals
              auto submat = thandle.get_submatrix(lvl);
              KokkosSparse::spmv(space, tran, one, submat, lhs, one, work);

              // solve with diagonals
              auto digmat = thandle.get_diagblock(lvl);
              KokkosSparse::spmv(space, tran, one, digmat, work, one, lhs);
            } else {
              std::cout << " ** invert_offdiag with U in CSR not supported **" << std::endl;
            }
          }
          // reinitialize workspace
          SparseTriSupernodalSpMVFunctor<LHSType> sptrsv_finalize_functor(0, node_count, nodes_grouped_by_level,
                                                                          supercols, supercols, lhs, work);
          Kokkos::parallel_for("parfor_lsolve_supernode",
                               Kokkos::Experimental::require(team_policy(space, lvl_nodes, Kokkos::AUTO),
                                                             Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                               sptrsv_finalize_functor);

#ifdef profile_supernodal_etree
          Kokkos::fence();
          double time_seconds = timer.seconds();
          std::cout << " > SUPERNODAL UpperTri: " << lvl << " " << time_seconds
                    << " kernel-type: " << kernel_type_host(lvl) << " # of supernodes: " << lvl_nodes << std::endl;
#endif
        }
#endif
        node_count += lvl_nodes;

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
        cudaProfilerStop();
#endif
      }  // end if
    }    // end for lvl
#ifdef profile_supernodal_etree
    Kokkos::fence();
    double sptrsv_time_seconds = sptrsv_timer.seconds();
    std::cout << " + SpTrsv(uppper) time: " << sptrsv_time_seconds << std::endl << std::endl;
    std::cout << "  + Execution space    : " << execution_space::name() << std::endl;
    std::cout << " + Memory space       : " << temp_mem_space::name() << std::endl;
#endif

  }  // end upper_tri_solve

  template <bool IsLower, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType>
  static void tri_solve_chain(execution_space &space, TriSolveHandle &thandle, const RowMapType row_map,
                              const EntriesType entries, const ValuesType values, const RHSType &rhs, LHSType &lhs) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSPSTRSV_SOLVE_IMPL_PROFILE)
    cudaProfilerStop();
#endif
    // Algorithm is checked before this function is called
    auto h_chain_ptr            = thandle.get_host_chain_ptr();
    size_type num_chain_entries = thandle.get_num_chain_entries();

    // Keep this a host View, create device version and copy to back to host
    // during scheduling This requires making sure the host view in the handle
    // is properly updated after the symbolic phase
    auto nodes_per_level  = thandle.get_nodes_per_level();
    auto hnodes_per_level = thandle.get_host_nodes_per_level();

    auto nodes_grouped_by_level = thandle.get_nodes_grouped_by_level();

    size_type node_count = 0;

    // REFACTORED to cleanup; next, need debug and timer routines
    using large_cutoff_policy_type = Kokkos::TeamPolicy<LargerCutoffTag, execution_space>;
    using SingleBlockFunctor =
        TriLvlSchedTP1SingleBlockFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, IsLower>;

    int team_size   = thandle.get_team_size();
    int vector_size = thandle.get_vector_size() > 0 ? thandle.get_vector_size() : 1;

    auto cutoff               = thandle.get_chain_threshold();
    int team_size_singleblock = team_size;

    // Enumerate options
    // ts -1,0 | cu 0 - select default ts == 1
    // ts -1,0 | cu > 0 - select default ts; restriction: ts <= tsmax (auto)
    // ts > 0 | cu 0 - set
    // ts > 0 | cu > 0 - set
    // Controls ts,cu > 0
    // co > ts  - not all rows can be mapped to a thread - must call
    // largercutoff impl co <= ts - okay, kernel must be careful not to access
    // out-of-bounds; some threads idol
    if (team_size_singleblock <= 0 && cutoff == 0) {
      team_size_singleblock = 1;
      // If cutoff == 0, no single-block calls will be made,
      // team_size_singleblock is unimportant
    }

    for (size_type chainlink = 0; chainlink < num_chain_entries; ++chainlink) {
      size_type schain = h_chain_ptr(chainlink);
      size_type echain = h_chain_ptr(chainlink + 1);

      if (echain - schain == 1) {
        // if team_size is -1 (unset), get recommended size from Kokkos
        TriLvlSchedTP1SolverFunctor<RowMapType, EntriesType, ValuesType, LHSType, RHSType, IsLower, false> tstf(
            row_map, entries, values, lhs, rhs, nodes_grouped_by_level, node_count);
        if (team_size == -1) {
          team_size = team_policy(space, 1, 1, vector_size).team_size_recommended(tstf, Kokkos::ParallelForTag());
        }

        size_type lvl_nodes = hnodes_per_level(schain);  // lvl == echain????
        Kokkos::parallel_for("parfor_l_team_chain1",
                             Kokkos::Experimental::require(team_policy(space, lvl_nodes, team_size, vector_size),
                                                           Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                             tstf);
        node_count += lvl_nodes;

      } else {
        size_type lvl_nodes = 0;

        for (size_type i = schain; i < echain; ++i) {
          lvl_nodes += hnodes_per_level(i);
        }

        if (team_size_singleblock <= 0) {
          team_size_singleblock =
              team_policy(space, 1, 1, vector_size)
                  .team_size_recommended(SingleBlockFunctor(row_map, entries, values, lhs, rhs, nodes_grouped_by_level,
                                                            nodes_per_level, node_count, schain, echain),
                                         Kokkos::ParallelForTag());
        }

        if (cutoff <= team_size_singleblock) {
          SingleBlockFunctor tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level,
                                  node_count, schain, echain);
          Kokkos::parallel_for("parfor_l_team_chainmulti",
                               Kokkos::Experimental::require(team_policy(space, 1, team_size_singleblock, vector_size),
                                                             Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                               tstf);
        } else {
          // team_size_singleblock < cutoff => kernel must allow for a
          // block-stride internally
          SingleBlockFunctor tstf(row_map, entries, values, lhs, rhs, nodes_grouped_by_level, nodes_per_level,
                                  node_count, schain, echain, 0, cutoff);
          Kokkos::parallel_for(
              "parfor_l_team_chainmulti_cutoff",
              Kokkos::Experimental::require(large_cutoff_policy_type(1, team_size_singleblock, vector_size),
                                            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
              tstf);
        }
        node_count += lvl_nodes;
      }
      // TODO: space.fence()
      Kokkos::fence();  // TODO - is this necessary? that is, can the
      // parallel_for launch before the s/echain values have
      // been updated?
    }
  }  // end tri_solve_chain

  // --------------------------------
  // Stream interfaces
  // --------------------------------
  template <bool IsLower, class RowMapType, class EntriesType, class ValuesType, class RHSType, class LHSType>
  static void tri_solve_streams(const std::vector<execution_space> &execspace_v,
                                const std::vector<TriSolveHandle *> &thandle_v,
                                const std::vector<RowMapType> &row_map_v, const std::vector<EntriesType> &entries_v,
                                const std::vector<ValuesType> &values_v, const std::vector<RHSType> &rhs_v,
                                std::vector<LHSType> &lhs_v) {
    // NOTE: Only support SEQLVLSCHD_RP and SEQLVLSCHD_TP1 at this moment
    using nodes_per_level_type        = typename TriSolveHandle::hostspace_nnz_lno_view_t;
    using nodes_grouped_by_level_type = typename TriSolveHandle::nnz_lno_view_t;
    using RPPointFunctor              = FunctorTypeMacro(TriLvlSchedRPSolverFunctor, IsLower, false);
    using TPPointFunctor              = FunctorTypeMacro(TriLvlSchedTP1SolverFunctor, IsLower, false);

    // Create vectors for handles' data in streams
    int nstreams = execspace_v.size();
    std::vector<size_type> nlevels_v(nstreams);
    std::vector<nodes_per_level_type> hnodes_per_level_v(nstreams);
    std::vector<nodes_grouped_by_level_type> nodes_grouped_by_level_v(nstreams);
    std::vector<size_type> node_count_v(nstreams);

    // Retrieve data from handles and find max. number of levels among streams
    size_type nlevels_max = 0;
    for (int i = 0; i < nstreams; i++) {
      nlevels_v[i]                = thandle_v[i]->get_num_levels();
      hnodes_per_level_v[i]       = thandle_v[i]->get_host_nodes_per_level();
      nodes_grouped_by_level_v[i] = thandle_v[i]->get_nodes_grouped_by_level();
      node_count_v[i]             = 0;
      if (nlevels_max < nlevels_v[i]) nlevels_max = nlevels_v[i];
    }

    // Main loop must be performed sequential
    for (size_type lvl = 0; lvl < nlevels_max; lvl++) {
      // 1. Launch work on all streams
      for (int i = 0; i < nstreams; i++) {
        // Only if stream i-th still has this level
        if (lvl < nlevels_v[i]) {
          size_type lvl_nodes = hnodes_per_level_v[i](lvl);
          if (lvl_nodes != 0) {
            if (thandle_v[i]->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_RP) {
              Kokkos::parallel_for("parfor_fixed_lvl",
                                   range_policy(execspace_v[i], node_count_v[i], node_count_v[i] + lvl_nodes),
                                   RPPointFunctor(row_map_v[i], entries_v[i], values_v[i], lhs_v[i], rhs_v[i],
                                                  nodes_grouped_by_level_v[i]));
            } else if (thandle_v[i]->get_algorithm() == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1) {
              int team_size = thandle_v[i]->get_team_size();
              auto tp       = team_size == -1 ? team_policy(execspace_v[i], lvl_nodes, Kokkos::AUTO)
                                              : team_policy(execspace_v[i], lvl_nodes, team_size);
              TPPointFunctor tstf(row_map_v[i], entries_v[i], values_v[i], lhs_v[i], rhs_v[i],
                                  nodes_grouped_by_level_v[i], node_count_v[i]);
              Kokkos::parallel_for("parfor_l_team", tp, tstf);
            }
            node_count_v[i] += lvl_nodes;
          }  // end if (lvl_nodes != 0)
        }    // end if (lvl < nlevels_v[i])
      }      // end for streams
    }        // end for lvl
  }          // end tri_solve_streams

};  // struct SptrsvWrap

}  // namespace Experimental
}  // namespace Impl
}  // namespace KokkosSparse

#undef FunctorTypeMacro

#endif
