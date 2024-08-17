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

#ifndef KOKKOSSPARSE_SPMV_IMPL_MERGE_HPP
#define KOKKOSSPARSE_SPMV_IMPL_MERGE_HPP

#include <sstream>

#include "KokkosKernels_Iota.hpp"
#include "KokkosKernels_AlwaysFalse.hpp"

#include "KokkosSparse_merge_matrix.hpp"

namespace KokkosSparse::Impl {

/*! \brief Merge-based SpMV
  Hierarchical GPU implementation
  Each team uses MergePath search to find the non-zeros and rows it is
  responsible for Each thread in the team similarly uses diagonal search within
  the team to determine which entries it will be responsible for
  The threads then atomically accumulate partial produces
*/
template <class ExecutionSpace, class AMatrix, class XVector, class YVector>
struct SpmvMergeHierarchical {
  using device_type                  = typename YVector::device_type;
  using exec_space                   = ExecutionSpace;
  using y_value_type                 = typename YVector::non_const_value_type;
  using x_value_type                 = typename XVector::non_const_value_type;
  using A_value_type                 = typename AMatrix::non_const_value_type;
  using A_ordinal_type               = typename AMatrix::non_const_ordinal_type;
  using A_size_type                  = typename AMatrix::non_const_size_type;
  using row_map_non_const_value_type = typename AMatrix::row_map_type::non_const_value_type;

  using policy_type = Kokkos::TeamPolicy<exec_space>;
  using team_member = typename policy_type::member_type;

  using um_row_map_type =
      Kokkos::View<typename AMatrix::row_map_type::data_type, typename AMatrix::row_map_type::device_type::memory_space,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using row_map_scratch_type =
      Kokkos::View<typename AMatrix::row_map_type::non_const_data_type, typename exec_space::scratch_memory_space,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using iota_type = KokkosKernels::Impl::Iota<A_size_type, A_size_type>;

  using DSR = typename KokkosSparse::Impl::MergeMatrixDiagonal<um_row_map_type, iota_type>::position_type;

  using KAT = Kokkos::ArithTraits<A_value_type>;

  // results of a lower-bound and upper-bound diagonal search
  struct Chunk {
    DSR lb;  // lower bound
    DSR ub;  // upper bound
  };

  template <bool NONZEROS_USE_SCRATCH, bool ROWENDS_USE_SCRATCH, bool Y_USE_SCRATCH, bool CONJ>
  struct SpmvMergeImplFunctor {
    SpmvMergeImplFunctor(const y_value_type& _alpha, const AMatrix& _A, const XVector& _x, const YVector& _y,
                         const A_size_type pathLengthThreadChunk)
        : alpha(_alpha), A(_A), x(_x), y(_y), pathLengthThreadChunk_(pathLengthThreadChunk) {}

    y_value_type alpha;
    AMatrix A;
    XVector x;
    YVector y;
    A_size_type pathLengthThreadChunk_;

    KOKKOS_INLINE_FUNCTION void operator()(const team_member& thread) const {
      const A_size_type pathLengthTeamChunk = thread.team_size() * pathLengthThreadChunk_;

      const A_size_type pathLength = A.numRows() + A.nnz();
      const A_size_type teamD      = thread.league_rank() * pathLengthTeamChunk;  // diagonal
      const A_size_type teamDEnd   = KOKKOSKERNELS_MACRO_MIN(teamD + pathLengthTeamChunk, pathLength);

      // iota(i) -> i
      iota_type iota(A.nnz());

      // remove leading 0 from row_map
      um_row_map_type rowEnds(&A.graph.row_map(1), A.graph.row_map.size() - 1);

      // compiler thinks these are "used" in team_broadcast below, so initialize
      // them with something to silence the warning
      DSR lb{};
      DSR ub{};

      // thread 0 does the lower bound, thread 1 does the upper bound
      if (0 == thread.team_rank() || 1 == thread.team_rank()) {
        const A_size_type d = thread.team_rank() ? teamDEnd : teamD;
        DSR dsr             = diagonal_search(rowEnds, iota, d);
        if (0 == thread.team_rank()) {
          lb = dsr;
        }
        if (1 == thread.team_rank()) {
          ub = dsr;
        }
      }
      thread.team_broadcast(lb, 0);
      thread.team_broadcast(ub, 1);
      const A_size_type teamNnzBegin    = lb.bi;  // the first nnz this team will handle
      const A_size_type teamNnzEnd      = ub.bi;  // one-past the last nnz this team will handle
      const A_ordinal_type teamRowBegin = lb.ai;  // <= the row than the first nnz is in
      const A_ordinal_type teamRowEnd   = ub.ai;  // >= the row than the last nnz is in

      // team-collaborative copy of matrix data into scratch
      A_size_type* rowEndsS{nullptr};
      A_ordinal_type* entriesS{nullptr};
      A_value_type* valuesS{nullptr};
      y_value_type* yS{nullptr};

      if constexpr (ROWENDS_USE_SCRATCH) {
        rowEndsS = (A_size_type*)thread.team_shmem().get_shmem(pathLengthTeamChunk * sizeof(A_size_type));

        // teamRowEnd may be equal to the row the team's last nnz is in
        // so in most cases we want to read it (teamRowEnd+1). However,
        // however, guard against reading off the end of the view
        Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, teamRowBegin, teamRowEnd + 1),
                             [&](const A_ordinal_type& i) {
                               if (i < A.numRows()) {
                                 rowEndsS[i - teamRowBegin] = rowEnds(i);
                               } else {
                                 rowEndsS[i - teamRowBegin] = A.nnz();
                               }
                             });
      } else {
        (void)(rowEndsS == rowEndsS);  // set but unused, expr has no effect
      }

      if constexpr (NONZEROS_USE_SCRATCH) {
        valuesS  = (A_value_type*)thread.team_shmem().get_shmem(pathLengthTeamChunk * sizeof(A_value_type));
        entriesS = (A_ordinal_type*)thread.team_shmem().get_shmem(pathLengthTeamChunk * sizeof(A_ordinal_type));
        Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, teamNnzBegin, teamNnzEnd), [&](const A_ordinal_type& i) {
          valuesS[i - teamNnzBegin]  = A.values(i);
          entriesS[i - teamNnzBegin] = A.graph.entries(i);
        });
      } else {
        (void)(entriesS == entriesS);  // set but unused, expr has no effect
        (void)(valuesS == valuesS);    // set but unused, expr has no effect
      }

      if constexpr (Y_USE_SCRATCH) {
        yS = (y_value_type*)thread.team_shmem().get_shmem(pathLengthTeamChunk * sizeof(y_value_type));
        Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, teamRowBegin, teamRowEnd + 1),
                             [&](const A_ordinal_type& i) {
                               if (i < A.numRows()) {
                                 yS[i - teamRowBegin] = 0;
                               }
                             });
      } else {
        (void)(yS == yS);  // set but unused, expr has no effect
      }

      if constexpr (ROWENDS_USE_SCRATCH || NONZEROS_USE_SCRATCH || Y_USE_SCRATCH) {
        thread.team_barrier();
      }

      // each thread determines its location within the team chunk

      // team's view of row map is either in scratch or global
      typename std::conditional<ROWENDS_USE_SCRATCH, row_map_scratch_type, um_row_map_type>::type teamRowEnds;
      if constexpr (ROWENDS_USE_SCRATCH) {
        teamRowEnds = row_map_scratch_type(rowEndsS, teamRowEnd - teamRowBegin);
      } else {
        teamRowEnds = um_row_map_type(&rowEnds(teamRowBegin), teamRowEnd - teamRowBegin);
      }

      iota_type teamIota(teamNnzEnd - teamNnzBegin,
                         teamNnzBegin);  // teamNnzBegin..<teamNnzEnd

      // diagonal is local to the team's path
      A_size_type threadD = KOKKOSKERNELS_MACRO_MIN(teamD + thread.team_rank() * pathLengthThreadChunk_, pathLength);
      threadD -= teamD;

      DSR threadLb;
      threadLb                            = diagonal_search(teamRowEnds, teamIota, threadD);
      const A_size_type threadNnzBegin    = threadLb.bi + teamNnzBegin;
      const A_ordinal_type threadRowBegin = threadLb.ai + teamRowBegin;

      // each thread does some accumulating
      y_value_type acc      = 0;
      A_ordinal_type curRow = threadRowBegin;
      A_size_type curNnz    = threadNnzBegin;
      for (A_size_type i = 0;
           i < pathLengthThreadChunk_ /*some threads have less work*/ && curNnz < A.nnz() + 1 && curRow < A.numRows();
           ++i) {
        A_size_type curRowEnd;
        if constexpr (ROWENDS_USE_SCRATCH) {
          curRowEnd = rowEndsS[curRow - teamRowBegin];
        } else {
          curRowEnd = rowEnds(curRow);
        }

        if (curNnz < curRowEnd) {
          typename AMatrix::non_const_ordinal_type col;
          A_value_type val;
          if constexpr (NONZEROS_USE_SCRATCH) {
            col = entriesS[curNnz - teamNnzBegin];
            val = (CONJ ? KAT::conj(valuesS[curNnz - teamNnzBegin]) : valuesS[curNnz - teamNnzBegin]);
          } else {
            col = A.graph.entries(curNnz);
            val = (CONJ ? KAT::conj(A.values(curNnz)) : A.values(curNnz));
          }

          acc += val * x(col);
          ++curNnz;
        } else {
          if constexpr (Y_USE_SCRATCH) {
            Kokkos::atomic_add(&yS[curRow - teamRowBegin], alpha * acc);
          } else {
            Kokkos::atomic_add(&y(curRow), alpha * acc);
          }
          acc = 0;
          ++curRow;
        }
      }
      // save the accumulated results of a partial last row.
      // might be 0 if last row was not partial
      if (curRow < A.numRows()) {
        if constexpr (Y_USE_SCRATCH) {
          Kokkos::atomic_add(&yS[curRow - teamRowBegin], alpha * acc);
        } else {
          Kokkos::atomic_add(&y(curRow), alpha * acc);
        }
      }

      if constexpr (Y_USE_SCRATCH) {
        thread.team_barrier();

        Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, teamRowBegin, teamRowEnd + 1),
                             [&](const A_ordinal_type& i) {
                               if (i < A.numRows()) {
                                 if (i > teamRowBegin && i < teamRowEnd) {
                                   y(i) += yS[i - teamRowBegin];
                                 } else {
                                   Kokkos::atomic_add(&y(i), yS[i - teamRowBegin]);
                                 }
                               }
                             });
      }
    }

    size_t team_shmem_size(int teamSize) const {
      const A_size_type pathLengthTeamChunk = pathLengthThreadChunk_ * teamSize;
      (void)pathLengthTeamChunk;  // silence declared but not referenced
      size_t val = 0;
      if constexpr (Y_USE_SCRATCH) {
        val += sizeof(y_value_type) * pathLengthTeamChunk;
      }
      if constexpr (ROWENDS_USE_SCRATCH) {
        val += sizeof(row_map_non_const_value_type) * pathLengthTeamChunk;
      }
      if constexpr (NONZEROS_USE_SCRATCH) {
        val += sizeof(A_ordinal_type) * pathLengthTeamChunk;
        val += sizeof(A_value_type) * pathLengthTeamChunk;
      }
      return val;
    }
  };  // struct SpmvMergeImplFunctor

  static void spmv(const ExecutionSpace& space, const char mode[], const y_value_type& alpha, const AMatrix& A,
                   const XVector& x, const y_value_type& beta, const YVector& y) {
    static_assert(XVector::rank == 1, "");
    static_assert(YVector::rank == 1, "");

    if (y_value_type(0) == beta) {
      Kokkos::deep_copy(space, y, y_value_type(0));
    } else {
      KokkosBlas::scal(space, y, beta, y);
    }

    /* determine launch parameters for different architectures
       On architectures where there is a natural execution hierarchy with true
       team scratch, we'll assign each team to use an appropriate amount of the
       scratch.
       On other architectures, just have each team do the maximal amount of work
       to amortize the cost of the diagonal search
    */
    const A_size_type pathLength = A.numRows() + A.nnz();
    A_size_type pathLengthThreadChunk;
    int teamSize;
    if constexpr (KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>()) {
      pathLengthThreadChunk = 4;
      teamSize              = 128;
    } else {
      teamSize              = 1;
      pathLengthThreadChunk = (pathLength + exec_space().concurrency() - 1) / exec_space().concurrency();
    }

    const size_t pathLengthTeamChunk = pathLengthThreadChunk * teamSize;
    const int leagueSize             = (pathLength + pathLengthTeamChunk - 1) / pathLengthTeamChunk;

    policy_type policy(space, leagueSize, teamSize);

    /* Currently:
       On GPU, assume atomics are fast, so don't accumuate into scratch.
       On CPU spaces, there's no real point to using scratch, just rely on the
       memory hierarchy. Using scratch just increases the number of required
       atomic operations
    */
    if (KokkosSparse::NoTranspose[0] == mode[0]) {
      constexpr bool CONJ = false;
      using GpuOp         = SpmvMergeImplFunctor<true, true, false, CONJ>;
      using CpuOp         = SpmvMergeImplFunctor<false, false, false, CONJ>;
      using Op =
          typename std::conditional<KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>(), GpuOp, CpuOp>::type;
      Op op(alpha, A, x, y, pathLengthThreadChunk);
      Kokkos::parallel_for("SpmvMergeHierarchical::spmv", policy, op);
    } else if (KokkosSparse::Conjugate[0] == mode[0]) {
      constexpr bool CONJ = true;
      using GpuOp         = SpmvMergeImplFunctor<true, true, false, CONJ>;
      using CpuOp         = SpmvMergeImplFunctor<false, false, false, CONJ>;
      using Op =
          typename std::conditional<KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>(), GpuOp, CpuOp>::type;
      Op op(alpha, A, x, y, pathLengthThreadChunk);
      Kokkos::parallel_for("SpmvMergeHierarchical::spmv", policy, op);
    } else {
      std::stringstream ss;
      ss << __FILE__ << ":" << __LINE__ << "SpmvMergeHierarchical::spmv() called with unsupported mode " << mode;
      throw std::logic_error(ss.str());
    }
  }
};

}  // namespace KokkosSparse::Impl

#endif  // KOKKOSSPARSE_SPMV_IMPL_MERGE_HPP
