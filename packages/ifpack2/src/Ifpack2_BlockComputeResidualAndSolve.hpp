// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKCOMPUTERESAND_SOLVE_IMPL_HPP
#define IFPACK2_BLOCKCOMPUTERESAND_SOLVE_IMPL_HPP

#include "Ifpack2_BlockHelper.hpp"
#include "Ifpack2_BlockComputeResidualVector.hpp"

namespace Ifpack2::BlockHelperDetails {

template <typename ExecSpace, typename DiagOffsets, typename Rowptrs,
          typename Entries>
DiagOffsets findDiagOffsets(const Rowptrs& rowptrs, const Entries& entries,
                            size_t nrows, int blocksize) {
  DiagOffsets diag_offsets(do_not_initialize_tag("btdm.diag_offsets"), nrows);
  int err1 = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecSpace>(0, nrows),
      KOKKOS_LAMBDA(size_t row, int& err2) {
        auto rowBegin = rowptrs(row);
        auto rowEnd   = rowptrs(row + 1);
        for (size_t j = rowBegin; j < rowEnd; j++) {
          if (size_t(entries(j)) == row) {
            diag_offsets(row) = j * blocksize * blocksize;
            return;
          }
        }
        err2++;
      },
      err1);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      err1, "Ifpack2 BTD: at least one row has no diagonal entry");
  return diag_offsets;
}

template <typename MatrixType, int B>
struct ComputeResidualAndSolve_1Pass {
  using impl_type        = BlockHelperDetails::ImplType<MatrixType>;
  using node_device_type = typename impl_type::node_device_type;
  using execution_space  = typename impl_type::execution_space;
  using memory_space     = typename impl_type::memory_space;

  using local_ordinal_type = typename impl_type::local_ordinal_type;
  using size_type          = typename impl_type::size_type;
  using impl_scalar_type   = typename impl_type::impl_scalar_type;
  using magnitude_type     = typename impl_type::magnitude_type;
  /// views
  using local_ordinal_type_1d_view =
      typename impl_type::local_ordinal_type_1d_view;
  using size_type_1d_view = typename impl_type::size_type_1d_view;
  using tpetra_block_access_view_type =
      typename impl_type::tpetra_block_access_view_type;  // block crs (layout
                                                          // right)
  using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view;
  using impl_scalar_type_2d_view_tpetra =
      typename impl_type::impl_scalar_type_2d_view_tpetra;  // block multivector
                                                            // (layout left)
  using btdm_scalar_type_3d_view = typename impl_type::btdm_scalar_type_3d_view;
  using btdm_scalar_type_4d_view = typename impl_type::btdm_scalar_type_4d_view;
  using i64_3d_view              = typename impl_type::i64_3d_view;

  /// team policy member type (used in cuda)
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  // enum for max blocksize and vector length
  enum : int { max_blocksize = 32 };

 private:
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> b;
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> x;  // x_owned
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> x_remote;
  Unmanaged<impl_scalar_type_2d_view_tpetra> y;

  // AmD information
  const ConstUnmanaged<impl_scalar_type_1d_view> tpetra_values;

  // blocksize
  const local_ordinal_type blocksize_requested;

  // block offsets
  const ConstUnmanaged<i64_3d_view> A_x_offsets;
  const ConstUnmanaged<i64_3d_view> A_x_offsets_remote;

  // diagonal block inverses
  const ConstUnmanaged<btdm_scalar_type_3d_view> d_inv;

  // squared update norms
  const Unmanaged<impl_scalar_type_1d_view> W;

  impl_scalar_type damping_factor;

 public:
  ComputeResidualAndSolve_1Pass(const AmD<MatrixType>& amd,
                                const btdm_scalar_type_3d_view& d_inv_,
                                const impl_scalar_type_1d_view& W_,
                                const local_ordinal_type& blocksize_requested_,
                                const impl_scalar_type& damping_factor_)
    : tpetra_values(amd.tpetra_values)
    , blocksize_requested(blocksize_requested_)
    , A_x_offsets(amd.A_x_offsets)
    , A_x_offsets_remote(amd.A_x_offsets_remote)
    , d_inv(d_inv_)
    , W(W_)
    , damping_factor(damping_factor_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& member) const {
    const local_ordinal_type blocksize      = (B == 0 ? blocksize_requested : B);
    const local_ordinal_type rowidx         = member.league_rank();
    const local_ordinal_type row            = rowidx * blocksize;
    const local_ordinal_type num_vectors    = b.extent(1);
    const local_ordinal_type num_local_rows = d_inv.extent(0);

    const impl_scalar_type* xx;
    auto A_block_cst = ConstUnmanaged<tpetra_block_access_view_type>(
        tpetra_values.data(), blocksize, blocksize);

    // Get shared allocation for a local copy of x, residual, and A
    impl_scalar_type* local_residual = reinterpret_cast<impl_scalar_type*>(
        member.team_scratch(0).get_shmem(blocksize * sizeof(impl_scalar_type)));
    impl_scalar_type* local_Dinv_residual = reinterpret_cast<impl_scalar_type*>(
        member.team_scratch(0).get_shmem(blocksize * sizeof(impl_scalar_type)));
    impl_scalar_type* local_x =
        reinterpret_cast<impl_scalar_type*>(member.thread_scratch(0).get_shmem(
            blocksize * sizeof(impl_scalar_type)));

    magnitude_type norm = 0;
    for (local_ordinal_type col = 0; col < num_vectors; ++col) {
      if (col) member.team_barrier();
      // y -= Rx
      // Initialize accumulation arrays
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, blocksize),
                           [&](const local_ordinal_type& i) {
                             local_Dinv_residual[i] = 0;
                             local_residual[i]      = b(row + i, col);
                           });
      member.team_barrier();

      int numEntries = A_x_offsets.extent(2);

      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, 0, numEntries), [&](const int k) {
            int64_t A_offset = A_x_offsets(rowidx, 0, k);
            int64_t x_offset = A_x_offsets(rowidx, 1, k);
#if KOKKOS_VERSION >= 40799
            if (A_offset != KokkosKernels::ArithTraits<int64_t>::min()) {
#else
            if (A_offset != Kokkos::ArithTraits<int64_t>::min()) {
#endif
              A_block_cst.assign_data(tpetra_values.data() + A_offset);
              // Pull x into local memory
              int64_t remote_cutoff = blocksize * num_local_rows;
              if (x_offset >= remote_cutoff)
                xx = &x_remote(x_offset - remote_cutoff, col);
              else
                xx = &x(x_offset, col);

              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(member, blocksize),
                  [&](const local_ordinal_type& i) { local_x[i] = xx[i]; });

              // matvec on block: local_residual -= A_block_cst * local_x
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize),
                                   [&](const int k0) {
                                     impl_scalar_type val = 0;
                                     for (int k1 = 0; k1 < blocksize; k1++)
                                       val += A_block_cst(k0, k1) * local_x[k1];
                                     Kokkos::atomic_add(local_residual + k0, -val);
                                   });
            }
          });
      member.team_barrier();
      // Compute local_Dinv_residual = D^-1 * local_residual
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, blocksize),
          [&](const local_ordinal_type& k0) {
            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(member, blocksize),
                [&](const local_ordinal_type& k1, impl_scalar_type& update) {
                  update += d_inv(rowidx, k0, k1) * local_residual[k1];
                },
                local_Dinv_residual[k0]);
          });
      member.team_barrier();
      // local_Dinv_residual is fully computed. Now compute the
      // squared y update norm and update y (using damping factor).
      magnitude_type colNorm;
      Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(member, blocksize),
          [&](const local_ordinal_type& k, magnitude_type& update) {
            // Compute the change in y (assuming damping_factor == 1) for this
            // entry.
            impl_scalar_type old_y    = x(row + k, col);
            impl_scalar_type y_update = local_Dinv_residual[k] - old_y;
#if KOKKOS_VERSION >= 40799
            if constexpr (KokkosKernels::ArithTraits<impl_scalar_type>::is_complex) {
#else
            if constexpr (Kokkos::ArithTraits<impl_scalar_type>::is_complex) {
#endif
              magnitude_type ydiff =
#if KOKKOS_VERSION >= 40799
                  KokkosKernels::ArithTraits<impl_scalar_type>::abs(y_update);
#else
                  Kokkos::ArithTraits<impl_scalar_type>::abs(y_update);
#endif
              update += ydiff * ydiff;
            } else {
              update += y_update * y_update;
            }
            y(row + k, col) = old_y + damping_factor * y_update;
          },
          colNorm);
      norm += colNorm;
    }
    Kokkos::single(Kokkos::PerTeam(member), [&]() { W(rowidx) = norm; });
  }

  // Launch SinglePass version (owned + nonowned residual, plus Dinv in a single
  // kernel)
  void run(const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& b_,
           const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& x_,
           const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& x_remote_,
           const Unmanaged<impl_scalar_type_2d_view_tpetra>& y_) {
    IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN;
    IFPACK2_BLOCKHELPER_TIMER_WITH_FENCE(
        "BlockTriDi::ComputeResidualAndSolve::RunSinglePass",
        ComputeResidualAndSolve0, execution_space);

    y        = y_;
    b        = b_;
    x        = x_;
    x_remote = x_remote_;

    const local_ordinal_type blocksize = blocksize_requested;
    const local_ordinal_type nrows     = d_inv.extent(0);

    const local_ordinal_type team_size   = 8;
    const local_ordinal_type vector_size = 8;
    // team: local_residual, local_Dinv_residual
    const size_t shmem_team_size = 2 * blocksize * sizeof(impl_scalar_type);
    // thread: local_x
    const size_t shmem_thread_size = blocksize * sizeof(impl_scalar_type);
    Kokkos::TeamPolicy<execution_space> policy(nrows, team_size, vector_size);
    policy.set_scratch_size(0, Kokkos::PerTeam(shmem_team_size),
                            Kokkos::PerThread(shmem_thread_size));
    Kokkos::parallel_for("ComputeResidualAndSolve::TeamPolicy::SinglePass",
                         policy, *this);
    IFPACK2_BLOCKHELPER_PROFILER_REGION_END;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(execution_space)
  }
};

template <typename MatrixType, int B>
struct ComputeResidualAndSolve_2Pass {
  using impl_type        = BlockHelperDetails::ImplType<MatrixType>;
  using node_device_type = typename impl_type::node_device_type;
  using execution_space  = typename impl_type::execution_space;
  using memory_space     = typename impl_type::memory_space;

  using local_ordinal_type = typename impl_type::local_ordinal_type;
  using size_type          = typename impl_type::size_type;
  using impl_scalar_type   = typename impl_type::impl_scalar_type;
  using magnitude_type     = typename impl_type::magnitude_type;
  /// views
  using local_ordinal_type_1d_view =
      typename impl_type::local_ordinal_type_1d_view;
  using size_type_1d_view = typename impl_type::size_type_1d_view;
  using tpetra_block_access_view_type =
      typename impl_type::tpetra_block_access_view_type;  // block crs (layout
                                                          // right)
  using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view;
  using impl_scalar_type_2d_view_tpetra =
      typename impl_type::impl_scalar_type_2d_view_tpetra;  // block multivector
                                                            // (layout left)
  using btdm_scalar_type_3d_view = typename impl_type::btdm_scalar_type_3d_view;
  using btdm_scalar_type_4d_view = typename impl_type::btdm_scalar_type_4d_view;
  using i64_3d_view              = typename impl_type::i64_3d_view;

  /// team policy member type (used in cuda)
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  // enum for max blocksize and vector length
  enum : int { max_blocksize = 32 };

  // Tag for computing residual with owned columns only (pass 1)
  struct OwnedTag {};

  // Tag for finishing the residual with nonowned columns, and solving/norming
  // (pass 2)
  struct NonownedTag {};

 private:
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> b;
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> x;  // x_owned
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> x_remote;
  Unmanaged<impl_scalar_type_2d_view_tpetra> y;

  // AmD information
  const ConstUnmanaged<impl_scalar_type_1d_view> tpetra_values;

  // blocksize
  const local_ordinal_type blocksize_requested;

  // block offsets
  const ConstUnmanaged<i64_3d_view> A_x_offsets;
  const ConstUnmanaged<i64_3d_view> A_x_offsets_remote;

  // diagonal block inverses
  const ConstUnmanaged<btdm_scalar_type_3d_view> d_inv;

  // squared update norms
  const Unmanaged<impl_scalar_type_1d_view> W;

  impl_scalar_type damping_factor;

 public:
  ComputeResidualAndSolve_2Pass(const AmD<MatrixType>& amd,
                                const btdm_scalar_type_3d_view& d_inv_,
                                const impl_scalar_type_1d_view& W_,
                                const local_ordinal_type& blocksize_requested_,
                                const impl_scalar_type& damping_factor_)
    : tpetra_values(amd.tpetra_values)
    , blocksize_requested(blocksize_requested_)
    , A_x_offsets(amd.A_x_offsets)
    , A_x_offsets_remote(amd.A_x_offsets_remote)
    , d_inv(d_inv_)
    , W(W_)
    , damping_factor(damping_factor_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const OwnedTag, const member_type& member) const {
    const local_ordinal_type blocksize   = (B == 0 ? blocksize_requested : B);
    const local_ordinal_type rowidx      = member.league_rank();
    const local_ordinal_type row         = rowidx * blocksize;
    const local_ordinal_type num_vectors = b.extent(1);

    auto A_block_cst = ConstUnmanaged<tpetra_block_access_view_type>(
        tpetra_values.data(), blocksize, blocksize);

    // Get shared allocation for a local copy of x, Ax, and A
    impl_scalar_type* local_residual = reinterpret_cast<impl_scalar_type*>(
        member.team_scratch(0).get_shmem(blocksize * sizeof(impl_scalar_type)));
    impl_scalar_type* local_x =
        reinterpret_cast<impl_scalar_type*>(member.thread_scratch(0).get_shmem(
            blocksize * sizeof(impl_scalar_type)));

    for (local_ordinal_type col = 0; col < num_vectors; ++col) {
      if (col) member.team_barrier();
      // y -= Rx
      // Initialize accumulation arrays
      Kokkos::parallel_for(
          Kokkos::TeamVectorRange(member, blocksize),
          [&](const local_ordinal_type& i) { local_residual[i] = b(row + i, col); });
      member.team_barrier();

      int numEntries = A_x_offsets.extent(2);

      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, 0, numEntries), [&](const int k) {
            int64_t A_offset = A_x_offsets(rowidx, 0, k);
            int64_t x_offset = A_x_offsets(rowidx, 1, k);
#if KOKKOS_VERSION >= 40799
            if (A_offset != KokkosKernels::ArithTraits<int64_t>::min()) {
#else
            if (A_offset != Kokkos::ArithTraits<int64_t>::min()) {
#endif
              A_block_cst.assign_data(tpetra_values.data() + A_offset);
              // Pull x into local memory
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize),
                                   [&](const local_ordinal_type& i) {
                                     local_x[i] = x(x_offset + i, col);
                                   });

              // MatVec op Ax += A*x
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize),
                                   [&](const local_ordinal_type& k0) {
                                     impl_scalar_type val = 0;
                                     for (int k1 = 0; k1 < blocksize; k1++)
                                       val += A_block_cst(k0, k1) * local_x[k1];
                                     Kokkos::atomic_add(local_residual + k0, -val);
                                   });
            }
          });
      member.team_barrier();
      // Write back the partial residual to y
      if (member.team_rank() == 0) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize),
                             [&](const local_ordinal_type& k) {
                               y(row + k, col) = local_residual[k];
                             });
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const NonownedTag, const member_type& member) const {
    const local_ordinal_type blocksize   = (B == 0 ? blocksize_requested : B);
    const local_ordinal_type rowidx      = member.league_rank();
    const local_ordinal_type row         = rowidx * blocksize;
    const local_ordinal_type num_vectors = b.extent(1);

    auto A_block_cst = ConstUnmanaged<tpetra_block_access_view_type>(
        tpetra_values.data(), blocksize, blocksize);

    // Get shared allocation for a local copy of x, Ax, and A
    impl_scalar_type* local_residual = reinterpret_cast<impl_scalar_type*>(
        member.team_scratch(0).get_shmem(blocksize * sizeof(impl_scalar_type)));
    impl_scalar_type* local_Dinv_residual = reinterpret_cast<impl_scalar_type*>(
        member.team_scratch(0).get_shmem(blocksize * sizeof(impl_scalar_type)));
    impl_scalar_type* local_x =
        reinterpret_cast<impl_scalar_type*>(member.thread_scratch(0).get_shmem(
            blocksize * sizeof(impl_scalar_type)));

    magnitude_type norm = 0;
    for (local_ordinal_type col = 0; col < num_vectors; ++col) {
      if (col) member.team_barrier();
      // y -= Rx
      // Initialize accumulation arrays.
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, blocksize),
                           [&](const local_ordinal_type& i) {
                             local_Dinv_residual[i] = 0;
                             local_residual[i]      = y(row + i, col);
                           });
      member.team_barrier();

      int numEntries = A_x_offsets_remote.extent(2);

      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, 0, numEntries), [&](const int k) {
            int64_t A_offset = A_x_offsets_remote(rowidx, 0, k);
            int64_t x_offset = A_x_offsets_remote(rowidx, 1, k);
#if KOKKOS_VERSION >= 40799
            if (A_offset != KokkosKernels::ArithTraits<int64_t>::min()) {
#else
            if (A_offset != Kokkos::ArithTraits<int64_t>::min()) {
#endif
              A_block_cst.assign_data(tpetra_values.data() + A_offset);
              // Pull x into local memory
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize),
                                   [&](const local_ordinal_type& i) {
                                     local_x[i] = x_remote(x_offset + i, col);
                                   });

              // matvec on block: local_residual -= A_block_cst * local_x
              Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize),
                                   [&](const int k0) {
                                     impl_scalar_type val = 0;
                                     for (int k1 = 0; k1 < blocksize; k1++)
                                       val += A_block_cst(k0, k1) * local_x[k1];
                                     Kokkos::atomic_add(local_residual + k0, -val);
                                   });
            }
          });
      member.team_barrier();
      // Compute local_Dinv_residual = D^-1 * local_residual
      Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, blocksize),
          [&](const local_ordinal_type& k0) {
            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(member, blocksize),
                [&](const local_ordinal_type& k1, impl_scalar_type& update) {
                  update += d_inv(rowidx, k0, k1) * local_residual[k1];
                },
                local_Dinv_residual[k0]);
          });
      member.team_barrier();
      // local_Dinv_residual is fully computed. Now compute the
      // squared y update norm and update y (using damping factor).
      magnitude_type colNorm;
      Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(member, blocksize),
          [&](const local_ordinal_type& k, magnitude_type& update) {
            // Compute the change in y (assuming damping_factor == 1) for this
            // entry.
            impl_scalar_type old_y    = x(row + k, col);
            impl_scalar_type y_update = local_Dinv_residual[k] - old_y;
#if KOKKOS_VERSION >= 40799
            if constexpr (KokkosKernels::ArithTraits<impl_scalar_type>::is_complex) {
#else
            if constexpr (Kokkos::ArithTraits<impl_scalar_type>::is_complex) {
#endif
              magnitude_type ydiff =
#if KOKKOS_VERSION >= 40799
                  KokkosKernels::ArithTraits<impl_scalar_type>::abs(y_update);
#else
                  Kokkos::ArithTraits<impl_scalar_type>::abs(y_update);
#endif
              update += ydiff * ydiff;
            } else {
              update += y_update * y_update;
            }
            y(row + k, col) = old_y + damping_factor * y_update;
          },
          colNorm);
      norm += colNorm;
    }
    Kokkos::single(Kokkos::PerTeam(member), [&]() { W(rowidx) = norm; });
  }

  // Launch pass 1 of the 2-pass version.
  // This computes just the owned part of residual and writes that back to y.
  void run_pass1(const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& b_,
                 const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& x_,
                 const Unmanaged<impl_scalar_type_2d_view_tpetra>& y_) {
    IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN;
    IFPACK2_BLOCKHELPER_TIMER_WITH_FENCE(
        "BlockTriDi::ComputeResidualAndSolve::RunPass1",
        ComputeResidualAndSolve0, execution_space);

    b = b_;
    x = x_;
    y = y_;

    const local_ordinal_type blocksize = blocksize_requested;
    const local_ordinal_type nrows     = d_inv.extent(0);

    const local_ordinal_type team_size   = 8;
    const local_ordinal_type vector_size = 8;
    const size_t shmem_team_size         = blocksize * sizeof(impl_scalar_type);
    const size_t shmem_thread_size       = blocksize * sizeof(impl_scalar_type);
    Kokkos::TeamPolicy<execution_space, OwnedTag> policy(nrows, team_size,
                                                         vector_size);
    policy.set_scratch_size(0, Kokkos::PerTeam(shmem_team_size),
                            Kokkos::PerThread(shmem_thread_size));
    Kokkos::parallel_for("ComputeResidualAndSolve::TeamPolicy::Pass1", policy,
                         *this);
    IFPACK2_BLOCKHELPER_PROFILER_REGION_END;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(execution_space)
  }

  // Launch pass 2 of the 2-pass version.
  // This finishes computing residual with x_remote,
  // and then applies Dinv and computes norm.
  void run_pass2(
      const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& x_,
      const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& x_remote_,
      const Unmanaged<impl_scalar_type_2d_view_tpetra>& y_) {
    IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN;
    IFPACK2_BLOCKHELPER_TIMER_WITH_FENCE(
        "BlockTriDi::ComputeResidualAndSolve::RunPass2",
        ComputeResidualAndSolve0, execution_space);

    x        = x_;
    x_remote = x_remote_;
    y        = y_;

    const local_ordinal_type blocksize = blocksize_requested;
    const local_ordinal_type nrows     = d_inv.extent(0);

    const local_ordinal_type team_size   = 8;
    const local_ordinal_type vector_size = 8;
    const size_t shmem_team_size         = 2 * blocksize * sizeof(impl_scalar_type);
    const size_t shmem_thread_size       = blocksize * sizeof(impl_scalar_type);
    Kokkos::TeamPolicy<execution_space, NonownedTag> policy(nrows, team_size,
                                                            vector_size);
    policy.set_scratch_size(0, Kokkos::PerTeam(shmem_team_size),
                            Kokkos::PerThread(shmem_thread_size));
    Kokkos::parallel_for("ComputeResidualAndSolve::TeamPolicy::Pass2", policy,
                         *this);
    IFPACK2_BLOCKHELPER_PROFILER_REGION_END;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(execution_space)
  }
};

template <typename MatrixType, int B>
struct ComputeResidualAndSolve_SolveOnly {
  using impl_type        = BlockHelperDetails::ImplType<MatrixType>;
  using node_device_type = typename impl_type::node_device_type;
  using execution_space  = typename impl_type::execution_space;
  using memory_space     = typename impl_type::memory_space;

  using local_ordinal_type = typename impl_type::local_ordinal_type;
  using size_type          = typename impl_type::size_type;
  using impl_scalar_type   = typename impl_type::impl_scalar_type;
  using magnitude_type     = typename impl_type::magnitude_type;
  /// views
  using local_ordinal_type_1d_view =
      typename impl_type::local_ordinal_type_1d_view;
  using size_type_1d_view = typename impl_type::size_type_1d_view;
  using tpetra_block_access_view_type =
      typename impl_type::tpetra_block_access_view_type;  // block crs (layout
                                                          // right)
  using impl_scalar_type_1d_view = typename impl_type::impl_scalar_type_1d_view;
  using impl_scalar_type_2d_view_tpetra =
      typename impl_type::impl_scalar_type_2d_view_tpetra;  // block multivector
                                                            // (layout left)
  using btdm_scalar_type_3d_view = typename impl_type::btdm_scalar_type_3d_view;
  using btdm_scalar_type_4d_view = typename impl_type::btdm_scalar_type_4d_view;
  using i64_3d_view              = typename impl_type::i64_3d_view;

  /// team policy member type (used in cuda)
  using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

  // enum for max blocksize and vector length
  enum : int { max_blocksize = 32 };

 private:
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> b;
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> x;  // x_owned
  ConstUnmanaged<impl_scalar_type_2d_view_tpetra> x_remote;
  Unmanaged<impl_scalar_type_2d_view_tpetra> y;

  // AmD information
  const ConstUnmanaged<impl_scalar_type_1d_view> tpetra_values;

  // blocksize
  const local_ordinal_type blocksize_requested;

  // block offsets
  const ConstUnmanaged<i64_3d_view> A_x_offsets;
  const ConstUnmanaged<i64_3d_view> A_x_offsets_remote;

  // diagonal block inverses
  const ConstUnmanaged<btdm_scalar_type_3d_view> d_inv;

  // squared update norms
  const Unmanaged<impl_scalar_type_1d_view> W;

  impl_scalar_type damping_factor;

 public:
  ComputeResidualAndSolve_SolveOnly(
      const AmD<MatrixType>& amd, const btdm_scalar_type_3d_view& d_inv_,
      const impl_scalar_type_1d_view& W_,
      const local_ordinal_type& blocksize_requested_,
      const impl_scalar_type& damping_factor_)
    : tpetra_values(amd.tpetra_values)
    , blocksize_requested(blocksize_requested_)
    , A_x_offsets(amd.A_x_offsets)
    , A_x_offsets_remote(amd.A_x_offsets_remote)
    , d_inv(d_inv_)
    , W(W_)
    , damping_factor(damping_factor_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const member_type& member) const {
    const local_ordinal_type blocksize = (B == 0 ? blocksize_requested : B);
    const local_ordinal_type rowidx =
        member.league_rank() * member.team_size() + member.team_rank();
    const local_ordinal_type row         = rowidx * blocksize;
    const local_ordinal_type num_vectors = b.extent(1);

    // Get shared allocation for a local copy of x, Ax, and A
    impl_scalar_type* local_Dinv_residual =
        reinterpret_cast<impl_scalar_type*>(member.thread_scratch(0).get_shmem(
            blocksize * sizeof(impl_scalar_type)));

    if (rowidx >= (local_ordinal_type)d_inv.extent(0)) return;

    magnitude_type norm = 0;
    for (local_ordinal_type col = 0; col < num_vectors; ++col) {
      // Compute local_Dinv_residual = D^-1 * local_residual
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, blocksize),
                           [&](const local_ordinal_type& k0) {
                             impl_scalar_type val = 0;
                             for (local_ordinal_type k1 = 0; k1 < blocksize;
                                  k1++) {
                               val += d_inv(rowidx, k0, k1) * b(row + k1, col);
                             }
                             local_Dinv_residual[k0] = val;
                           });

      magnitude_type colNorm;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(member, blocksize),
          [&](const local_ordinal_type& k, magnitude_type& update) {
            // Compute the change in y (assuming damping_factor == 1) for this
            // entry.
            impl_scalar_type y_update = local_Dinv_residual[k];
#if KOKKOS_VERSION >= 40799
            if constexpr (KokkosKernels::ArithTraits<impl_scalar_type>::is_complex) {
#else
            if constexpr (Kokkos::ArithTraits<impl_scalar_type>::is_complex) {
#endif
              magnitude_type ydiff =
#if KOKKOS_VERSION >= 40799
                  KokkosKernels::ArithTraits<impl_scalar_type>::abs(y_update);
#else
                  Kokkos::ArithTraits<impl_scalar_type>::abs(y_update);
#endif
              update += ydiff * ydiff;
            } else {
              update += y_update * y_update;
            }
            y(row + k, col) = damping_factor * y_update;
          },
          colNorm);
      norm += colNorm;
    }
    Kokkos::single(Kokkos::PerThread(member), [&]() { W(rowidx) = norm; });
  }

  // ComputeResidualAndSolve_SolveOnly::run does the solve for the first
  // iteration, when the initial guess for y is zero. This means the residual
  // vector is just b. The kernel applies the inverse diags to b to find y, and
  // also puts the partial squared update norms (1 per row) into W.
  void run(const ConstUnmanaged<impl_scalar_type_2d_view_tpetra>& b_,
           const Unmanaged<impl_scalar_type_2d_view_tpetra>& y_) {
    IFPACK2_BLOCKHELPER_PROFILER_REGION_BEGIN;
    IFPACK2_BLOCKHELPER_TIMER_WITH_FENCE(
        "BlockTriDi::ComputeResidualAndSolve::Run_Y_Zero",
        ComputeResidualAndSolve0, execution_space);

    y = y_;
    b = b_;

    const local_ordinal_type blocksize = blocksize_requested;
    const local_ordinal_type nrows     = d_inv.extent(0);

    const local_ordinal_type team_size   = 8;
    const local_ordinal_type vector_size = 8;
    const size_t shmem_thread_size       = blocksize * sizeof(impl_scalar_type);
    Kokkos::TeamPolicy<execution_space> policy(
        (nrows + team_size - 1) / team_size, team_size, vector_size);
    policy.set_scratch_size(0, Kokkos::PerThread(shmem_thread_size));
    Kokkos::parallel_for("ComputeResidualAndSolve::TeamPolicy::y_zero", policy,
                         *this);
    IFPACK2_BLOCKHELPER_PROFILER_REGION_END;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(execution_space)
  }
};

}  // namespace Ifpack2::BlockHelperDetails

#endif
