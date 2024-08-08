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

#include <Kokkos_Core.hpp>
#include <iostream>
#include <string>

#ifndef _PAR_ILUTHANDLE_HPP
#define _PAR_ILUTHANDLE_HPP

namespace KokkosSparse {
namespace Experimental {

/**
 * Handle for par_ilut. Contains useful types, par_ilut configuration settings,
 * symbolic settings and scalar output info.
 *
 * For more info, see KokkosSparse_par_ilut.hpp doxygen
 */
template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class PAR_ILUTHandle {
 public:
  using HandleExecSpace             = ExecutionSpace;
  using HandleTempMemorySpace       = TemporaryMemorySpace;
  using HandlePersistentMemorySpace = PersistentMemorySpace;

  using execution_space = ExecutionSpace;
  using memory_space    = HandlePersistentMemorySpace;
  using TeamPolicy      = Kokkos::TeamPolicy<execution_space>;
  using RangePolicy     = Kokkos::RangePolicy<execution_space>;

  using size_type       = typename std::remove_const<size_type_>::type;
  using const_size_type = const size_type;

  using nnz_lno_t       = typename std::remove_const<lno_t_>::type;
  using const_nnz_lno_t = const nnz_lno_t;

  using nnz_scalar_t       = typename std::remove_const<scalar_t_>::type;
  using const_nnz_scalar_t = const nnz_scalar_t;

  using float_t = typename Kokkos::ArithTraits<nnz_scalar_t>::mag_type;

  using nnz_row_view_t = typename Kokkos::View<size_type *, HandlePersistentMemorySpace>;

  using nnz_lno_view_t = typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace>;

  using nnz_value_view_t = typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace>;

  using signed_integral_t = typename std::make_signed<typename nnz_row_view_t::non_const_value_type>::type;

  using signed_nnz_lno_view_t =
      Kokkos::View<signed_integral_t *, typename nnz_row_view_t::array_layout, typename nnz_row_view_t::device_type,
                   typename nnz_row_view_t::memory_traits>;

 private:
  // User inputs
  size_type max_iter;                /// Hard cap on the number of par_ilut iterations
  float_t residual_norm_delta_stop;  /// When the change in residual from
                                     /// iteration to iteration drops below
                                     /// this, the algorithm will stop (even if
                                     /// max_iters has not been hit)
  float_t fill_in_limit;             /// The threshold for removing candidates
                                     /// from the intermediate L and U is set such
                                     /// that the resulting sparsity pattern has
                                     /// at most `fill_in_limit` times the number
                                     /// of non-zeros of the ILU(0)
                                     /// factorization. This selection is executed
                                     /// separately for both factors L and U.
  bool async_update;                 /// Whether compute LU factors should do asychronous
                                     /// updates. When ON, the algorithm will usually converge
                                     /// faster but it makes the algorithm non-deterministic.
  bool verbose;                      /// Print information while executing par_ilut

  // Stored by parent KokkosKernelsHandle
  int team_size;    /// Kokkos team size. Set by the parent handle. -1 implies
                    /// AUTO
  int vector_size;  /// Kokkos vector size. Set by the parent handle.

  // Stored by symbolic phase
  size_type nrows;         /// Number of rows in the CSRs given to the symbolic par_ilut
  size_type nnzL;          /// Number of non-zero entries in the L part of A in the CSRs
                           /// given to the symbolic par_ilut
  size_type nnzU;          /// Number of non-zero entries in the U part of A in the CSRs
                           /// given to the symbolic par_ilut
  bool symbolic_complete;  /// Whether symbolic par_ilut has been called

  // Outputs
  int num_iters;             /// The number of iterations par_ilut took to finish
  nnz_scalar_t end_rel_res;  /// The A - LU residual norm at the time the
                             /// algorithm finished

 public:
  // See KokkosKernelsHandle::create_par_ilut_handle for default user input
  // values
  PAR_ILUTHandle(const size_type max_iter_, const float_t residual_norm_delta_stop_, const float_t fill_in_limit_,
                 const bool async_update_, const bool verbose_)
      : max_iter(max_iter_),
        residual_norm_delta_stop(residual_norm_delta_stop_),
        fill_in_limit(fill_in_limit_),
        async_update(async_update_),
        verbose(verbose_),
        team_size(-1),
        vector_size(-1),
        nrows(0),
        nnzL(0),
        nnzU(0),
        symbolic_complete(false),
        num_iters(-1),
        end_rel_res(-1) {}

  KOKKOS_INLINE_FUNCTION
  ~PAR_ILUTHandle() {}

  KOKKOS_INLINE_FUNCTION
  size_type get_nrows() const { return nrows; }

  KOKKOS_INLINE_FUNCTION
  void set_nrows(const size_type nrows_) { this->nrows = nrows_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nnzL() const { return nnzL; }

  KOKKOS_INLINE_FUNCTION
  void set_nnzL(const size_type nnzL_) { this->nnzL = nnzL_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nnzU() const { return nnzU; }

  KOKKOS_INLINE_FUNCTION
  void set_nnzU(const size_type nnzU_) { this->nnzU = nnzU_; }

  bool is_symbolic_complete() const { return symbolic_complete; }

  void set_symbolic_complete() { this->symbolic_complete = true; }
  void reset_symbolic_complete() { this->symbolic_complete = false; }

  void set_team_size(const int ts) { this->team_size = ts; }
  int get_team_size() const { return this->team_size; }

  void set_vector_size(const int vs) { this->vector_size = vs; }
  int get_vector_size() const { return this->vector_size; }

  void set_max_iter(const size_type max_iter_) { this->max_iter = max_iter_; }
  int get_max_iter() const { return this->max_iter; }

  void set_residual_norm_delta_stop(const float_t residual_norm_delta_stop_) {
    this->residual_norm_delta_stop = residual_norm_delta_stop_;
  }
  float_t get_residual_norm_delta_stop() const { return this->residual_norm_delta_stop; }

  void set_fill_in_limit(const float_t fill_in_limit_) { this->fill_in_limit = fill_in_limit_; }
  float_t get_fill_in_limit() const { return this->fill_in_limit; }

  bool get_verbose() const { return verbose; }

  void set_verbose(const bool verbose_) { this->verbose = verbose_; }

  bool get_async_update() const { return async_update; }

  void set_async_update(const bool async_update_) { this->async_update = async_update_; }

  TeamPolicy get_default_team_policy() const {
    if (team_size == -1) {
      return TeamPolicy(nrows, Kokkos::AUTO);
    } else {
      return TeamPolicy(nrows, team_size);
    }
  }

  int get_num_iters() const { return num_iters; }

  nnz_scalar_t get_end_rel_res() const { return end_rel_res; }

  void set_stats(int num_iters_, nnz_scalar_t end_rel_res_) {
    num_iters   = num_iters_;
    end_rel_res = end_rel_res_;
  }
};

}  // namespace Experimental
}  // namespace KokkosSparse

#endif
