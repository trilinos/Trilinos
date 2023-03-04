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

template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace,
          class TemporaryMemorySpace, class PersistentMemorySpace>
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

  using nnz_row_view_t =
      typename Kokkos::View<size_type *, HandlePersistentMemorySpace>;

  using nnz_lno_view_t =
      typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace>;

  using nnz_value_view_t =
      typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace>;

  using signed_integral_t = typename std::make_signed<
      typename nnz_row_view_t::non_const_value_type>::type;

  using signed_nnz_lno_view_t =
      Kokkos::View<signed_integral_t *, typename nnz_row_view_t::array_layout,
                   typename nnz_row_view_t::device_type,
                   typename nnz_row_view_t::memory_traits>;

 private:
  size_type nrows;
  size_type nnzL;
  size_type nnzU;
  size_type max_iter;
  nnz_scalar_t residual_norm_delta_stop;

  bool symbolic_complete;

  int team_size;
  int vector_size;

  float_t fill_in_limit;

 public:
  PAR_ILUTHandle(const size_type nrows_, const size_type nnzL_ = 0,
                 const size_type nnzU_ = 0, const size_type max_iter_ = 1)
      : nrows(nrows_),
        nnzL(nnzL_),
        nnzU(nnzU_),
        max_iter(max_iter_),
        residual_norm_delta_stop(0.),
        symbolic_complete(false),
        team_size(-1),
        vector_size(-1),
        fill_in_limit(0.75) {}

  void reset_handle(const size_type nrows_, const size_type nnzL_,
                    const size_type nnzU_) {
    set_nrows(nrows_);
    set_nnzL(nnzL_);
    set_nnzU(nnzU_);
    set_residual_norm_delta_stop(0.);
    reset_symbolic_complete();
    set_fill_in_limit(0.75);
  }

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

  void set_residual_norm_delta_stop(
      const nnz_scalar_t residual_norm_delta_stop_) {
    this->residual_norm_delta_stop = residual_norm_delta_stop_;
  }
  nnz_scalar_t get_residual_norm_delta_stop() const {
    return this->residual_norm_delta_stop;
  }

  void set_fill_in_limit(const float_t fill_in_limit_) {
    this->fill_in_limit = fill_in_limit_;
  }
  float_t get_fill_in_limit() const { return this->fill_in_limit; }

  TeamPolicy get_default_team_policy() const {
    if (team_size == -1) {
      return TeamPolicy(nrows, Kokkos::AUTO);
    } else {
      return TeamPolicy(nrows, team_size);
    }
  }
};

}  // namespace Experimental
}  // namespace KokkosSparse

#endif
