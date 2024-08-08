/*
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
*/

#include <Kokkos_Core.hpp>
#include <KokkosSparse_Preconditioner.hpp>
#include <iostream>
#include <string>

#ifndef _GMRESHANDLE_HPP
#define _GMRESHANDLE_HPP

namespace KokkosSparse {
namespace Experimental {

/**
 * The handle class for GMRES. Used to store some input parameters and
 * results.
 *
 * For more info, see KokkosSparse_gmres.hpp doxygen
 */
template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class GMRESHandle {
 public:
  using HandleExecSpace             = ExecutionSpace;
  using HandleTempMemorySpace       = TemporaryMemorySpace;
  using HandlePersistentMemorySpace = PersistentMemorySpace;

  using execution_space = ExecutionSpace;
  using memory_space    = HandlePersistentMemorySpace;
  using device_t        = Kokkos::Device<execution_space, memory_space>;

  using size_type       = typename std::remove_const<size_type_>::type;
  using const_size_type = const size_type;

  using nnz_lno_t       = typename std::remove_const<lno_t_>::type;
  using const_nnz_lno_t = const nnz_lno_t;

  using nnz_scalar_t       = typename std::remove_const<scalar_t_>::type;
  using const_nnz_scalar_t = const nnz_scalar_t;

  using float_t = typename Kokkos::ArithTraits<nnz_scalar_t>::mag_type;

  using nnz_row_view_t = typename Kokkos::View<size_type *, device_t>;

  using nnz_lno_view_t = typename Kokkos::View<nnz_lno_t *, device_t>;

  using nnz_value_view_t = typename Kokkos::View<nnz_scalar_t *, device_t>;

  using nnz_value_view2d_t = typename Kokkos::View<nnz_scalar_t **, Kokkos::LayoutLeft, device_t>;

  using signed_integral_t = typename std::make_signed<typename nnz_row_view_t::non_const_value_type>::type;

  using signed_nnz_lno_view_t =
      Kokkos::View<signed_integral_t *, typename nnz_row_view_t::array_layout, typename nnz_row_view_t::device_type,
                   typename nnz_row_view_t::memory_traits>;

  /**
   * The orthogonalization type
   */
  enum Ortho {
    CGS2,  // Two iterations of Classical Gram-Schmidt
    MGS
  };  // One iteration of Modified Gram-Schmidt

  /**
   * The result of the run
   */
  enum Flag {
    Conv,    // Converged
    NoConv,  // Did not converge
    LOA,     // Solver had loss of accuracy
    NotRun
  };  // GMRES was never run

 private:
  // Inputs

  size_type m;            /// Max subspace size before restarting
  float_t tol;            /// Relative residual convergence tolerance
  size_type max_restart;  /// Maximum number of times to restart the solver
  Ortho ortho;            /// The orthogonalization type
  bool verbose;           /// Print extra info to stdout

  // Outputs
  int num_iters;        /// Number of iterations the sovler took
  float_t end_rel_res;  /// Residual from solver
  Flag conv_flag_val;   /// Denotes end result of the run

 public:
  // Use set methods to control ortho, and verbose
  GMRESHandle(const size_type m_ = 50, const float_t tol_ = 1e-8, const size_type max_restart_ = 50)
      : m(m_),
        tol(tol_),
        max_restart(max_restart_),
        ortho(CGS2),
        verbose(false),
        num_iters(-1),
        end_rel_res(-1),
        conv_flag_val(NotRun) {
    if (m <= 0) {
      throw std::invalid_argument("gmres: Please choose restart size m greater than zero.");
    }
  }

  void reset_handle(const size_type m_ = 50, const float_t tol_ = 1e-8, const size_type max_restart_ = 50) {
    set_m(m_);
    set_tol(tol_);
    set_max_restart(max_restart_);
    set_ortho(CGS2);
    set_verbose(false);
    num_iters     = -1;
    end_rel_res   = -1;
    conv_flag_val = NotRun;
  }

  KOKKOS_INLINE_FUNCTION
  ~GMRESHandle() {}

  KOKKOS_INLINE_FUNCTION
  size_type get_m() const { return m; }

  KOKKOS_INLINE_FUNCTION
  void set_m(const size_type m_) { this->m = m_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_max_restart() const { return max_restart; }

  KOKKOS_INLINE_FUNCTION
  void set_max_restart(const size_type max_restart_) { this->max_restart = max_restart_; }

  KOKKOS_INLINE_FUNCTION
  float_t get_tol() const { return tol; }

  KOKKOS_INLINE_FUNCTION
  void set_tol(const float_t tol_) { this->tol = tol_; }

  KOKKOS_INLINE_FUNCTION
  Ortho get_ortho() const { return ortho; }

  KOKKOS_INLINE_FUNCTION
  void set_ortho(const Ortho ortho_) { this->ortho = ortho_; }

  KOKKOS_INLINE_FUNCTION
  bool get_verbose() const { return verbose; }

  KOKKOS_INLINE_FUNCTION
  void set_verbose(const bool verbose_) { this->verbose = verbose_; }

  int get_num_iters() const {
    assert(get_conv_flag_val() != NotRun);
    return num_iters;
  }
  float_t get_end_rel_res() const {
    assert(get_conv_flag_val() != NotRun);
    return end_rel_res;
  }
  Flag get_conv_flag_val() const { return conv_flag_val; }

  void set_stats(int num_iters_, float_t end_rel_res_, Flag conv_flag_val_) {
    assert(conv_flag_val_ != NotRun);
    num_iters     = num_iters_;
    end_rel_res   = end_rel_res_;
    conv_flag_val = conv_flag_val_;
  }
};

}  // namespace Experimental
}  // namespace KokkosSparse

#endif
