// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_NUMERIC_TOOLS_SERIAL_HPP__
#define __TACHO_NUMERIC_TOOLS_SERIAL_HPP__

/// \file Tacho_NumericTools_Serial.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_NumericTools_Base.hpp"

#include "Tacho_DenseFlopCount.hpp"
#include "Tacho_SupernodeInfo.hpp"

#include "Tacho_CholSupernodes.hpp"
#include "Tacho_CholSupernodes_Serial.hpp"
#include "Tacho_CholSupernodes_SerialPanel.hpp"

#include "Tacho_LDL_Supernodes.hpp"
#include "Tacho_LDL_Supernodes_Serial.hpp"

#include "Tacho_LU_Supernodes.hpp"
#include "Tacho_LU_Supernodes_Serial.hpp"

namespace Tacho {

template <typename ValueType, typename DeviceType>
class NumericToolsSerial : public NumericToolsBase<ValueType, DeviceType> {
public:
  ///
  /// types
  ///
  using base_type = NumericToolsBase<ValueType, DeviceType>;
  using typename base_type::device_type;
  using typename base_type::exec_memory_space;
  using typename base_type::exec_space;
  using typename base_type::host_memory_space;
  using typename base_type::ordinal_type_array;
  using typename base_type::ordinal_type_array_host;
  using typename base_type::size_type_array;
  using typename base_type::value_type;
  using typename base_type::value_type_array;
  using typename base_type::value_type_matrix;

  using base_type::base_type;

private:
  using base_type::_aj;
  using base_type::_ap;
  using base_type::_ax;
  using base_type::_diag;
  using base_type::_info;
  using base_type::_m;
  using base_type::_nsupernodes;
  using base_type::_peri;
  using base_type::_perm;
  using base_type::_piv;
  using base_type::_stree_roots;
  using base_type::_supernodes;
  using base_type::_superpanel_buf;

  using base_type::reset_stat;
  using base_type::stat;
  using base_type::track_alloc;
  using base_type::track_free;
  // using base_type::print_stat_factor;
  using base_type::print_stat_memory;
  using base_type::print_stat_solve;

public:
  inline void print_stat_factor() override {
    base_type::print_stat_factor();
    double flop = 0;
    auto h_supernodes = Kokkos::create_mirror_view_and_copy(host_memory_space(), _supernodes);
    switch (this->getSolutionMethod()) {
    case 1: {
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid) {
        auto &s = h_supernodes(sid);
        const ordinal_type m = s.m, n = s.n - s.m;
        flop += DenseFlopCount<value_type>::Chol(m);
        flop += DenseFlopCount<value_type>::Trsm(true, m, n);
        flop += DenseFlopCount<value_type>::Syrk(m, n);
      }
      break;
    }
    case 2: {
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid) {
        auto &s = h_supernodes(sid);
        const ordinal_type m = s.m, n = s.n - s.m;
        flop += DenseFlopCount<value_type>::LDL(m);
        flop += DenseFlopCount<value_type>::Trsm(true, m, n);
        flop += DenseFlopCount<value_type>::Syrk(m, n);
      }
      break;
    }
    case 3: {
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid) {
        auto &s = h_supernodes(sid);
        const ordinal_type m = s.m, n = s.n - s.m;
        flop += DenseFlopCount<value_type>::LU(m, m);
        flop += 2 * DenseFlopCount<value_type>::Trsm(true, m, n);
        flop += DenseFlopCount<value_type>::Gemm(n, n, m);
      }
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(false, std::logic_error, "The solution method is not supported");
    }
    }
    const double kilo(1024);
    printf("  FLOPs\n");
    printf("             gflop   for numeric factorization:               %10.3f GFLOP\n", flop / kilo / kilo / kilo);
    printf("             gflop/s for numeric factorization:               %10.3f GFLOP/s\n",
           flop / stat.t_factor / kilo / kilo / kilo);
    printf("\n");
  }

  ///
  /// Choleksy
  ///
  inline void factorizeCholesky(const value_type_array &ax, const ordinal_type verbose) {
    Kokkos::Timer timer;
    {
      timer.reset();
      {
        /// matrix values
        _ax = ax;

        /// copy the input matrix into super panels
        const bool copy_to_l_buf(false);
        _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, _perm, _peri);
      }
      stat.t_copy = timer.seconds();
    }

    timer.reset();
    {
      /// valgrind reports the following buf array as uninitialized even if it is initialized
      /// while the task is executed. to remove the valgrind error, we initialize the array with zero.
      /// value_type_array buf(do_not_initialize_tag("buf"), _info.max_schur_size*(_info.max_schur_size + 1));
      value_type_array buf("buf", _info.max_schur_size * (_info.max_schur_size + 1));
      const size_t bufsize = buf.span() * sizeof(value_type);
      track_alloc(bufsize);

      /// recursive tree traversal
      const ordinal_type member = 0, nroots = _stree_roots.extent(0);
      for (ordinal_type i = 0; i < nroots; ++i)
        CholSupernodes<Algo::Workflow::Serial>::factorize_recursive_serial(member, _info, _stree_roots(i), true,
                                                                           buf.data(), bufsize);

      track_free(bufsize);
    }
    stat.t_factor = timer.seconds();

    if (verbose) {
      printf("Summary: NumericTools, Cholesky (SerialFactorization)\n");
      printf("=====================================================\n");

      print_stat_factor();
    }
  }

  inline void factorizeCholesky(const value_type_array &ax, const ordinal_type panelsize, const ordinal_type verbose) {
    Kokkos::Timer timer;
    {
      timer.reset();
      {
        /// matrix values
        _ax = ax;

        /// copy the input matrix into super panels
        const bool copy_to_l_buf(false);
        _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, _perm, _peri);
      }
      stat.t_copy = timer.seconds();
    }

    const ordinal_type nb = panelsize > 0 ? panelsize : _info.max_schur_size;
    timer.reset();
    {
      /// valgrind reports the following buf array as uninitialized even if it is initialized
      /// while the task is executed. to remove the valgrind error, we initialize the array with zero.
      /// value_type_array buf(do_not_initialize_tag("buf"), _info.max_schur_size*(nb + 1));
      value_type_array buf("buf", _info.max_schur_size * (nb + 1));
      const size_t bufsize = buf.span() * sizeof(value_type);
      track_alloc(bufsize);

      /// recursive tree traversal
      const ordinal_type member = 0, nroots = _stree_roots.extent(0);
      for (ordinal_type i = 0; i < nroots; ++i)
        CholSupernodes<Algo::Workflow::SerialPanel>::factorize_recursive_serial(member, _info, _stree_roots(i), true,
                                                                                buf.data(), bufsize, nb);

      track_free(bufsize);
    }
    stat.t_factor = timer.seconds();

    if (verbose) {
      printf("Summary: NumericTools, Cholesky (SerialPanelFactorization: %3d)\n", nb);
      printf("===============================================================\n");

      print_stat_factor();
    }
  }

  inline void solveCholesky(const value_type_matrix &x, // solution
                            const value_type_matrix &b, // right hand side
                            const value_type_matrix &t, // temporary workspace (store permuted vectors)
                            const ordinal_type verbose) {
    Kokkos::Timer timer;

    _info.x = t;

    // copy b -> t
    timer.reset();
    const auto exec_instance = exec_space();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, b, _perm, t);
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      value_type_array buf(do_not_initialize_tag("buf"), _info.max_schur_size * x.extent(1));
      const size_t bufsize = buf.span() * sizeof(value_type);
      track_alloc(bufsize);

      /// recursive tree traversal
      const ordinal_type member = 0, nroots = _stree_roots.extent(0);
      for (ordinal_type i = 0; i < nroots; ++i)
        CholSupernodes<Algo::Workflow::Serial>::solve_lower_recursive_serial(member, _info, _stree_roots(i), true,
                                                                             buf.data(), bufsize);
      for (ordinal_type i = 0; i < nroots; ++i)
        CholSupernodes<Algo::Workflow::Serial>::solve_upper_recursive_serial(member, _info, _stree_roots(i), true,
                                                                             buf.data(), bufsize);

      track_free(bufsize);
    }
    stat.t_solve = timer.seconds();

    // copy t -> x
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, t, _peri, x);
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: NumericTools, Cholesky (SerialSolve: %3d)\n", ordinal_type(x.extent(1)));
      printf("==================================================\n");

      print_stat_solve();
    }
  }

  ///
  /// LDL
  ///
  inline void factorizeLDL(const value_type_array &ax, const ordinal_type verbose) {
    {
      const bool test = !std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
      TACHO_TEST_FOR_EXCEPTION(test, std::logic_error, "Serial interface works on host device only");
    }

    Kokkos::Timer timer;
    {
      timer.reset();
      {
        /// matrix values
        _ax = ax;

        /// copy the input matrix into super panels
        const bool copy_to_l_buf(false);
        _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, _perm, _peri);
      }
      stat.t_copy = timer.seconds();
    }

    timer.reset();
    {
      /// valgrind reports the following buf array as uninitialized even if it is initialized
      /// while the task is executed. to remove the valgrind error, we initialize the array with zero.
      /// value_type_array buf(do_not_initialize_tag("buf"), _info.max_schur_size*(_info.max_schur_size + 1));
      value_type_array buf("buf",
                           _info.max_schur_size * (_info.max_schur_size + 1) + // ABR
                               _info.max_supernode_size *
                                   std::max(32, _info.max_schur_size)); // ATR copy and workspace for LDL
      const size_t bufsize = buf.span() * sizeof(value_type);
      track_alloc(bufsize);

      /// recursive tree traversal
      const ordinal_type member = 0, nroots = _stree_roots.extent(0);
      for (ordinal_type i = 0; i < nroots; ++i)
        LDL_Supernodes<Algo::Workflow::Serial>::factorize_recursive_serial(
            member, _info, _stree_roots(i), true, _piv.data(), _diag.data(), buf.data(), bufsize);

      track_free(bufsize);
    }
    stat.t_factor = timer.seconds();

    if (verbose) {
      printf("Summary: NumericTools, LDL (SerialFactorization)\n");
      printf("================================================\n");

      print_stat_factor();
    }
  }

  inline void solveLDL(const value_type_matrix &x, // solution
                       const value_type_matrix &b, // right hand side
                       const value_type_matrix &t, // temporary workspace (store permuted vectors)
                       const ordinal_type verbose) {
    {
      const bool test = !std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
      TACHO_TEST_FOR_EXCEPTION(test, std::logic_error, "Serial interface works on host device only");
      TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) || x.extent(1) != b.extent(1) || x.extent(0) != t.extent(0) ||
                                   x.extent(1) != t.extent(1),
                               std::logic_error, "supernode data structure is not allocated");
      TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() || x.data() == t.data() || t.data() == b.data(), std::logic_error,
                               "x, b and t have the same data pointer");
    }

    Kokkos::Timer timer;

    _info.x = t;

    // copy b -> t
    timer.reset();
    const auto exec_instance = exec_space();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, b, _perm, t);
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      value_type_array buf(do_not_initialize_tag("buf"), _info.max_schur_size * x.extent(1));
      const size_t bufsize = buf.span() * sizeof(value_type);
      track_alloc(bufsize);

      /// recursive tree traversal
      const ordinal_type member = 0, nroots = _stree_roots.extent(0);
      for (ordinal_type i = 0; i < nroots; ++i)
        LDL_Supernodes<Algo::Workflow::Serial>::solve_lower_recursive_serial(member, _info, _stree_roots(i), true,
                                                                             _piv.data(), buf.data(), bufsize);
      for (ordinal_type i = 0; i < nroots; ++i)
        LDL_Supernodes<Algo::Workflow::Serial>::solve_upper_recursive_serial(
            member, _info, _stree_roots(i), true, _piv.data(), _diag.data(), buf.data(), bufsize);

      track_free(bufsize);
    }
    stat.t_solve = timer.seconds();

    // copy t -> x
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, t, _peri, x);
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: NumericTools, LDL (SerialSolve: %3d)\n", ordinal_type(x.extent(1)));
      printf("=============================================\n");

      print_stat_solve();
    }
  }

  ///
  /// LU
  ///
  inline void factorizeLU(const value_type_array &ax, const ordinal_type verbose) {
    {
      const bool test = !std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
      TACHO_TEST_FOR_EXCEPTION(test, std::logic_error, "Serial interface works on host device only");
    }

    Kokkos::Timer timer;
    {
      timer.reset();
      {
        /// matrix values
        _ax = ax;

        /// copy the input matrix into super panels
        const bool copy_to_l_buf(true);
        _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, _perm, _peri);
      }
      stat.t_copy = timer.seconds();
    }

    timer.reset();
    {
      /// valgrind reports the following buf array as uninitialized even if it is initialized
      /// while the task is executed. to remove the valgrind error, we initialize the array with zero.
      /// value_type_array buf(do_not_initialize_tag("buf"), _info.max_schur_size*(_info.max_schur_size + 1));
      value_type_array buf("buf",
                           _info.max_schur_size * (_info.max_schur_size + 1)); // ABR
      const size_t bufsize = buf.span() * sizeof(value_type);
      track_alloc(bufsize);

      /// recursive tree traversal
      const ordinal_type member = 0, nroots = _stree_roots.extent(0);
      for (ordinal_type i = 0; i < nroots; ++i)
        LU_Supernodes<Algo::Workflow::Serial>::factorize_recursive_serial(member, _info, _stree_roots(i), true,
                                                                          _piv.data(), buf.data(), bufsize);

      track_free(bufsize);
    }
    stat.t_factor = timer.seconds();

    if (verbose) {
      printf("Summary: NumericTools, LU (SerialFactorization)\n");
      printf("===============================================\n");

      print_stat_factor();
    }
  }

  inline void solveLU(const value_type_matrix &x, // solution
                      const value_type_matrix &b, // right hand side
                      const value_type_matrix &t, // temporary workspace (store permuted vectors)
                      const ordinal_type verbose) {
    {
      const bool test = !std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
      TACHO_TEST_FOR_EXCEPTION(test, std::logic_error, "Serial interface works on host device only");
      TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) || x.extent(1) != b.extent(1) || x.extent(0) != t.extent(0) ||
                                   x.extent(1) != t.extent(1),
                               std::logic_error, "supernode data structure is not allocated");
      TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() || x.data() == t.data() || t.data() == b.data(), std::logic_error,
                               "x, b and t have the same data pointer");
    }

    Kokkos::Timer timer;

    _info.x = t;

    // copy b -> t
    timer.reset();
    const auto exec_instance = exec_space();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, b, _perm, t);
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      value_type_array buf(do_not_initialize_tag("buf"), _info.max_schur_size * x.extent(1));
      const size_t bufsize = buf.span() * sizeof(value_type);
      track_alloc(bufsize);

      /// recursive tree traversal
      const ordinal_type member = 0, nroots = _stree_roots.extent(0);
      for (ordinal_type i = 0; i < nroots; ++i)
        LU_Supernodes<Algo::Workflow::Serial>::solve_lower_recursive_serial(member, _info, _stree_roots(i), true,
                                                                            _piv.data(), buf.data(), bufsize);
      for (ordinal_type i = 0; i < nroots; ++i)
        LU_Supernodes<Algo::Workflow::Serial>::solve_upper_recursive_serial(member, _info, _stree_roots(i), true,
                                                                            _piv.data(), buf.data(), bufsize);

      track_free(bufsize);
    }
    stat.t_solve = timer.seconds();

    // copy t -> x
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, t, _peri, x);
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: NumericTools, LU (SerialSolve: %3d)\n", ordinal_type(x.extent(1)));
      printf("============================================\n");

      print_stat_solve();
    }
  }

  ///
  /// main interface
  ///
  inline void factorize(const value_type_array &ax, const ordinal_type verbose = 0) override {
    {
      const bool test = !std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
      TACHO_TEST_FOR_EXCEPTION(test, std::logic_error, "Serial interface works on host device only");
    }
    /// reset the supernode buffer for potential reuse cases
    Kokkos::deep_copy(_superpanel_buf, value_type(0));
    switch (this->getSolutionMethod()) {
    case 1: { /// Cholesky
      // if (_nb > 0) {
      //   factorizeCholesky(ax, _nb, verbose);
      // } else {
      //   factorizeCholesky(ax, verbose);
      // }
      factorizeCholesky(ax, verbose);
      break;
    }
    case 2: { /// LDL
      {
        const ordinal_type rlen = 4 * _m, plen = _piv.span();
        if (plen < rlen) {
          track_free(this->_piv.span() * sizeof(ordinal_type));
          this->_piv = ordinal_type_array("piv", rlen);
          track_alloc(this->_piv.span() * sizeof(ordinal_type));
        }
      }
      {
        const ordinal_type rlen = 2 * _m, dlen = _diag.span();
        if (dlen < rlen) {
          track_free(this->_diag.span() * sizeof(value_type));
          this->_diag = value_type_array("diag", rlen);
          track_alloc(this->_diag.span() * sizeof(value_type));
        }
      }
      factorizeLDL(ax, verbose);
      break;
    }
    case 3: { /// LU
      {
        const ordinal_type rlen = 4 * _m, plen = _piv.span();
        if (plen < rlen) {
          track_free(this->_piv.span() * sizeof(ordinal_type));
          this->_piv = ordinal_type_array("piv", rlen);
          track_alloc(this->_piv.span() * sizeof(ordinal_type));
        }
      }
      factorizeLU(ax, verbose);
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(false, std::logic_error, "The solution method is not supported");
      break;
    }
    }
  }

  inline void solve(const value_type_matrix &x, // solution
                    const value_type_matrix &b, // right hand side
                    const value_type_matrix &t, // temporary workspace (store permuted vectors)
                    const ordinal_type verbose = 0) override {
    {
      const bool test = !std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
      TACHO_TEST_FOR_EXCEPTION(test, std::logic_error, "Serial interface works on host device only");
      TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) || x.extent(1) != b.extent(1) || x.extent(0) != t.extent(0) ||
                                   x.extent(1) != t.extent(1),
                               std::logic_error, "Input x, b and t dimensions are not compatible");
      TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() || x.data() == t.data() || t.data() == b.data(), std::logic_error,
                               "Input x, b and t have the same data pointer");
    }

    switch (this->getSolutionMethod()) {
    case 1: {
      solveCholesky(x, b, t, verbose);
      break;
    }
    case 2: {
      solveLDL(x, b, t, verbose);
      break;
    }
    case 3: {
      solveLU(x, b, t, verbose);
      break;
    }
    default: {
    }
    }
  }
};

} // namespace Tacho
#endif
