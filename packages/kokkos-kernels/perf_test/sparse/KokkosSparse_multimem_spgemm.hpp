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

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_run_spgemm.hpp"
#include "KokkosSparse_IOUtils.hpp"

namespace KokkosKernels {

namespace Experiment {

template <typename size_type, typename lno_t, typename scalar_t,
          typename exec_space, typename hbm_mem_space, typename sbm_mem_space>
void run_multi_mem_spgemm(Parameters params) {
  typedef exec_space myExecSpace;
  typedef Kokkos::Device<exec_space, hbm_mem_space> myFastDevice;
  typedef Kokkos::Device<exec_space, sbm_mem_space> mySlowExecSpace;

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, myFastDevice, void,
                                           size_type>
      fast_crstmat_t;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, mySlowExecSpace,
                                           void, size_type>
      slow_crstmat_t;

  char *a_mat_file = params.a_mtx_bin_file;
  char *b_mat_file = params.b_mtx_bin_file;
  char *c_mat_file = params.c_mtx_bin_file;

  slow_crstmat_t a_slow_crsmat, b_slow_crsmat, c_slow_crsmat;
  fast_crstmat_t a_fast_crsmat, b_fast_crsmat, c_fast_crsmat;

  // read a and b matrices and store them on slow or fast memory.

  if (params.a_mem_space == 1) {
    a_fast_crsmat =
        KokkosSparse::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(a_mat_file);
  } else {
    a_slow_crsmat =
        KokkosSparse::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(a_mat_file);
  }

  if ((b_mat_file == NULL || strcmp(b_mat_file, a_mat_file) == 0) &&
      params.b_mem_space == params.a_mem_space) {
    std::cout << "Using A matrix for B as well" << std::endl;
    b_fast_crsmat = a_fast_crsmat;
    b_slow_crsmat = a_slow_crsmat;
  } else if (params.b_mem_space == 1) {
    if (b_mat_file == NULL) b_mat_file = a_mat_file;
    b_fast_crsmat =
        KokkosSparse::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(b_mat_file);
  } else {
    if (b_mat_file == NULL) b_mat_file = a_mat_file;
    b_slow_crsmat =
        KokkosSparse::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(b_mat_file);
  }

  if (params.a_mem_space == 1) {
    if (params.b_mem_space == 1) {
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, fast_crstmat_t, fast_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_fast_crsmat, b_fast_crsmat,
                                            params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, fast_crstmat_t, fast_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_fast_crsmat, b_fast_crsmat,
                                            params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, fast_crstmat_t, slow_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_fast_crsmat, b_fast_crsmat,
                                            params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, fast_crstmat_t, slow_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_fast_crsmat, b_fast_crsmat,
                                            params);
        }
      }
    } else {
      // B is in slow memory
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, slow_crstmat_t, fast_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_fast_crsmat, b_slow_crsmat,
                                            params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, slow_crstmat_t, fast_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_fast_crsmat, b_slow_crsmat,
                                            params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, slow_crstmat_t, slow_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_fast_crsmat, b_slow_crsmat,
                                            params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, fast_crstmat_t, slow_crstmat_t, slow_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_fast_crsmat, b_slow_crsmat,
                                            params);
        }
      }
    }
  } else {
    // A is in slow memory
    if (params.b_mem_space == 1) {
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, fast_crstmat_t, fast_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_slow_crsmat, b_fast_crsmat,
                                            params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, fast_crstmat_t, fast_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_slow_crsmat, b_fast_crsmat,
                                            params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, fast_crstmat_t, slow_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_slow_crsmat, b_fast_crsmat,
                                            params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, fast_crstmat_t, slow_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_slow_crsmat, b_fast_crsmat,
                                            params);
        }
      }
    } else {
      // B is in slow memory
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, slow_crstmat_t, fast_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_slow_crsmat, b_slow_crsmat,
                                            params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, slow_crstmat_t, fast_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_slow_crsmat, b_slow_crsmat,
                                            params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, slow_crstmat_t, slow_crstmat_t,
              hbm_mem_space, hbm_mem_space>(a_slow_crsmat, b_slow_crsmat,
                                            params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<
              myExecSpace, slow_crstmat_t, slow_crstmat_t, slow_crstmat_t,
              sbm_mem_space, sbm_mem_space>(a_slow_crsmat, b_slow_crsmat,
                                            params);
        }
      }
    }
  }

  if (c_mat_file != NULL) {
    if (params.c_mem_space == 1) {
      KokkosSparse::sort_crs_matrix(c_fast_crsmat);

      KokkosSparse::Impl::write_graph_bin(
          (lno_t)(c_fast_crsmat.numRows()),
          (size_type)(c_fast_crsmat.graph.entries.extent(0)),
          c_fast_crsmat.graph.row_map.data(),
          c_fast_crsmat.graph.entries.data(), c_fast_crsmat.values.data(),
          c_mat_file);
    } else {
      KokkosSparse::sort_crs_matrix(c_slow_crsmat);

      KokkosSparse::Impl::write_graph_bin(
          (lno_t)c_slow_crsmat.numRows(),
          (size_type)c_slow_crsmat.graph.entries.extent(0),
          c_slow_crsmat.graph.row_map.data(),
          c_slow_crsmat.graph.entries.data(), c_slow_crsmat.values.data(),
          c_mat_file);
    }
  }
}

}  // namespace Experiment
}  // namespace KokkosKernels
