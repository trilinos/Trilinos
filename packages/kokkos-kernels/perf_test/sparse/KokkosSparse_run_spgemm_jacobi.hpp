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
#include "KokkosKernels_TestParameters.hpp"
#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "KokkosSparse_IOUtils.hpp"

#define TRANSPOSEFIRST false
#define TRANSPOSESECOND false

namespace KokkosKernels {

namespace Experiment {

template <typename crsMat_t, typename device>
bool is_same_matrix(crsMat_t output_mat1, crsMat_t output_mat2) {
  using graph_t        = typename crsMat_t::StaticCrsGraphType;
  using lno_view_t     = typename graph_t::row_map_type::non_const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t  = typename crsMat_t::values_type::non_const_type;

  size_t nrows1    = output_mat1.graph.row_map.extent(0);
  size_t nentries1 = output_mat1.graph.entries.extent(0);
  size_t nvals1    = output_mat1.values.extent(0);

  size_t nrows2    = output_mat2.graph.row_map.extent(0);
  size_t nentries2 = output_mat2.graph.entries.extent(0);
  size_t nvals2    = output_mat2.values.extent(0);

  KokkosSparse::sort_crs_matrix(output_mat1);

  if (nrows1 != nrows2) {
    std::cerr << "row count is different" << std::endl;
    return false;
  }
  if (nentries1 != nentries2) {
    std::cerr << "nentries2 is different" << std::endl;
    return false;
  }
  if (nvals1 != nvals2) {
    std::cerr << "nvals1 is different" << std::endl;
    return false;
  }

  KokkosSparse::sort_crs_matrix(output_mat2);

  bool is_identical = true;
  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<typename graph_t::row_map_type, typename graph_t::row_map_type,
                                                typename lno_view_t::value_type, typename device::execution_space>(
          output_mat1.graph.row_map, output_mat2.graph.row_map, 0);
  if (!is_identical) {
    std::cerr << "rowmaps are different" << std::endl;
    return false;
  }

  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<lno_nnz_view_t, lno_nnz_view_t, typename lno_nnz_view_t::value_type,
                                                typename device::execution_space>(output_mat1.graph.entries,
                                                                                  output_mat2.graph.entries, 0);

  if (!is_identical) {
    std::cerr << "entries are different" << std::endl;
    return false;
  }

  is_identical =
      KokkosKernels::Impl::kk_is_identical_view<scalar_view_t, scalar_view_t, typename scalar_view_t::value_type,
                                                typename device::execution_space>(output_mat1.values,
                                                                                  output_mat2.values, 0.00001);
  if (!is_identical) {
    std::cerr << "values are different" << std::endl;
  }
  return true;
}

template <typename ExecSpace, typename crsMat_t, typename crsMat_t2, typename crsMat_t3, typename TempMemSpace,
          typename PersistentMemSpace>
crsMat_t3 run_experiment(crsMat_t crsMat, crsMat_t2 crsMat2, Parameters params) {
  using namespace KokkosSparse;
  using namespace KokkosSparse::Experimental;
  using device_t       = Kokkos::Device<ExecSpace, PersistentMemSpace>;
  using scalar_view_t  = typename crsMat_t3::values_type::non_const_type;
  using lno_view_t     = typename crsMat_t3::row_map_type::non_const_type;
  using lno_nnz_view_t = typename crsMat_t3::index_type::non_const_type;
  using lno_t          = typename lno_nnz_view_t::value_type;
  using size_type      = typename lno_view_t::value_type;
  using scalar_t       = typename scalar_view_t::value_type;
  using KernelHandle   = KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, ExecSpace,
                                                                        TempMemSpace, PersistentMemSpace>;

  int algorithm                 = params.algorithm;
  int repeat                    = params.repeat;
  int chunk_size                = params.chunk_size;
  int shmemsize                 = params.shmemsize;
  int team_size                 = params.team_size;
  int use_dynamic_scheduling    = params.use_dynamic_scheduling;
  int verbose                   = params.verbose;
  int calculate_read_write_cost = params.calculate_read_write_cost;
  int vector_size               = params.vector_size;
  int check_output              = params.check_output;

  lno_view_t row_mapC;
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  KernelHandle kh;
  kh.set_team_work_size(chunk_size);
  kh.set_shmem_size(shmemsize);
  kh.set_suggested_team_size(team_size);
  kh.set_suggested_vector_size(vector_size);

  if (use_dynamic_scheduling) {
    kh.set_dynamic_scheduling(true);
  }
  if (verbose) {
    kh.set_verbose(true);
  }

  const lno_t m = crsMat.numRows();
  const lno_t n = crsMat2.numRows();
  const lno_t k = crsMat2.numCols();

  if (verbose) std::cout << "m:" << m << " n:" << n << " k:" << k << std::endl;
  if (m != n) {
    std::cerr << "left.numCols():" << n << " left.numRows():" << m << std::endl;
    exit(1);
  }
  if (n < crsMat.numCols()) {
    std::cerr << "left.numCols():" << crsMat.numCols() << " right.numRows():" << crsMat2.numRows() << std::endl;
    exit(1);
  }

  typedef typename Kokkos::View<scalar_t **,
                                typename KokkosKernels::Impl::GetUnifiedLayout<scalar_view_t>::array_layout, device_t>
      view_t;

  view_t dinv("Dinv", m, 1);
  Kokkos::deep_copy(dinv, 2.0);
  scalar_t omega = 3;

  // The reference product (for verifying correctness)
  // Don't allocate them if they won't be used, but they must be declared here.
  lno_view_t row_mapC_ref;
  lno_nnz_view_t entriesC_ref;
  scalar_view_t valuesC_ref;
  crsMat_t3 Ccrsmat_ref;

  if (check_output) {
    if (verbose) std::cout << "Running a reference algorithm" << std::endl;

    row_mapC_ref = lno_view_t("non_const_lnow_row", m + 1);

    KernelHandle sequential_kh;
    sequential_kh.set_team_work_size(chunk_size);
    sequential_kh.set_shmem_size(shmemsize);
    sequential_kh.set_suggested_team_size(team_size);
    sequential_kh.create_spgemm_handle(KokkosSparse::SPGEMM_SERIAL);

    if (use_dynamic_scheduling) {
      sequential_kh.set_dynamic_scheduling(true);
    }

    spgemm_symbolic(&sequential_kh, m, n, k, crsMat.graph.row_map, crsMat.graph.entries, TRANSPOSEFIRST,
                    crsMat2.graph.row_map, crsMat2.graph.entries, TRANSPOSESECOND, row_mapC_ref);

    ExecSpace().fence();

    size_type c_nnz_size = sequential_kh.get_spgemm_handle()->get_c_nnz();
    entriesC_ref         = lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
    valuesC_ref          = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);

    spgemm_jacobi(&sequential_kh, m, n, k, crsMat.graph.row_map, crsMat.graph.entries, crsMat.values, TRANSPOSEFIRST,
                  crsMat2.graph.row_map, crsMat2.graph.entries, crsMat2.values, TRANSPOSESECOND, row_mapC_ref,
                  entriesC_ref, valuesC_ref, omega, dinv);

    ExecSpace().fence();

    Ccrsmat_ref = crsMat_t3("CorrectC", m, k, valuesC_ref.extent(0), valuesC_ref, row_mapC_ref, entriesC_ref);
  }

  for (int i = 0; i < repeat; ++i) {
    kh.create_spgemm_handle(KokkosSparse::SPGEMMAlgorithm(algorithm));

    // 250000 default. if cache-mode is used on KNL can increase to 1M.
    kh.get_spgemm_handle()->MaxColDenseAcc = params.MaxColDenseAcc;

    if (i == 0) {
      kh.get_spgemm_handle()->set_read_write_cost_calc(calculate_read_write_cost);
    }
    // do the compression whether in 2 step, or 1 step.
    kh.get_spgemm_handle()->set_compression_steps(!params.compression2step);
    // whether to scale the hash more. default is 1, so no scale.
    kh.get_spgemm_handle()->set_min_hash_size_scale(params.minhashscale);
    // max occupancy in 1-level LP hashes. LL hashes can be 100%
    kh.get_spgemm_handle()->set_first_level_hash_cut_off(params.first_level_hash_cut_off);
    // min reduction on FLOPs to run compression
    kh.get_spgemm_handle()->set_compression_cut_off(params.compression_cut_off);

    row_mapC = lno_view_t("non_const_lnow_row", m + 1);
    entriesC = lno_nnz_view_t("entriesC (empty)", 0);
    valuesC  = scalar_view_t("valuesC (empty)", 0);

    Kokkos::Timer timer1;
    spgemm_symbolic(&kh, m, n, k, crsMat.graph.row_map, crsMat.graph.entries, TRANSPOSEFIRST, crsMat2.graph.row_map,
                    crsMat2.graph.entries, TRANSPOSESECOND, row_mapC);

    ExecSpace().fence();
    double symbolic_time = timer1.seconds();

    Kokkos::Timer timer2;
    size_type c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (verbose) std::cout << "C SIZE:" << c_nnz_size << std::endl;
    if (c_nnz_size) {
      entriesC = lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
      valuesC  = scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);
    }

    spgemm_jacobi(&kh, m, n, k, crsMat.graph.row_map, crsMat.graph.entries, crsMat.values, TRANSPOSEFIRST,
                  crsMat2.graph.row_map, crsMat2.graph.entries, crsMat2.values, TRANSPOSESECOND, row_mapC, entriesC,
                  valuesC, omega, dinv);

    ExecSpace().fence();
    double numeric_time = timer2.seconds();

    std::cout << "mm_time:" << symbolic_time + numeric_time << " symbolic_time:" << symbolic_time
              << " numeric_time:" << numeric_time << std::endl;
  }
  if (verbose) {
    std::cout << "row_mapC:" << row_mapC.extent(0) << std::endl;
    std::cout << "entriesC:" << entriesC.extent(0) << std::endl;
    std::cout << "valuesC:" << valuesC.extent(0) << std::endl;
    KokkosKernels::Impl::print_1Dview(valuesC);
    KokkosKernels::Impl::print_1Dview(entriesC);
    KokkosKernels::Impl::print_1Dview(row_mapC);
  }
  crsMat_t3 Ccrsmat_result("CrsMatrixC", m, k, valuesC.extent(0), valuesC, row_mapC, entriesC);
  if (check_output) {
    bool is_identical = is_same_matrix<crsMat_t3, device_t>(Ccrsmat_result, Ccrsmat_ref);
    if (!is_identical) {
      std::cerr << "Result differs. If values are differing, might be floating "
                   "point order error."
                << std::endl;
      exit(1);
    }
  }
  return Ccrsmat_result;
}

template <typename size_type, typename lno_t, typename scalar_t, typename exec_space, typename hbm_mem_space,
          typename sbm_mem_space>
void run_spgemm_jacobi(Parameters params) {
  typedef exec_space myExecSpace;
  typedef Kokkos::Device<exec_space, hbm_mem_space> myFastDevice;
  typedef Kokkos::Device<exec_space, sbm_mem_space> mySlowExecSpace;

  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, myFastDevice, void, size_type> fast_crstmat_t;
  typedef typename KokkosSparse::CrsMatrix<scalar_t, lno_t, mySlowExecSpace, void, size_type> slow_crstmat_t;

  const char *a_mat_file = params.a_mtx_bin_file.c_str();
  const char *b_mat_file = params.b_mtx_bin_file.c_str();
  const char *c_mat_file = params.c_mtx_bin_file.c_str();

  slow_crstmat_t a_slow_crsmat, b_slow_crsmat, c_slow_crsmat;
  fast_crstmat_t a_fast_crsmat, b_fast_crsmat, c_fast_crsmat;

  // read a and b matrices and store them on slow or fast memory.

  if (params.a_mem_space == 1) {
    a_fast_crsmat = KokkosSparse::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(a_mat_file);
  } else {
    a_slow_crsmat = KokkosSparse::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(a_mat_file);
  }

  if ((b_mat_file == NULL || strcmp(b_mat_file, a_mat_file) == 0) && params.b_mem_space == params.a_mem_space) {
    std::cout << "Using A matrix for B as well" << std::endl;
    b_fast_crsmat = a_fast_crsmat;
    b_slow_crsmat = a_slow_crsmat;
  } else if (params.b_mem_space == 1) {
    if (b_mat_file == NULL) b_mat_file = a_mat_file;
    b_fast_crsmat = KokkosSparse::Impl::read_kokkos_crst_matrix<fast_crstmat_t>(b_mat_file);
  } else {
    if (b_mat_file == NULL) b_mat_file = a_mat_file;
    b_slow_crsmat = KokkosSparse::Impl::read_kokkos_crst_matrix<slow_crstmat_t>(b_mat_file);
  }

  if (params.a_mem_space == 1) {
    if (params.b_mem_space == 1) {
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, fast_crstmat_t,
                                                                    fast_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_fast_crsmat, b_fast_crsmat, params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, fast_crstmat_t,
                                                                    fast_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_fast_crsmat, b_fast_crsmat, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, fast_crstmat_t,
                                                                    slow_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_fast_crsmat, b_fast_crsmat, params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, fast_crstmat_t,
                                                                    slow_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_fast_crsmat, b_fast_crsmat, params);
        }
      }
    } else {
      // B is in slow memory
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, slow_crstmat_t,
                                                                    fast_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_fast_crsmat, b_slow_crsmat, params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, slow_crstmat_t,
                                                                    fast_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_fast_crsmat, b_slow_crsmat, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, slow_crstmat_t,
                                                                    slow_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_fast_crsmat, b_slow_crsmat, params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, fast_crstmat_t, slow_crstmat_t,
                                                                    slow_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_fast_crsmat, b_slow_crsmat, params);
        }
      }
    }
  } else {
    // A is in slow memory
    if (params.b_mem_space == 1) {
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, fast_crstmat_t,
                                                                    fast_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_slow_crsmat, b_fast_crsmat, params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, fast_crstmat_t,
                                                                    fast_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_slow_crsmat, b_fast_crsmat, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, fast_crstmat_t,
                                                                    slow_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_slow_crsmat, b_fast_crsmat, params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, fast_crstmat_t,
                                                                    slow_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_slow_crsmat, b_fast_crsmat, params);
        }
      }
    } else {
      // B is in slow memory
      if (params.c_mem_space == 1) {
        if (params.work_mem_space == 1) {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, slow_crstmat_t,
                                                                    fast_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_slow_crsmat, b_slow_crsmat, params);
        } else {
          c_fast_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, slow_crstmat_t,
                                                                    fast_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_slow_crsmat, b_slow_crsmat, params);
        }

      } else {
        // C is in slow memory.
        if (params.work_mem_space == 1) {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, slow_crstmat_t,
                                                                    slow_crstmat_t, hbm_mem_space, hbm_mem_space>(
              a_slow_crsmat, b_slow_crsmat, params);
        } else {
          c_slow_crsmat = KokkosKernels::Experiment::run_experiment<myExecSpace, slow_crstmat_t, slow_crstmat_t,
                                                                    slow_crstmat_t, sbm_mem_space, sbm_mem_space>(
              a_slow_crsmat, b_slow_crsmat, params);
        }
      }
    }
  }

  if (c_mat_file != NULL) {
    if (params.c_mem_space == 1) {
      KokkosSparse::sort_crs_matrix(c_fast_crsmat);

      KokkosSparse::Impl::write_graph_bin((lno_t)(c_fast_crsmat.numRows()),
                                          (size_type)(c_fast_crsmat.graph.entries.extent(0)),
                                          c_fast_crsmat.graph.row_map.data(), c_fast_crsmat.graph.entries.data(),
                                          c_fast_crsmat.values.data(), c_mat_file);
    } else {
      KokkosSparse::sort_crs_matrix(c_slow_crsmat);

      KokkosSparse::Impl::write_graph_bin((lno_t)c_slow_crsmat.numRows(),
                                          (size_type)c_slow_crsmat.graph.entries.extent(0),
                                          c_slow_crsmat.graph.row_map.data(), c_slow_crsmat.graph.entries.data(),
                                          c_slow_crsmat.values.data(), c_mat_file);
    }
  }
}

}  // namespace Experiment
}  // namespace KokkosKernels
