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

#include <KokkosKernels_config.h>
#include "KokkosSparse_pcg.hpp"

#include "KokkosKernels_Utils.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_TestUtils.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include <iostream>

#define MAXVAL 1

template <typename scalar_view_t>
scalar_view_t create_x_vector(default_lno_t nv, default_scalar max_value = 1.0) {
  scalar_view_t kok_x("X", nv);

  typename scalar_view_t::HostMirror h_x = Kokkos::create_mirror_view(kok_x);

  for (default_lno_t i = 0; i < nv; ++i) {
    default_scalar r = static_cast<default_scalar>(rand()) / static_cast<default_scalar>(RAND_MAX / max_value);
    h_x(i)           = r;
  }
  Kokkos::deep_copy(kok_x, h_x);
  return kok_x;
}

template <typename crsMat_t, typename vector_t>
vector_t create_y_vector(crsMat_t crsMat, vector_t x_vector) {
  vector_t y_vector("Y VECTOR", crsMat.numRows());
  KokkosSparse::spmv("N", 1, crsMat, x_vector, 1, y_vector);
  return y_vector;
}

template <typename ExecSpace, typename crsMat_t>
void run_experiment(crsMat_t crsmat, int clusterSize, bool useSequential) {
  typedef typename crsMat_t::values_type::non_const_type scalar_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;

  typedef typename lno_nnz_view_t::value_type lno_t;
  typedef typename lno_view_t::value_type size_type;
  typedef typename scalar_view_t::value_type scalar_t;

  default_lno_t nv             = crsmat.numRows();
  scalar_view_t kok_x_original = create_x_vector<scalar_view_t>(nv, MAXVAL);
  scalar_view_t kok_b_vector   = create_y_vector(crsmat, kok_x_original);

  // create X vector
  scalar_view_t kok_x_vector("kok_x_vector", nv);

  double solve_time                   = 0;
  const unsigned cg_iteration_limit   = 100000;
  const double cg_iteration_tolerance = 1e-7;

  KokkosKernels::Experimental::Example::CGSolveResult cg_result;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle<size_type, lno_t, scalar_t, ExecSpace, ExecSpace, ExecSpace>
      KernelHandle;

  KernelHandle kh;

  if (clusterSize == 1)
    kh.create_gs_handle();
  else
    kh.create_gs_handle(KokkosSparse::CLUSTER_BALLOON, clusterSize);
  Kokkos::Timer timer1;
  KokkosKernels::Experimental::Example::pcgsolve(kh, crsmat, kok_b_vector, kok_x_vector, cg_iteration_limit,
                                                 cg_iteration_tolerance, &cg_result, true, clusterSize, useSequential);
  Kokkos::fence();

  solve_time = timer1.seconds();

  std::string algoSummary;
  if (useSequential)
    algoSummary = "SEQUENTIAL SGS";
  else {
    if (clusterSize == 1)
      algoSummary = "POINT-COLORING SGS";
    else
      algoSummary = "CLUSTER-COLORING SGS (CLUSTER SIZE " + std::to_string(clusterSize) + ")";
  }

  std::cout << "DEFAULT SOLVE: " << algoSummary << " PRECONDITIONER"
            << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
            << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
            << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
            << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
            << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
            << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time << "]"
            << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time / (cg_result.iteration + 1) << "]"
            << "\n\tSOLVE_TIME                  [" << solve_time << "]" << std::endl;

  /*
  kh.destroy_gs_handle();
  kh.create_gs_handle(KokkosKernels::Experimental::Graph::GS_PERMUTED);

  kok_x_vector = scalar_view_t("kok_x_vector", nv);
  timer1.reset();
  KokkosKernels::Experimental::Example::pcgsolve(
        kh
      , crsmat
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , true
  );

  Kokkos::fence();
  solve_time = timer1.seconds();
  std::cout  << "\nPERMUTED SGS SOLVE:"
      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time <<
  "]"
      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time /
  (cg_result.iteration  + 1) << "]"
      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
      << std::endl ;


  kh.destroy_gs_handle();
  kh.create_gs_handle(KokkosKernels::Experimental::Graph::GS_TEAM);

  kok_x_vector = scalar_view_t("kok_x_vector", nv);
  timer1.reset();
  KokkosKernels::Experimental::Example::pcgsolve(
        kh
      , crsmat
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , true
  );
  Kokkos::fence();

  solve_time = timer1.seconds();
  std::cout  << "\nTEAM SGS SOLVE:"
      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time <<
  "]"
      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time /
  (cg_result.iteration  + 1) << "]"
      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
      << std::endl ;




  kok_x_vector = scalar_view_t("kok_x_vector", nv);
  timer1.reset();
  KokkosKernels::Experimental::Example::pcgsolve(
        kh
      , crsmat
      , kok_b_vector
      , kok_x_vector
      , cg_iteration_limit
      , cg_iteration_tolerance
      , & cg_result
      , false
  );
  Kokkos::fence();

  solve_time = timer1.seconds();
  std::cout  << "\nCG SOLVE (With no Preconditioner):"
      << "\n\t(P)CG_NUM_ITER              [" << cg_result.iteration << "]"
      << "\n\tMATVEC_TIME                 [" << cg_result.matvec_time << "]"
      << "\n\tCG_RESIDUAL                 [" << cg_result.norm_res << "]"
      << "\n\tCG_ITERATION_TIME           [" << cg_result.iter_time << "]"
      << "\n\tPRECONDITIONER_TIME         [" << cg_result.precond_time << "]"
      << "\n\tPRECONDITIONER_INIT_TIME    [" << cg_result.precond_init_time <<
  "]"
      << "\n\tPRECOND_APPLY_TIME_PER_ITER [" << cg_result.precond_time /
  (cg_result.iteration  + 1) << "]"
      << "\n\tSOLVE_TIME                  [" << solve_time<< "]"
      << std::endl ;
  */
}

enum {
  CMD_USE_THREADS = 0,
  CMD_USE_NUMA,
  CMD_USE_CORE_PER_NUMA,
  CMD_USE_CUDA,
  CMD_USE_HIP,
  CMD_USE_OPENMP,
  CMD_DEVICE,
  CMD_BIN_MTX,
  CMD_CLUSTER_SIZE,
  CMD_USE_SEQUENTIAL_SGS,
  CMD_ERROR,
  CMD_COUNT
};

template <typename execution_space>
void run_pcg(int *cmdline, const char *mtx_file) {
  default_lno_t nv = 0, ne = 0;
  default_lno_t *xadj, *adj;
  default_scalar *ew;

  KokkosSparse::Impl::read_matrix<default_lno_t, default_lno_t, default_scalar>(&nv, &ne, &xadj, &adj, &ew, mtx_file);

  typedef typename KokkosSparse::CrsMatrix<default_scalar, default_lno_t, execution_space, void, default_size_type>
      crsMat_t;

  typedef typename crsMat_t::StaticCrsGraphType graph_t;
  typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
  typedef typename crsMat_t::index_type::non_const_type cols_view_t;
  typedef typename crsMat_t::values_type::non_const_type values_view_t;

  row_map_view_t rowmap_view("rowmap_view", nv + 1);
  cols_view_t columns_view("colsmap_view", ne);
  values_view_t values_view("values_view", ne);

  {
    typename row_map_view_t::HostMirror hr = Kokkos::create_mirror_view(rowmap_view);
    typename cols_view_t::HostMirror hc    = Kokkos::create_mirror_view(columns_view);
    typename values_view_t::HostMirror hv  = Kokkos::create_mirror_view(values_view);

    for (default_lno_t i = 0; i <= nv; ++i) {
      hr(i) = xadj[i];
    }

    for (default_lno_t i = 0; i < ne; ++i) {
      hc(i) = adj[i];
      hv(i) = ew[i];
    }
    Kokkos::deep_copy(rowmap_view, hr);
    Kokkos::deep_copy(columns_view, hc);
    Kokkos::deep_copy(values_view, hv);
  }
  graph_t static_graph(columns_view, rowmap_view);
  crsMat_t crsmat("CrsMatrix", nv, values_view, static_graph);

  delete[] xadj;
  delete[] adj;
  delete[] ew;

  run_experiment<execution_space, crsMat_t>(crsmat, cmdline[CMD_CLUSTER_SIZE], cmdline[CMD_USE_SEQUENTIAL_SGS]);
}

int main(int argc, char **argv) {
  int cmdline[CMD_COUNT];
  char *mtx_file = NULL;
  for (int i = 0; i < CMD_COUNT; ++i) cmdline[i] = 0;

  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      cmdline[CMD_USE_THREADS] = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      cmdline[CMD_USE_OPENMP] = atoi(argv[++i]);
    }
    /*
    else if ( 0 == Test::string_compare_no_case( argv[i] , "--cores" ) ) {
      //Note BMK: specifying #NUMA regions isn't supported by initialize
      sscanf( argv[++i] , "%dx%d" ,
              cmdline + CMD_USE_NUMA ,
              cmdline + CMD_USE_CORE_PER_NUMA );
    }
    */
    else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      cmdline[CMD_USE_CUDA] = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--hip")) {
      cmdline[CMD_USE_HIP] = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--device-id")) {
      cmdline[CMD_DEVICE] = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--cluster-size")) {
      cmdline[CMD_CLUSTER_SIZE] = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--seq-gs")) {
      cmdline[CMD_USE_SEQUENTIAL_SGS] = 1;
    }

    else if (0 == Test::string_compare_no_case(argv[i], "--mtx")) {
      mtx_file = argv[++i];
    } else {
      cmdline[CMD_ERROR] = 1;
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      std::cerr << "OPTIONS\n\t--threads [numThreads]\n\t--openmp "
                   "[numThreads]\n\t--cuda\n\t--hip\n\t--device-id[DeviceIndex]"
                   "\n\t--mtx[binary_mtx_file]"
                << std::endl;

      return 0;
    }
  }
  // default cluster size is always 1 (this runs point coloring GS)
  if (cmdline[CMD_CLUSTER_SIZE] == 0) cmdline[CMD_CLUSTER_SIZE] = 1;

  if (mtx_file == NULL) {
    std::cerr << "Provide a matrix file" << std::endl;
    std::cerr << "OPTIONS\n\t--threads [numThreads]\n\t--openmp "
                 "[numThreads]\n\t--cuda\n\t--hip\n\t--device-id[DeviceIndex]"
                 "\n\t--mtx[matrix]"
              << std::endl;

    return 0;
  }

  // Construct with default args, change members based on exec space
  Kokkos::InitializationSettings init_args;

  init_args.set_device_id(cmdline[CMD_DEVICE]);
  init_args.set_num_threads(std::max(cmdline[CMD_USE_THREADS], cmdline[CMD_USE_OPENMP]));
  if (cmdline[CMD_USE_NUMA] && cmdline[CMD_USE_CORE_PER_NUMA]) {
    KokkosKernels::Impl::throw_runtime_exception("NUMA init arg is no longer supported by Kokkos");
    // init_args.num_numa = cmdline[CMD_USE_NUMA];
  }

  Kokkos::initialize(init_args);
  {
#if defined(KOKKOS_ENABLE_THREADS)
    if (cmdline[CMD_USE_THREADS]) run_pcg<Kokkos::Threads>(cmdline, mtx_file);
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
    if (cmdline[CMD_USE_OPENMP]) run_pcg<Kokkos::OpenMP>(cmdline, mtx_file);
#endif
#if defined(KOKKOS_ENABLE_CUDA)
    if (cmdline[CMD_USE_CUDA]) run_pcg<Kokkos::Cuda>(cmdline, mtx_file);
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (cmdline[CMD_USE_HIP]) run_pcg<Kokkos::HIP>(cmdline, mtx_file);
#endif
  }
  Kokkos::finalize();
  return 0;
}
