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

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>
#include <unordered_map>
#include <KokkosSparse_spmv_test.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_IOUtils.hpp>
#include <KokkosSparse_spmv.hpp>
#include "KokkosKernels_default_types.hpp"
#include <spmv/KokkosKernels_spmv_data.hpp>
#include <spmv/Kokkos_SPMV.hpp>
#include <spmv/Kokkos_SPMV_Inspector.hpp>

#ifdef HAVE_CUSPARSE
#include <CuSparse_SPMV.hpp>
#endif

#ifdef HAVE_MKL
#include <MKL_SPMV.hpp>
#endif

#ifdef KOKKOS_ENABLE_OPENMP
#include <OpenMPStatic_SPMV.hpp>
#include <OpenMPDynamic_SPMV.hpp>
#include <OpenMPSmartStatic_SPMV.hpp>
#endif

int test_crs_matrix_singlevec(Ordinal numRows, Ordinal numCols, int test, const char* filename, Ordinal rows_per_thread,
                              int team_size, int vector_length, int schedule, int loop) {
  typedef KokkosSparse::CrsMatrix<Scalar, Ordinal, Kokkos::DefaultExecutionSpace, void, Offset> matrix_type;

  spmv_additional_data data(test);

  std::cout << "running CRS matrix single vec" << std::endl;

  srand(17312837);
  matrix_type A;
  if (filename)
    A = KokkosSparse::Impl::read_kokkos_crst_matrix<matrix_type>(filename);
  else {
    Offset nnz = 10 * numRows;
    // note: the help text says the bandwidth is fixed at 0.01 * numRows
    // CAVEAT:  small problem sizes are problematic, b/c of 0.01*numRows
    A = KokkosSparse::Impl::kk_generate_sparse_matrix<matrix_type>(numRows, numCols, nnz, 0, 0.01 * numRows);
  }
  SPMVTestData test_data = setup_test(&data, A, rows_per_thread, team_size, vector_length, schedule, loop);
  for (int i = 0; i < loop; i++) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_ARMPL
    if (test == ARMPL) {
      if (std::is_same<Scalar, double>::value || std::is_same<Scalar, float>::value) {
        data.set_armpl_spmat(test_data.numRows, test_data.numCols, test_data.A.graph.row_map.data(),
                             test_data.A.graph.entries.data(), test_data.A.values.data());
      } else {
        throw std::runtime_error(
            "Can't use ArmPL mat-vec for scalar types other than double and "
            "float.");
      }
    }
#endif
    run_benchmark(test_data);
  }

  // Performance Output
  double matrix_size =
      1.0 * ((test_data.nnz * (sizeof(Scalar) + sizeof(Ordinal)) + numRows * sizeof(Offset))) / 1024 / 1024;
  double vector_size      = 2.0 * numRows * sizeof(Scalar) / 1024 / 1024;
  double vector_readwrite = (test_data.nnz + numCols) * sizeof(Scalar) / 1024 / 1024;

  double problem_size = matrix_size + vector_size;
  printf(
      "NNZ NumRows NumCols ProblemSize(MB) AveBandwidth(GB/s) "
      "MinBandwidth(GB/s) MaxBandwidth(GB/s) AveGFlop MinGFlop MaxGFlop "
      "aveTime(ms) maxTime(ms) minTime(ms) numErrors\n");
  printf(
      "%i %i %i %6.2lf ( %6.2lf %6.2lf %6.2lf ) ( %6.3lf %6.3lf %6.3lf ) ( "
      "%6.3lf %6.3lf %6.3lf ) %i RESULT\n",
      test_data.nnz, numRows, numCols, problem_size,
      (matrix_size + vector_readwrite) / test_data.ave_time * loop / 1024,
      (matrix_size + vector_readwrite) / test_data.max_time / 1024,
      (matrix_size + vector_readwrite) / test_data.min_time / 1024,
      2.0 * test_data.nnz * loop / test_data.ave_time / 1e9, 2.0 * test_data.nnz / test_data.max_time / 1e9,
      2.0 * test_data.nnz / test_data.min_time / 1e9, test_data.ave_time / loop * 1000, test_data.max_time * 1000,
      test_data.min_time * 1000, test_data.num_errors);
  return (int)test_data.total_error;
}

void print_help() {
  printf("SPMV benchmark code written by Christian Trott.\n");
  printf(
      "OpenMP implementations written by Simon Hammond (Sandia National "
      "Laboratories).\n\n");
  printf("Options:\n");
  printf(
      "  -s [N]          : generate a semi-random banded (band size 0.01xN) "
      "NxN matrix\n");
  printf("                    with average of 10 entries per row.\n");
  printf("  --test [OPTION] : Use different kernel implementations\n");
  printf("                    Options:\n");
  printf("                      kk,kk-kernels          (Kokkos/Trilinos)\n");
  printf(
      "                      kk-insp                (Kokkos Structure "
      "Inspection)\n");
#ifdef KOKKOS_ENABLE_OPENMP
  printf("                      omp-dynamic,omp-static (Standard OpenMP)\n");
  printf(
      "                      omp-insp               (OpenMP Structure "
      "Inspection)\n");
#endif
  printf("                      mkl, armpl,cusparse    (Vendor Libraries)\n\n");
  printf(
      "  --schedule [SCH]: Set schedule for kk variant (static,dynamic,auto [ "
      "default ]).\n");
  printf(
      "  -f [file]       : Read in Matrix Market formatted text file "
      "'file'.\n");
  printf("  -fb [file]      : Read in binary Matrix files 'file'.\n");
  printf("  --write-binary  : In combination with -f, generate binary files.\n");
  printf("  --offset [O]    : Subtract O from every index.\n");
  printf(
      "                    Useful in case the matrix market file is not 0 "
      "based.\n\n");
  printf("  -rpt [K]        : Number of Rows assigned to a thread.\n");
  printf("  -ts [T]         : Number of threads per team.\n");
  printf(
      "  -vl [V]         : Vector-length (i.e. how many Cuda threads are a "
      "Kokkos 'thread').\n");
  printf("  -l [LOOP]       : How many spmv to run to aggregate average time. \n");
}

int main(int argc, char** argv) {
  long long int size = 110503;  // a prime number
  // int numVecs = 4;
  int test = KOKKOS;
  // int type=-1;
  char* filename = NULL;

  int rows_per_thread = -1;
  int vector_length   = -1;
  int team_size       = -1;
  int schedule        = AUTO;
  int loop            = 100;

  if (argc == 1) {
    print_help();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-s") == 0)) {
      size = atoi(argv[++i]);
      continue;
    }

    if ((strcmp(argv[i], "--test") == 0)) {
      i++;
      if (i == argc) {
        std::cerr << "Must pass algorithm name after '--test'";
        exit(1);
      }

      if ((strcmp(argv[i], "mkl") == 0)) test = MKL;
      if ((strcmp(argv[i], "armpl") == 0)) test = ARMPL;
      if ((strcmp(argv[i], "kk") == 0)) test = KOKKOS;
      if ((strcmp(argv[i], "cusparse") == 0)) test = CUSPARSE;
      if ((strcmp(argv[i], "kk-kernels") == 0)) test = KK_KERNELS;
      if ((strcmp(argv[i], "kk-kernels-insp") == 0)) test = KK_KERNELS_INSP;
      if ((strcmp(argv[i], "kk-insp") == 0)) test = KK_INSP;
#ifdef KOKKOS_ENABLE_OPENMP
      if ((strcmp(argv[i], "omp-static") == 0)) test = OMP_STATIC;
      if ((strcmp(argv[i], "omp-dynamic") == 0)) test = OMP_DYNAMIC;
      if ((strcmp(argv[i], "omp-insp") == 0)) test = OMP_INSP;
#endif
      continue;
    }
    // if((strcmp(argv[i],"--type")==0)) {type=atoi(argv[++i]); continue;}
    if ((strcmp(argv[i], "-f") == 0)) {
      filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-fb") == 0)) {
      filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-rpt") == 0)) {
      rows_per_thread = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-ts") == 0)) {
      team_size = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-vl") == 0)) {
      vector_length = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-l") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--schedule") == 0)) {
      i++;
      if ((strcmp(argv[i], "auto") == 0)) schedule = AUTO;
      if ((strcmp(argv[i], "dynamic") == 0)) schedule = DYNAMIC;
      if ((strcmp(argv[i], "static") == 0)) schedule = STATIC;
      continue;
    }
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help();
      return 0;
    }
  }

  Kokkos::initialize(argc, argv);

  int total_errors =
      test_crs_matrix_singlevec(size, size, test, filename, rows_per_thread, team_size, vector_length, schedule, loop);

  if (total_errors == 0)
    printf("Kokkos::MultiVector Test: Passed\n");
  else
    printf("Kokkos::MultiVector Test: Failed\n");

  Kokkos::finalize();
}
