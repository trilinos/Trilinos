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
#include <random>
#include <unordered_map>

#include <sstream>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_spmv.hpp>
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_IOUtils.hpp"

template <class matrix_type>
matrix_type generate_unbalanced_matrix(const typename matrix_type::ordinal_type numRows,
                                       const typename matrix_type::ordinal_type numEntries,
                                       const typename matrix_type::ordinal_type numLongRows,
                                       const typename matrix_type::ordinal_type numLongEntries) {
  using Scalar = typename matrix_type::value_type;
  using lno_t  = typename matrix_type::ordinal_type;

  using row_map_type = typename matrix_type::row_map_type::non_const_type;
  using entries_type = typename matrix_type::index_type::non_const_type;
  using values_type  = typename matrix_type::values_type::non_const_type;

  // Structure of the matrix:
  // the last numLongRows will contain the highly connected rows
  // the rest of the matrix will have a more classic sparse structure
  // with numEntries per row.

  // Randomly pick the length of the rows using a normal distribution
  std::mt19937 rand_generator(42);  // Seed with 42 for reproducibility
  std::normal_distribution<Scalar> row_dist{static_cast<Scalar>(numEntries),
                                            static_cast<Scalar>(std::sqrt(numEntries))};

  std::vector<lno_t> permutation(numRows - numLongRows);
  std::vector<lno_t> row_map_vec(numRows + 1);
  row_map_vec[0] = 0;
  for (lno_t rowIdx = 0; rowIdx < numRows - numLongRows; ++rowIdx) {
    row_map_vec[rowIdx + 1] = row_map_vec[rowIdx] + static_cast<lno_t>(row_dist(rand_generator));

    // also filling the permutation vector that will be used to construct long
    // rows
    permutation[rowIdx] = rowIdx;
  }

  std::normal_distribution<Scalar> long_row_dist{static_cast<Scalar>(numLongEntries),
                                                 static_cast<Scalar>(numLongEntries / 2)};
  lno_t rand_number;
  for (lno_t rowIdx = numRows - numLongRows; rowIdx < numRows; ++rowIdx) {
    rand_number             = static_cast<lno_t>(long_row_dist(rand_generator));
    row_map_vec[rowIdx + 1] = row_map_vec[rowIdx] + rand_number;
  }
  const lno_t numNNZ = row_map_vec[numRows];

  std::vector<lno_t> colind_vec(row_map_vec[numRows]);
  // We loop over the first part of the matrix and assume that the bandwidth is
  // 0.01*numRows i.e. highly concentrated around digaonal
  std::normal_distribution<Scalar> entry_dist{static_cast<Scalar>(0.0), static_cast<Scalar>(numRows / 100)};
  for (lno_t rowIdx = 0; rowIdx < numRows - numLongRows; ++rowIdx) {
    const lno_t rowLength = row_map_vec[rowIdx + 1] - row_map_vec[rowIdx];
    // Making the stencil symmetric because it looks a bit more like a regular
    // discretization
    for (lno_t entryIdx = 0; entryIdx < (rowLength / 2); ++entryIdx) {
      colind_vec[row_map_vec[rowIdx] + entryIdx]         = rowIdx - static_cast<lno_t>(entry_dist(rand_generator));
      colind_vec[row_map_vec[rowIdx + 1] - entryIdx - 1] = rowIdx + static_cast<lno_t>(entry_dist(rand_generator));
    }
    // Add diagonal entry if row length is an odd number
    if ((rowLength % 2) == 1) {
      colind_vec[row_map_vec[rowIdx] + rowLength / 2] = rowIdx;
    }
  }

  for (lno_t rowIdx = numRows - numLongRows; rowIdx < numRows; ++rowIdx) {
    // Generate a random permutation
    std::shuffle(permutation.begin(), permutation.end(), rand_generator);

    lno_t rowLength = row_map_vec[rowIdx + 1] - row_map_vec[rowIdx];
    for (lno_t entryIdx = 0; entryIdx < rowLength; ++entryIdx) {
      colind_vec[row_map_vec[rowIdx] + entryIdx] = permutation[entryIdx];
    }
  }

  row_map_type row_map("row map", numRows + 1);
  entries_type entries("entries", numNNZ);
  values_type values("values", numNNZ);

  // Copy row map values to view
  typename row_map_type::HostMirror row_map_h = Kokkos::create_mirror_view(row_map);
  row_map_h(0)                                = 0;
  for (lno_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
    row_map_h(rowIdx + 1) = row_map_vec[rowIdx + 1];
  }
  Kokkos::deep_copy(row_map, row_map_h);

  // Copy column indices to view
  typename row_map_type::HostMirror entries_h = Kokkos::create_mirror_view(entries);
  entries_h(0)                                = 0;
  for (lno_t entryIdx = 0; entryIdx < numNNZ; ++entryIdx) {
    entries_h(entryIdx) = colind_vec[entryIdx];
  }
  Kokkos::deep_copy(entries, entries_h);

  // Fill the values view with 1.0
  Kokkos::deep_copy(values, 1.0);

  matrix_type unbalanced_matrix("unbalanced matrix", numRows, numRows, numNNZ, values, row_map, entries);

  std::cout << std::endl;
  std::cout << "Matrix statistics:" << std::endl
            << "  - average nnz per row: " << row_map_vec[numRows - numLongRows] / (numRows - numLongRows) << std::endl;

  return unbalanced_matrix;
}

void print_help() {
  printf("SPMV merge benchmark code written by Luc Berger-Vergiat.\n");
  printf("The goal is to compare the merge path algorithm vs.\n");
  printf("TPLs and the KK native algorithm on imbalanced matrices.\n");
  printf("Options:\n");
  printf(
      "  --compare       : Compare the performance of the merge algo with the "
      "default algo.\n");
  printf("  -l [LOOP]       : How many spmv to run to aggregate average time. \n");
  printf("  -numRows        : Number of rows the matrix will contain.\n");
  printf("  -numEntries     : Number of entries per row.\n");
  printf(
      "  -numLongRows    : Number of rows that will contain more entries than "
      "the average.\n");
  printf(
      "  -numLongEntries : Number of entries per row in the unbalanced "
      "rows.\n");
}

int main(int argc, char** argv) {
  using Scalar = default_scalar;
  using lno_t  = default_lno_t;

  bool compare         = false;
  lno_t loop           = 100;
  lno_t numRows        = 175000;
  lno_t numEntries     = 15;
  lno_t numLongRows    = 4;
  lno_t numLongEntries = 30000;

  if (argc == 1) {
    print_help();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "--compare") == 0)) {
      compare = true;
      continue;
    }
    if ((strcmp(argv[i], "-l") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-numRows") == 0)) {
      numRows = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-numEntries") == 0)) {
      numEntries = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-numLongRows") == 0)) {
      numLongRows = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-numLongEntries") == 0)) {
      numLongEntries = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help();
      return 0;
    }
  }

  // We want an odd number of entries in all rows to generate a symmetric matrix
  if ((numEntries / 2) == 0) {
    ++numEntries;
  }
  if ((numLongEntries / 2) == 0) {
    ++numLongEntries;
  }

  std::cout << "Test parameters:" << std::endl
            << "  - loop:           " << loop << std::endl
            << "  - compare:        " << compare << std::endl
            << "  - numRows:        " << numRows << std::endl
            << "  - numEntries:     " << numEntries << std::endl
            << "  - numLongRows:    " << numLongRows << std::endl
            << "  - numLongEntries: " << numLongEntries << std::endl;

  Kokkos::initialize(argc, argv);

  {
    // Note that we template the matrix with entries=lno_t and offsets=lno_t
    // so that TPLs can be used
    using matrix_type = KokkosSparse::CrsMatrix<Scalar, lno_t, Kokkos::DefaultExecutionSpace, void, lno_t>;
    using values_type = typename matrix_type::values_type::non_const_type;
    using handle_type = KokkosSparse::SPMVHandle<Kokkos::DefaultExecutionSpace, matrix_type, values_type, values_type>;
    const Scalar SC_ONE = Kokkos::ArithTraits<Scalar>::one();
    const Scalar alpha  = SC_ONE + SC_ONE;
    const Scalar beta   = alpha + SC_ONE;

    matrix_type test_matrix = generate_unbalanced_matrix<matrix_type>(numRows, numEntries, numLongRows, numLongEntries);

    values_type y("right hand side", test_matrix.numRows());
    values_type x("left hand side", test_matrix.numCols());
    Kokkos::deep_copy(x, SC_ONE);
    Kokkos::deep_copy(y, SC_ONE);

    handle_type handleMerge(KokkosSparse::SPMV_MERGE_PATH);

    // Perform a so called "warm-up" run
    KokkosSparse::spmv(&handleMerge, "N", alpha, test_matrix, x, beta, y);

    double min_time = 1.0e32, max_time = 0.0, avg_time = 0.0;
    for (int iterIdx = 0; iterIdx < loop; ++iterIdx) {
      Kokkos::Timer timer;
      KokkosSparse::spmv(&handleMerge, "N", alpha, test_matrix, x, beta, y);
      Kokkos::fence();
      double time = timer.seconds();
      avg_time += time;
      if (time > max_time) max_time = time;
      if (time < min_time) min_time = time;
    }

    std::cout << "KK Merge alg  ---  min: " << min_time << " max: " << max_time << " avg: " << avg_time / loop
              << std::endl;

    // Run the cusparse default algorithm and native kokkos-kernels algorithm
    // then output timings for comparison
    if (compare) {
      handle_type handleDefault;
      // Warm up
      KokkosSparse::spmv(&handleDefault, "N", alpha, test_matrix, x, beta, y);

      min_time = 1.0e32;
      max_time = 0.0;
      avg_time = 0.0;
      for (int iterIdx = 0; iterIdx < loop; ++iterIdx) {
        Kokkos::Timer timer;
        KokkosSparse::spmv(&handleDefault, "N", alpha, test_matrix, x, beta, y);
        Kokkos::fence();
        double time = timer.seconds();
        avg_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }

      std::cout << "Default alg   ---  min: " << min_time << " max: " << max_time << " avg: " << avg_time / loop
                << std::endl;

      handle_type handleNative(KokkosSparse::SPMV_NATIVE);
      KokkosSparse::spmv(&handleNative, "N", alpha, test_matrix, x, beta, y);

      min_time = 1.0e32;
      max_time = 0.0;
      avg_time = 0.0;
      for (int iterIdx = 0; iterIdx < loop; ++iterIdx) {
        Kokkos::Timer timer;
        KokkosSparse::spmv(&handleNative, "N", alpha, test_matrix, x, beta, y);
        Kokkos::fence();
        double time = timer.seconds();
        avg_time += time;
        if (time > max_time) max_time = time;
        if (time < min_time) min_time = time;
      }

      std::cout << "KK Native alg ---  min: " << min_time << " max: " << max_time << " avg: " << avg_time / loop
                << std::endl;
    }
  }

  Kokkos::finalize();
}  // main
