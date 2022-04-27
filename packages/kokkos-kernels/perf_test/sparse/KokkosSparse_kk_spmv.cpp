/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>
#include <unordered_map>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_spmv.hpp>
#include "KokkosKernels_default_types.hpp"

typedef default_scalar Scalar;
typedef default_lno_t Ordinal;
typedef default_size_type Offset;

template <typename Layout>
void run_spmv(Ordinal numRows, Ordinal numCols, const char* filename, int loop,
              int num_vecs, char mode, Scalar beta) {
  typedef KokkosSparse::CrsMatrix<Scalar, Ordinal,
                                  Kokkos::DefaultExecutionSpace, void, Offset>
      matrix_type;
  typedef typename Kokkos::View<Scalar**, Layout> mv_type;
  typedef typename mv_type::HostMirror h_mv_type;

  srand(17312837);
  matrix_type A;
  if (filename)
    A = KokkosKernels::Impl::read_kokkos_crst_matrix<matrix_type>(filename);
  else {
    Offset nnz = 10 * numRows;
    // note: the help text says the bandwidth is fixed at 0.01 * numRows
    A = KokkosKernels::Impl::kk_generate_sparse_matrix<matrix_type>(
        numRows, numCols, nnz, 0, 0.01 * numRows);
  }
  numRows = A.numRows();
  numCols = A.numCols();
  mv_type x("X", numCols, num_vecs);
  mv_type y("Y", numRows, num_vecs);
  h_mv_type h_x         = Kokkos::create_mirror_view(x);
  h_mv_type h_y         = Kokkos::create_mirror_view(y);
  h_mv_type h_y_compare = Kokkos::create_mirror(y);

  for (int v = 0; v < num_vecs; v++) {
    for (int i = 0; i < numCols; i++) {
      h_x(i, v) = (Scalar)(1.0 * (rand() % 40) - 20.);
    }
  }

  Kokkos::deep_copy(x, h_x);

  // Benchmark
  auto x0 = Kokkos::subview(x, Kokkos::ALL(), 0);
  auto y0 = Kokkos::subview(y, Kokkos::ALL(), 0);
  Kokkos::Timer timer;
  for (int i = 0; i < loop; i++) {
    if (num_vecs == 1) {
      // run the rank-1 version
      KokkosSparse::spmv(&mode, 1.0, A, x0, beta, y0);
    } else {
      // rank-2
      KokkosSparse::spmv(&mode, 1.0, A, x, beta, y);
    }
    Kokkos::DefaultExecutionSpace().fence();
  }
  double avg_time = timer.seconds() / loop;
  std::cout << avg_time << " s\n";
}

void print_help() {
  printf("  -s [nrows]            : matrix dimension (square)\n");
  printf(
      "  --nv n                : number of columns in x/y multivector (default "
      "1).\n");
  printf(
      "  --layout left|right   : memory layout of x/y. Default depends on "
      "build's default execution space\n");
  printf(
      "  -m N|T                : matrix apply mode: N (normal, default), T "
      "(transpose)\n");
  printf(
      "  -f [file],-fb [file]  : Read in Matrix Market (.mtx), or binary "
      "(.bin) matrix file.\n");
  printf(
      "  -l [LOOP]             : How many spmv to run to aggregate average "
      "time. \n");
  printf("  -b beta               : beta, as in y := Ax + (beta)y\n");
}

int main(int argc, char** argv) {
  long long int size = 110503;  // a prime number
  char* filename     = NULL;

  char mode = 'N';
  char layout;
  if (std::is_same<default_layout, Kokkos::LayoutLeft>::value)
    layout = 'L';
  else
    layout = 'R';
  int loop     = 100;
  int num_vecs = 1;
  Scalar beta  = 0.0;

  if (argc == 1) {
    print_help();
    return 0;
  }

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-s") == 0)) {
      size = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-fb") == 0)) {
      filename = argv[++i];
      continue;
    }
    if ((strcmp(argv[i], "-l") == 0)) {
      loop = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-m") == 0)) {
      mode = toupper(argv[++i][0]);
      continue;
    }
    if ((strcmp(argv[i], "--nv") == 0)) {
      num_vecs = atoi(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "-b") == 0)) {
      beta = atof(argv[++i]);
      continue;
    }
    if ((strcmp(argv[i], "--layout") == 0)) {
      i++;
      if (toupper(argv[i][0]) == 'L')
        layout = 'L';
      else if (toupper(argv[i][0]) == 'R')
        layout = 'R';
      else
        throw std::runtime_error("Invalid layout");
    }
    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      print_help();
      return 0;
    }
  }

  Kokkos::initialize(argc, argv);

  std::cout << size << " rows/cols, mode " << mode << ", " << num_vecs
            << " vectors, beta = " << beta << ", layout " << layout << ": ";

  if (layout == 'L')
    run_spmv<Kokkos::LayoutLeft>(size, size, filename, loop, num_vecs, mode,
                                 beta);
  else
    run_spmv<Kokkos::LayoutRight>(size, size, filename, loop, num_vecs, mode,
                                  beta);

  Kokkos::finalize();
}
