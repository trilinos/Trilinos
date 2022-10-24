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
// Questions? Contact Luc Berger-Vergiat (lberge@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cmath>

#include <limits>
#include <unordered_map>

#include <Kokkos_Core.hpp>
#include <KokkosSparse_BlockCrsMatrix.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_spmv.hpp>
#include "KokkosKernels_default_types.hpp"
#include <KokkosKernels_Test_Structured_Matrix.hpp>

namespace details {

enum class Implementation : int { KokkosKernels = 0, Cuda = 1, MKL = 2 };

///
/// Define default types
///
typedef double Scalar;
typedef int Ordinal;
///
//////////////////////////

/// Random generator
template <typename scalar_t>
inline scalar_t random() {
  auto const max = static_cast<scalar_t>(RAND_MAX) + static_cast<scalar_t>(1);
  return static_cast<scalar_t>(std::rand()) / max;
}

template <typename scalar_t>
inline void set_random_value(scalar_t &v) {
  v = random<scalar_t>();
}

template <typename scalar_t>
inline void set_random_value(Kokkos::complex<scalar_t> &v) {
  Scalar vre = random<scalar_t>();
  Scalar vim = random<scalar_t>();
  v          = Kokkos::complex<scalar_t>(vre, vim);
}

template <typename scalar_t>
inline void set_random_value(std::complex<scalar_t> &v) {
  scalar_t vre = random<scalar_t>();
  scalar_t vim = random<scalar_t>();
  v            = std::complex<scalar_t>(vre, vim);
}

template <typename scalar_t>
void make_block_entries(
    const KokkosSparse::CrsMatrix<scalar_t, Ordinal, Kokkos::HostSpace, void,
                                  size_t> &mat_b1,
    int blockSize, std::vector<Ordinal> &mat_rowmap,
    std::vector<Ordinal> &mat_colidx, std::vector<scalar_t> &mat_val) {
  Ordinal nRow = blockSize * mat_b1.numRows();
  size_t nnz = static_cast<size_t>(blockSize) * static_cast<size_t>(blockSize) *
               mat_b1.nnz();

  mat_val.resize(nnz);
  for (size_t ii = 0; ii < nnz; ++ii) set_random_value(mat_val[ii]);

  //
  // Create graph for CrsMatrix
  //

  mat_rowmap.assign(nRow + 1, 0);
  mat_colidx.assign(nnz, 0);

  for (Ordinal ir = 0; ir < mat_b1.numRows(); ++ir) {
    const auto jbeg = mat_b1.graph.row_map(ir);
    const auto jend = mat_b1.graph.row_map(ir + 1);
    for (Ordinal ib = 0; ib < blockSize; ++ib) {
      const Ordinal my_row   = ir * blockSize + ib;
      mat_rowmap[my_row + 1] = mat_rowmap[my_row] + (jend - jbeg) * blockSize;
      for (auto ijk = jbeg; ijk < jend; ++ijk) {
        const auto col0 = mat_b1.graph.entries(ijk);
        for (Ordinal jb = 0; jb < blockSize; ++jb) {
          mat_colidx[mat_rowmap[my_row] + (ijk - jbeg) * blockSize + jb] =
              col0 * blockSize + jb;
        }
      }
    }
  }  // for (Ordinal ir = 0; ir < mat_b1.numRows(); ++ir)
}

template <typename scalar_t>
int test_blockcrs_matrix_single_vec(
    const char fOp[],
    KokkosSparse::CrsMatrix<Scalar, Ordinal, Kokkos::HostSpace, void, size_t>
        mat_b1,
    int test, int loop, const scalar_t alpha, const scalar_t beta,
    const int bMax) {
  typedef typename KokkosSparse::CrsMatrix<
      scalar_t, Ordinal, Kokkos::DefaultExecutionSpace, void, size_t>
      crsMat_type;

  typedef typename crsMat_type::values_type::non_const_type scalar_view_t;
  typedef scalar_view_t x_vector_type;
  typedef scalar_view_t y_vector_type;

  srand(17312837);

  int num_errors    = 0;
  const auto bMax_o = static_cast<Ordinal>(bMax);
  for (Ordinal blockSize = 1; blockSize <= bMax_o; ++blockSize) {
    Ordinal nRow = blockSize * mat_b1.numRows();
    Ordinal nCol = nRow;
    std::vector<Ordinal> mat_rowmap;
    std::vector<Ordinal> mat_colidx;
    std::vector<scalar_t> mat_val;

    // Create the entries
    make_block_entries<scalar_t>(mat_b1, blockSize, mat_rowmap, mat_colidx,
                                 mat_val);

    // Create the CrsMatrix for the reference computation
    crsMat_type Acrs("new_crs_matr", nRow, nCol, mat_val.size(), &mat_val[0],
                     &mat_rowmap[0], &mat_colidx[0]);

    x_vector_type xref("new_right_hand_side", nRow);
    auto h_xref = Kokkos::create_mirror_view(xref);
    for (Ordinal ir = 0; ir < nRow; ++ir) {
      set_random_value(h_xref(ir));
    }
    Kokkos::deep_copy(xref, h_xref);

    y_vector_type y0("y_init", nRow);
    auto h_y0 = Kokkos::create_mirror_view(y0);
    for (Ordinal ir = 0; ir < nRow; ++ir) set_random_value(h_y0(ir));
    Kokkos::deep_copy(y0, h_y0);

    y_vector_type ycrs("crs_product_result", nRow);
    auto h_ycrs = Kokkos::create_mirror_view(ycrs);

    // Time a series of multiplications with the CrsMatrix
    double time_crs = 0.0;
    for (int jr = 0; jr < loop; ++jr) {
      for (Ordinal ir = 0; ir < nRow; ++ir) h_ycrs(ir) = h_y0(ir);
      Kokkos::deep_copy(ycrs, h_ycrs);
      Kokkos::Timer timer;
      KokkosSparse::spmv(fOp, alpha, Acrs, xref, beta, ycrs);
      time_crs += timer.seconds();
    }

    // Create the output vector
    y_vector_type yblockcrs("product_result", nRow);
    auto h_yblockcrs = Kokkos::create_mirror_view(yblockcrs);

    double time_blockcrs = 0.0;
    // Create the BlockCrsMatrix
    KokkosSparse::Experimental::BlockCrsMatrix<
        scalar_t, Ordinal, Kokkos::DefaultExecutionSpace, void, size_t>
        Ablockcrs(Acrs, blockSize);

    switch (static_cast<details::Implementation>(test)) {
      default:
      case Implementation::KokkosKernels: {
        // Time a series of multiplications with the BlockCrsMatrix
        for (int jr = 0; jr < loop; ++jr) {
          for (Ordinal ir = 0; ir < nRow; ++ir) h_yblockcrs(ir) = h_y0(ir);
          Kokkos::deep_copy(yblockcrs, h_yblockcrs);
          Kokkos::Timer timer;
          KokkosSparse::spmv(fOp, alpha, Ablockcrs, xref, beta, yblockcrs);
          time_blockcrs += timer.seconds();
        }
      } break;
    }

    // Check that the numerical result is matching
    Kokkos::deep_copy(h_ycrs, ycrs);
    Kokkos::deep_copy(h_yblockcrs, yblockcrs);
    double error = 0.0, maxNorm = 0.0;
    for (size_t ir = 0; ir < h_ycrs.extent(0); ++ir) {
      maxNorm = std::max(maxNorm, Kokkos::ArithTraits<Scalar>::abs(h_ycrs(ir)));
      error   = std::max(error, Kokkos::ArithTraits<Scalar>::abs(
                                  h_ycrs(ir) - h_yblockcrs(ir)));
    }

    double tol =
        (mat_val.size() / nRow) * std::numeric_limits<double>::epsilon();
    if (error > tol * maxNorm) {
      num_errors += 1;
      std::cout << static_cast<int>(test) << " ";
      std::cout << fOp << ", " << blockSize << " : "
                << " error " << error << " maxNorm " << maxNorm << " tol "
                << tol << " tol * maxNorm " << tol * maxNorm << "\n";
    }

    //-- Print the number of Gflops for both products
    if (blockSize == 1) {
      printf("Op, blockSize: AvgGFlop(CrsMatrix) AvgGFlop(BlockCrsMatrix) \n");
    }
    double num_flops     = mat_val.size() * 2 * loop;
    double crs_flop      = (num_flops / time_crs) * 1.0e-09;
    double blockcrs_flop = (num_flops / time_blockcrs) * 1.0e-09;
    std::cout << fOp << ", " << blockSize << "         : ";
    if (crs_flop < blockcrs_flop) {
      std::cout << crs_flop << "        <" << blockcrs_flop << ">";
    } else {
      std::cout << "<" << crs_flop << ">         " << blockcrs_flop;
    }
    std::cout << std::endl;

  }  // for (Ordinal blockSize = 1; blockSize < bMax; ++blockSize)

  return int(num_errors);
}

template <typename scalar_t>
int test_blockcrs_matrix_vec(
    const char fOp[],
    KokkosSparse::CrsMatrix<Scalar, Ordinal, Kokkos::HostSpace, void, size_t>
        mat_b1,
    int nvec, int test, int loop, const scalar_t alpha, const scalar_t beta,
    const int bMax) {
  typedef typename KokkosSparse::CrsMatrix<
      scalar_t, Ordinal, Kokkos::DefaultExecutionSpace, void, size_t>
      crsMat_type;

  typedef Kokkos::View<scalar_t **, Kokkos::LayoutLeft,
                       Kokkos::DefaultExecutionSpace>
      block_vector_t;

  srand(17312837);

  int num_errors    = 0;
  const auto bMax_o = static_cast<Ordinal>(bMax);
  for (Ordinal blockSize = 1; blockSize <= bMax_o; ++blockSize) {
    Ordinal nRow = blockSize * mat_b1.numRows();
    Ordinal nCol = nRow;
    std::vector<Ordinal> mat_rowmap;
    std::vector<Ordinal> mat_colidx;
    std::vector<scalar_t> mat_val;

    make_block_entries<scalar_t>(mat_b1, blockSize, mat_rowmap, mat_colidx,
                                 mat_val);

    // Create the CrsMatrix for the reference computation
    crsMat_type Acrs("new_crs_matr", nRow, nCol, mat_val.size(), &mat_val[0],
                     &mat_rowmap[0], &mat_colidx[0]);

    block_vector_t xref("new_right_hand_side", nRow, nvec);
    auto h_xref = Kokkos::create_mirror_view(xref);
    for (Ordinal jc = 0; jc < nvec; ++jc) {
      for (Ordinal ir = 0; ir < nRow; ++ir) {
        set_random_value(h_xref(ir, jc));
      }
    }
    Kokkos::deep_copy(xref, h_xref);

    block_vector_t y0("y_init", nRow, nvec);
    auto h_y0 = Kokkos::create_mirror_view(y0);
    for (Ordinal jc = 0; jc < nvec; ++jc)
      for (Ordinal ir = 0; ir < nRow; ++ir) set_random_value(h_y0(ir, jc));
    Kokkos::deep_copy(y0, h_y0);

    block_vector_t ycrs("crs_product_result", nRow, nvec);
    auto h_ycrs = Kokkos::create_mirror_view(ycrs);

    // Time a series of multiplications with the CrsMatrix format
    double time_crs = 0.0;
    for (int jr = 0; jr < loop; ++jr) {
      for (Ordinal jc = 0; jc < nvec; ++jc)
        for (Ordinal ir = 0; ir < nRow; ++ir) h_ycrs(ir, jc) = h_y0(ir, jc);
      Kokkos::deep_copy(ycrs, h_ycrs);
      Kokkos::Timer timer;
      KokkosSparse::spmv(fOp, alpha, Acrs, xref, beta, ycrs);
      time_crs += timer.seconds();
    }

    // Create the BlockCrsMatrix variable
    KokkosSparse::Experimental::BlockCrsMatrix<
        scalar_t, Ordinal, Kokkos::DefaultExecutionSpace, void, size_t>
        Ablockcrs(Acrs, blockSize);

    block_vector_t yblockcrs("blockcrs_product_result", nRow, nvec);
    auto h_yblockcrs = Kokkos::create_mirror_view(yblockcrs);

    // Time a series of multiplications with the BlockCrsMatrix
    double time_blockcrs = 0.0;
    switch (static_cast<details::Implementation>(test)) {
      default:
      case Implementation::KokkosKernels: {
        // Time a series of multiplications with the BlockCrsMatrix
        for (int jr = 0; jr < loop; ++jr) {
          for (Ordinal jc = 0; jc < nvec; ++jc) {
            for (Ordinal ir = 0; ir < nRow; ++ir)
              h_yblockcrs(ir, jc) = h_y0(ir, jc);
          }
          Kokkos::deep_copy(yblockcrs, h_yblockcrs);
          Kokkos::Timer timer;
          KokkosSparse::spmv(fOp, alpha, Ablockcrs, xref, beta, yblockcrs);
          time_blockcrs += timer.seconds();
        }
      } break;
    }

    // Check that the result is matching
    Kokkos::deep_copy(h_ycrs, ycrs);
    Kokkos::deep_copy(h_yblockcrs, yblockcrs);
    double tol =
        (mat_val.size() / nRow) * std::numeric_limits<double>::epsilon();
    for (int jc = 0; jc < nvec; ++jc) {
      double error = 0.0, maxNorm = 0.0;
      for (size_t ir = 0; ir < h_ycrs.extent(0); ++ir) {
        maxNorm =
            std::max(maxNorm, Kokkos::ArithTraits<Scalar>::abs(h_ycrs(ir, jc)));
        error = std::max(error, Kokkos::ArithTraits<Scalar>::abs(
                                    h_ycrs(ir, jc) - h_yblockcrs(ir, jc)));
      }
      if (error > tol * maxNorm) {
        num_errors += 1;
        std::cout << fOp << ", " << blockSize << " : rhs " << jc << " error "
                  << error << " maxNorm " << maxNorm << " tol " << tol
                  << " tol * maxNorm " << tol * maxNorm << "\n";
      }
    }

    // Print the number of Gflops
    if (blockSize == 1) {
      printf("Op, blockSize: AvgGFlop(CrsMatrix) AvgGFlop(BlockCrsMatrix) \n");
    }
    double num_flops     = mat_val.size() * 2 * loop * nvec;
    double crs_flop      = (num_flops / time_crs) * 1.0e-09;
    double blockcrs_flop = (num_flops / time_blockcrs) * 1.0e-09;
    std::cout << fOp << ", " << blockSize << "         ";
    if (crs_flop < blockcrs_flop) {
      // std::cout << crs_flop << "        <" << blockcrs_flop << ">";
      std::cout << crs_flop << "        " << blockcrs_flop << " ";
    } else {
      // std::cout << "<" << crs_flop << ">         " << blockcrs_flop;
      std::cout << " " << crs_flop << "         " << blockcrs_flop;
    }
    std::cout << std::endl;
  }

  return int(num_errors);
}

void print_help() {
  printf("BlockCrsMatrix SPMV benchmark code \n");
  printf("Options:\n");
  printf(
      "  -bs             : Maximum blocksize for the sparse matrix (default "
      "= "
      "16). \n");
  printf("  -h              : Help. \n");
  printf(
      "  -l [LOOP]       : How many spmv to run to aggregate average time "
      "(default = 512). \n");
  printf(
      "  -nx             : Number of points in the x-direction (default = "
      "32).\n");
  printf(
      "                    The matrix will be of dimension nx (nx - 1) (nx + "
      "1).\n");
  printf(
      "  -nv             : Number of vectors to multiply with (default = 1). "
      "\n");
  printf("  --op            : Use different operation \n");
  printf("                    Options: \n");
  printf("                    N = normal (default)  y <- alpha A x + beta y\n");
  printf(
      "                    C = conjugate         y <- alpha conj(A) x + beta "
      "y\n");
  printf(
      "                    T = transpose         y <- alpha A^T x + beta "
      "y\n");
  printf(
      "                    H = hermitian         y <- alpha A^H x + beta "
      "y\n");
}
}  // namespace details

int main(int argc, char **argv) {
  int loop = 512;
  int bMax = 16;
  int nvec = 1;
  int nx   = 32;

  char fOp[] = "N";

  int test = static_cast<int>(details::Implementation::KokkosKernels);

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-bs") == 0)) {
      int tmp = atoi(argv[++i]);
      bMax    = (tmp > 0) ? tmp : bMax;
      continue;
    }

    if ((strcmp(argv[i], "--tpl") == 0)) {
      i++;
      if ((strcmp(argv[i], "cuda") == 0))
        test = static_cast<int>(details::Implementation::Cuda);
      if ((strcmp(argv[i], "mkl") == 0))
        test = static_cast<int>(details::Implementation::MKL);
      continue;
    }

    if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
      details::print_help();
      return 0;
    }

    if ((strcmp(argv[i], "-l") == 0)) {
      int tmp = atoi(argv[++i]);
      loop    = (tmp > 0) ? tmp : loop;
      continue;
    }

    if ((strcmp(argv[i], "-nx") == 0)) {
      int tmp = atoi(argv[++i]);
      nx      = (tmp > 0) ? tmp : nx;
      continue;
    }

    if ((strcmp(argv[i], "-nv") == 0)) {
      int tmp = atoi(argv[++i]);
      nvec    = (tmp > 0) ? tmp : nvec;
      continue;
    }

    if ((strcmp(argv[i], "--op") == 0)) {
      i++;
      if ((strcmp(argv[i], "N") == 0)) strcpy(fOp, "N");
      if ((strcmp(argv[i], "C") == 0)) strcpy(fOp, "C");
      if ((strcmp(argv[i], "T") == 0)) strcpy(fOp, "T");
      if ((strcmp(argv[i], "H") == 0)) strcpy(fOp, "H");
      continue;
    }
  }

  Kokkos::initialize(argc, argv);
  {
    // The mat_structure view is used to generate a matrix using
    // finite difference (FD) or finite element (FE) discretization
    // on a cartesian grid.
    Kokkos::View<details::Ordinal * [3], Kokkos::HostSpace> mat_structure(
        "Matrix Structure", 3);
    mat_structure(0, 0) = nx;      // Request 8 grid point in 'x' direction
    mat_structure(0, 1) = 0;       // Add BC to the left
    mat_structure(0, 2) = 0;       // Add BC to the right
    mat_structure(1, 0) = nx - 1;  // Request 7 grid point in 'y' direction
    mat_structure(1, 1) = 0;       // Add BC to the bottom
    mat_structure(1, 2) = 0;       // Add BC to the top
    mat_structure(2, 0) = nx + 1;  // Request 9 grid point in 'z' direction
    mat_structure(2, 1) = 0;       // Add BC to the bottom
    mat_structure(2, 2) = 0;       // Add BC to the top

    typedef typename KokkosSparse::CrsMatrix<details::Scalar, details::Ordinal,
                                             Kokkos::HostSpace, void, size_t>
        h_crsMat_type;

    h_crsMat_type mat_b1 =
        Test::generate_structured_matrix3D<h_crsMat_type>("FD", mat_structure);

    int total_errors = 0;

    if (nvec == 1)
      total_errors = details::test_blockcrs_matrix_single_vec(
          fOp, mat_b1, test, loop, details::Scalar(3.1), details::Scalar(-2.4),
          bMax);
    else
      total_errors = details::test_blockcrs_matrix_vec(
          fOp, mat_b1, nvec, test, loop, details::Scalar(3.1),
          details::Scalar(-2.4), bMax);

    if (total_errors != 0) {
      printf("Kokkos::BlockCrsMatrix SpMV Test: Failed\n");
    }
  }
  Kokkos::finalize();
}
