// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include "Kokkos_Core.hpp"

#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_SortCrs.hpp"
#include "KokkosSparse_spiluk.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_nrm2.hpp"

int main(int argc, char *argv[]) {
  using Scalar  = KokkosKernels::default_scalar;
  using Ordinal = KokkosKernels::default_lno_t;
  using Offset  = KokkosKernels::default_size_type;

  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using MemSpace  = typename ExecSpace::memory_space;
  using Device    = Kokkos::Device<ExecSpace, MemSpace>;
  using Matrix    = KokkosSparse::CrsMatrix<Scalar, Ordinal, Device, void, Offset>;
  using RowMap    = typename Matrix::StaticCrsGraphType::row_map_type::non_const_type;
  using Entries   = typename Matrix::StaticCrsGraphType::entries_type::non_const_type;
  using Values    = typename Matrix::values_type::non_const_type;
  using Handle =
      KokkosKernels::Experimental::KokkosKernelsHandle<Offset, Ordinal, Scalar, ExecSpace, MemSpace, MemSpace>;

  std::string filename;
  Ordinal fillLevel = 2;

  for (int i = 1; i < argc; ++i) {
    const std::string token = argv[i];
    if (token == "--filename" || token == "-f") {
      filename = argv[++i];
    } else if (token == "--fill-level" || token == "-l") {
      fillLevel = std::atoi(argv[++i]);
    } else if (token == "--help" || token == "-h") {
      std::cout << "KokkosSparse spiluk example options:\n"
                << "  --filename,   -f <file>  : Path to a MatrixMarket (.mtx) sparse matrix file\n"
                << "                             (if omitted, a synthetic matrix is generated)\n"
                << "  --fill-level, -l <level> : ILU(k) fill level (default: 2)\n"
                << "  --help,       -h         : Print this message\n"
                << "Example: ./KokkosSparse_example_spiluk -f mymatrix.mtx -l 3\n";
      return 0;
    }
  }

  const Scalar zero = KokkosKernels::ArithTraits<Scalar>::zero();
  const Scalar one  = KokkosKernels::ArithTraits<Scalar>::one();
  int return_value  = 0;

  Kokkos::initialize();
  {
    Matrix A;
    if (filename.empty()) {
      const Ordinal numRows      = 5000;
      Offset nnz                 = numRows * 10;
      const Ordinal bandwidth    = static_cast<Ordinal>(0.05 * numRows);
      const Scalar diagDominance = 2 * one;
      std::cout << "No filename provided; generating a " << numRows << "x" << numRows
                << " diagonally dominant sparse matrix.\n";
      A = KokkosSparse::Impl::kk_generate_diagonally_dominant_sparse_matrix<Matrix>(numRows, numRows, nnz, 0, bandwidth,
                                                                                    diagDominance);
    } else {
      // Read the sparse matrix A from a MatrixMarket file.
      std::cout << "Reading matrix from file: " << filename << "\n";
      A = KokkosSparse::Impl::read_kokkos_crst_matrix<Matrix>(filename.c_str());
    }

    // spiluk requires sorted column indices within each row.
    KokkosSparse::sort_crs_matrix(A);

    const Ordinal numRows = A.numRows();

    const Offset nnz = A.nnz();
    std::cout << "Matrix: " << numRows << " x " << A.numCols() << ", nnz=" << nnz << "\n";

    // Conservative upper bound for nnz in L and U.
    constexpr Offset expand_fact = 5;
    const Offset nnzLU           = expand_fact * nnz * (fillLevel + 1);

    // SEQLVLSCHD_TP1 uses a team-parallel level-scheduling algorithm.
    Handle kh;
    // kh.create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_RP, numRows, nnzLU, nnzLU);
    kh.create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1, numRows, nnzLU, nnzLU);

    // Allocate output row maps and initial entry arrays for L and U.
    RowMap L_row_map("L_row_map", numRows + 1);
    RowMap U_row_map("U_row_map", numRows + 1);
    Entries L_entries("L_entries", kh.get_spiluk_handle()->get_nnzL());
    Entries U_entries("U_entries", kh.get_spiluk_handle()->get_nnzU());

    // Determines the sparsity pattern (row maps and entries) of L and U.
    KokkosSparse::spiluk_symbolic(&kh, fillLevel, A.graph.row_map, A.graph.entries, L_row_map, L_entries, U_row_map,
                                  U_entries);

    Kokkos::fence();

    // Resize entry and value arrays to the exact sizes found by symbolic.
    Kokkos::resize(L_entries, kh.get_spiluk_handle()->get_nnzL());
    Kokkos::resize(U_entries, kh.get_spiluk_handle()->get_nnzU());
    Values L_values("L_values", kh.get_spiluk_handle()->get_nnzL());
    Values U_values("U_values", kh.get_spiluk_handle()->get_nnzU());

    // Computes the values of L and U given the sparsity pattern from symbolic.
    KokkosSparse::spiluk_numeric(&kh, fillLevel, A.graph.row_map, A.graph.entries, A.values, L_row_map, L_entries,
                                 L_values, U_row_map, U_entries, U_values);

    Kokkos::fence();

    // Save nonzero counts before destroying the handle.
    const Offset nnzL = kh.get_spiluk_handle()->get_nnzL();
    const Offset nnzU = kh.get_spiluk_handle()->get_nnzU();

    kh.destroy_spiluk_handle();

    // Verify the factorization: compute ||L*U*e - A*e|| / ||A*e||
    // where e is the all-ones vector.
    const Matrix L("L", numRows, numRows, nnzL, L_values, L_row_map, L_entries);
    const Matrix U("U", numRows, numRows, nnzU, U_values, U_row_map, U_entries);

    Values e("e", numRows);
    Values Ae("Ae", numRows);
    Values Ue("Ue", numRows);
    Kokkos::deep_copy(e, one);

    // Ae = A*e
    KokkosSparse::spmv("N", one, A, e, zero, Ae);
    using mag_t   = typename KokkosKernels::ArithTraits<Scalar>::mag_type;
    mag_t Ae_norm = KokkosBlas::nrm2(Ae);

    // Ue = U*e, then Ae = L*(U*e) - A*e (residual)
    KokkosSparse::spmv("N", one, U, e, zero, Ue);
    KokkosSparse::spmv("N", one, L, Ue, -one, Ae);
    mag_t res_norm = KokkosBlas::nrm2(Ae);
    mag_t rel_res  = res_norm / Ae_norm;

    std::cout << "spiluk (fill level=" << fillLevel << "): "
              << "nnzL=" << nnzL << ", nnzU=" << nnzU << "\n";
    std::cout << "  ||L*U*e - A*e|| / ||A*e|| = " << rel_res << "\n";

    const bool success = rel_res < 1e-3;
    if (success) {
      std::cout << "  SUCCESS: factorization residual is within tolerance.\n";
    } else {
      std::cout << "  FAILURE: factorization residual exceeds tolerance.\n";
    }

    return_value = success ? 0 : 1;
  }

  Kokkos::finalize();
  return return_value;
}
