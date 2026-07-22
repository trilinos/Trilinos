// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include "Kokkos_Core.hpp"

#include "KokkosKernels_default_types.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"

#include "KokkosKernels_Test_Structured_Matrix.hpp"

using Scalar  = KokkosKernels::default_scalar;
using Ordinal = KokkosKernels::default_lno_t;
using Offset  = KokkosKernels::default_size_type;
using Layout  = KokkosKernels::default_layout;

template <class Yvector>
struct check_spmv_functor {
  Yvector y;
  const Scalar SC_ONE = KokkosKernels::ArithTraits<Scalar>::one();

  check_spmv_functor(Yvector y_) : y(y_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, Ordinal& update) const {
    if (y(i) != (SC_ONE + SC_ONE)) {
      ++update;
    }
  }
};

int main() {
  Kokkos::initialize();

  using device_type =
      typename Kokkos::Device<Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space>;
  using matrix_type = typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void, Offset>;
  using values_type = typename matrix_type::values_type;

  int return_value = 0;

  {
    const Scalar SC_ONE = KokkosKernels::ArithTraits<Scalar>::one();

    // The mat_structure view is used to generate a matrix using
    // finite difference (FD) or finite element (FE) discretization
    // on a cartesian grid.
    // Each row corresponds to an axis (x, y and z)
    // In each row the first entry is the number of grid point in
    // that direction, the second and third entries are used to apply
    // BCs in that direction, BC=0 means Neumann BC is applied,
    // BC=1 means Dirichlet BC is applied by zeroing out the row and putting
    // one on the diagonal.
    Kokkos::View<Ordinal* [3], Kokkos::HostSpace> mat_structure("Matrix Structure", 2);
    mat_structure(0, 0) = 10;  // Request 10 grid point in 'x' direction
    mat_structure(0, 1) = 0;   // Add BC to the left
    mat_structure(0, 2) = 0;   // Add BC to the right
    mat_structure(1, 0) = 10;  // Request 10 grid point in 'y' direction
    mat_structure(1, 1) = 0;   // Add BC to the bottom
    mat_structure(1, 2) = 0;   // Add BC to the top

    matrix_type myMatrix = Test::generate_structured_matrix2D<matrix_type>("FD", mat_structure);

    const Ordinal numRows = myMatrix.numRows();

    const Scalar alpha = SC_ONE;
    const Scalar beta  = SC_ONE;

    typename values_type::non_const_type x("lhs", numRows);
    typename values_type::non_const_type y("rhs", numRows);
    Kokkos::deep_copy(x, SC_ONE);
    Kokkos::deep_copy(y, SC_ONE + SC_ONE);

    KokkosSparse::spmv("N", alpha, myMatrix, x, beta, y);

    Ordinal count_errors = 0;
    check_spmv_functor<values_type> check_spmv(y);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Ordinal>(0, numRows), check_spmv, count_errors);
    if (count_errors > 0) {
      return_value = 1;
      std::cout << "Found " << count_errors << " errors in y vector!" << std::endl;
    } else {
      std::cout << "spmv was performed correctly: y = beta*y + alpha*A*x" << std::endl;
    }
  }

  Kokkos::finalize();

  return return_value;
}
