// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include "Kokkos_Core.hpp"

#include "KokkosKernels_default_types.hpp"
#include "KokkosSparse_spadd.hpp"

#include "KokkosKernels_Test_Structured_Matrix.hpp"

using Scalar  = KokkosKernels::default_scalar;
using Ordinal = KokkosKernels::default_lno_t;
using Offset  = KokkosKernels::default_size_type;
using Layout  = KokkosKernels::default_layout;

int main() {
  Kokkos::initialize();

  using device_type =
      typename Kokkos::Device<Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space>;
  using execution_space = typename device_type::execution_space;
  using memory_space    = typename device_type::memory_space;
  using matrix_type     = typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void, Offset>;

  int return_value = 0;

  {
    // The mat_structure view is used to generate a matrix using
    // finite difference (FD) or finite element (FE) discretization
    // on a cartesian grid.
    // Each row corresponds to an axis (x, y and z)
    // In each row the first entry is the number of grid point in
    // that direction, the second and third entries are used to apply
    // BCs in that direction.
    Kokkos::View<Ordinal* [3], Kokkos::HostSpace> mat_structure("Matrix Structure", 2);
    mat_structure(0, 0) = 10;  // Request 10 grid point in 'x' direction
    mat_structure(0, 1) = 1;   // Add BC to the left
    mat_structure(0, 2) = 1;   // Add BC to the right
    mat_structure(1, 0) = 10;  // Request 10 grid point in 'y' direction
    mat_structure(1, 1) = 1;   // Add BC to the bottom
    mat_structure(1, 2) = 1;   // Add BC to the top

    matrix_type A = Test::generate_structured_matrix2D<matrix_type>("FD", mat_structure);
    matrix_type B = Test::generate_structured_matrix2D<matrix_type>("FE", mat_structure);
    matrix_type C;

    // Create KokkosKernelHandle
    using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<Offset, Ordinal, Scalar, execution_space,
                                                                          memory_space, memory_space>;
    KernelHandle kh;
    kh.create_spadd_handle(false);

    const Scalar alpha = 2.5;
    const Scalar beta  = 1.2;

    KokkosSparse::spadd_symbolic(&kh, A, B, C);
    KokkosSparse::spadd_numeric(&kh, alpha, A, beta, B, C);
    kh.destroy_spadd_handle();

    std::cout << "spadd was performed correctly!" << std::endl;
  }

  Kokkos::finalize();

  return return_value;
}
