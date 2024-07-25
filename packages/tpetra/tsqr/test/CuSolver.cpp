// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Impl_CuBlas.hpp"
#include "Tsqr_Impl_CuSolver.hpp"
#include "Tsqr_Impl_CuTypes.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>
#include <type_traits>

namespace { // (anonymous)

template<class RealType>
void
verifyReal (std::ostream& out, bool& success)
{
  using TSQR::Impl::CuSolver;
  using TSQR::Impl::CuSolverHandle;
  using TSQR::Impl::getCuSolverHandleSingleton;
  using TSQR::Impl::CudaValue;
  using std::endl;

  std::shared_ptr<CuSolverHandle> s = getCuSolverHandleSingleton ();
  TEST_ASSERT( s.get () != nullptr );

  Kokkos::View<int, Kokkos::CudaSpace> info ("info");
  CuSolver<RealType> solver (s, info.data ());

  using IST = typename CudaValue<RealType>::type;
  static_assert (std::is_same<RealType, IST>::value,
                 "CudaValue::type is wrong.");
  const RealType x (666.0);
  out << "Original x: " << x << ": Converted x: "
      << CudaValue<RealType>::makeValue (x) << endl;

  using TSQR::Impl::CuBlasHandle;
  using TSQR::Impl::getCuBlasHandleSingleton;
  std::shared_ptr<CuBlasHandle> b = getCuBlasHandleSingleton ();
  TEST_ASSERT( b.get () != nullptr );

  using TSQR::Impl::CuBlas;
  CuBlas<RealType> blas1 (b);
  CuBlas<RealType> blas2;

  // Trivial test with 1x1 matrices, just to make sure that cuBLAS
  // actually works.
  using matrix_type =
    Kokkos::View<RealType**, Kokkos::LayoutLeft, Kokkos::CudaSpace>;
  matrix_type A ("A", 1, 1);
  matrix_type B ("B", 1, 1);
  matrix_type C ("C", 1, 1);
  const RealType alpha (1.0);
  const RealType beta (1.0);

  Kokkos::deep_copy (A, RealType (3.0));
  Kokkos::deep_copy (B, RealType (4.0));
  Kokkos::deep_copy (C, RealType (-5.0));
  const RealType C_expected_1 (7.0);
  blas1.gemm ('N', 'N', 1, 1, 1, alpha, A.data (), A.stride (1),
              B.data (), B.stride (1), beta, C.data (), C.stride (1));
  auto C_h = Kokkos::create_mirror_view (Kokkos::HostSpace (), C);
  Kokkos::deep_copy (C_h, C);
  TEST_EQUALITY_CONST( C_h(0,0), C_expected_1 );

  Kokkos::deep_copy (A, RealType (6.0));
  Kokkos::deep_copy (B, RealType (6.0));
  Kokkos::deep_copy (C, RealType (-5.0));
  const RealType C_expected_2 (31.0);
  blas2.gemm ('N', 'N', 1, 1, 1, alpha, A.data (), A.stride (1),
              B.data (), B.stride (1), beta, C.data (), C.stride (1));
  Kokkos::deep_copy (C_h, C);
  TEST_EQUALITY_CONST( C_h(0,0), C_expected_2 );
}

#ifdef HAVE_TPETRATSQR_COMPLEX
template<class ComplexType>
void
verifyComplex (std::ostream& out, bool& success)
{
  using TSQR::Impl::CuSolver;
  using TSQR::Impl::CuSolverHandle;
  using TSQR::Impl::getCuSolverHandleSingleton;
  using TSQR::Impl::CudaValue;
  using std::endl;

  std::shared_ptr<CuSolverHandle> s = getCuSolverHandleSingleton ();
  TEST_ASSERT( s.get () != nullptr );

  Kokkos::View<int, Kokkos::CudaSpace> info ("info");
  CuSolver<ComplexType> solver (s, info.data ());

  using IST = typename CudaValue<ComplexType>::type;

  using expected_z_IST = cuDoubleComplex;
  using expected_c_IST = cuFloatComplex;
  constexpr bool is_z =
    std::is_same<ComplexType, std::complex<double>>::value;
  using expected_IST = typename std::conditional<
    is_z,
    expected_z_IST,
    expected_c_IST>::type;
  static_assert (std::is_same<expected_IST, IST>::value,
                 "CudaValue::type is wrong.");
  const ComplexType x (666.0, 418.0);
  const IST x_out = CudaValue<ComplexType>::makeValue (x);
  out << "Original x: " << x << ": Converted x: ("
      << x_out.x << "," << x_out.y << ")" << endl;

  using TSQR::Impl::CuBlas;
  using TSQR::Impl::CuBlasHandle;
  using TSQR::Impl::getCuBlasHandleSingleton;
  std::shared_ptr<CuBlasHandle> b = getCuBlasHandleSingleton ();
  TEST_ASSERT( b.get () != nullptr );

  CuBlas<ComplexType> blas (b);
}
#endif // HAVE_TPETRATSQR_COMPLEX

void
verify (std::ostream& out, bool& success)
{
  verifyReal<double> (out, success);
  verifyReal<float> (out, success);

#ifdef HAVE_TPETRATSQR_COMPLEX
  verifyComplex<std::complex<double>> (out, success);
  verifyComplex<std::complex<float>> (out, success);
#endif // HAVE_TPETRATSQR_COMPLEX
}

} // namespace (anonymous)

int
main (int argc, char *argv[])
{
  using std::cout;
  using std::endl;

  cout << "Test cuBLAS and cuSOLVER handle creation" << endl;

  bool success = true;
  try {
    Kokkos::ScopeGuard kokkosScope (argc, argv);
    verify (cout, success);
    // The Trilinos test framework expects a message like this.
    if (success) {
      cout << "\nEnd Result: TEST PASSED" << endl;
    }
    else {
      cout << "\nEnd Result: TEST FAILED" << endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
