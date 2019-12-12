//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER

#include "Tsqr_Impl_CuBlasHandle.hpp"
#include "Tsqr_Impl_CuSolverHandle.hpp"
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
  using TSQR::Impl::CudaValue;
  using std::endl;

  CuSolverHandle s = CuSolverHandle::getSingleton ();
  TEST_ASSERT( s.getHandle () != nullptr );

  Kokkos::View<int, Kokkos::CudaSpace> info ("info");
  CuSolver<RealType> solver (s, info.data ());

  using IST = typename CudaValue<RealType>::type;
  static_assert (std::is_same<RealType, IST>::value,
                 "CudaValue::type is wrong.");
  const RealType x (666.0);
  out << "Original x: " << x << ": Converted x: "
      << CudaValue<RealType>::makeValue (x) << endl;

  using TSQR::Impl::CuBlas;
  using TSQR::Impl::CuBlasHandle;
  CuBlasHandle b = CuBlasHandle::getSingleton ();
  TEST_ASSERT( b.getHandle () != nullptr );

  CuBlas<RealType> blas (b);
}

#ifdef HAVE_TPETRATSQR_COMPLEX
template<class ComplexType>
void
verifyComplex (std::ostream& out, bool& success)
{
  using TSQR::Impl::CuSolver;
  using TSQR::Impl::CuSolverHandle;
  using TSQR::Impl::CudaValue;
  using std::endl;

  CuSolverHandle s = CuSolverHandle::getSingleton ();
  TEST_ASSERT( s.getHandle () != nullptr );

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
  CuBlasHandle b = CuBlasHandle::getSingleton ();
  TEST_ASSERT( b.getHandle () != nullptr );

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
