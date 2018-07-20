/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/// \file Ifpack2_UnitTestDenseSolver.cpp
// \brief Unit test for Ifpack2::Details::DenseSolver.

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Details_DenseSolver.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <Teuchos_LAPACK.hpp>


namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(DenseSolver, LapackComparison, ScalarType, LocalOrdinalType, GlobalOrdinalType)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::toString;
  using std::endl;

  typedef ScalarType scalar_type;
  typedef LocalOrdinalType local_ordinal_type;
  typedef GlobalOrdinalType global_ordinal_type;
  typedef Tpetra::Map<>::node_type node_type;

  typedef Tpetra::global_size_t GST;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Tpetra::CrsMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> crs_matrix_type;
  typedef Tpetra::RowMatrix<scalar_type,
                            local_ordinal_type,
                            global_ordinal_type,
                            node_type> row_matrix_type;
  typedef Tpetra::Vector<scalar_type,
                         local_ordinal_type,
                         global_ordinal_type,
                         node_type> vec_type;
  typedef Tpetra::Map<local_ordinal_type,
                      global_ordinal_type,
                      node_type> map_type;
  typedef Ifpack2::Details::DenseSolver<row_matrix_type> solver_type;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  // We are now in a class method declared by the above macro.
  // The method has these input arguments:
  // (Teuchos::FancyOStream& out, bool& success)

  out << "Ifpack2::Version(): " << Ifpack2::Version () << endl;

  out << "Creating Maps" << endl;

  // Create a square matrix with 3 rows per process.
  // Give it diagonal blocks that are easy for LAPACK to solve.

  const size_t localNumRows = 3;
  const GST globalNumRows = comm->getSize () * localNumRows;
  const global_ordinal_type indexBase = 0;

  RCP<const map_type> rowMap (new map_type (globalNumRows, localNumRows, indexBase, comm));
  RCP<const map_type> colMap = rowMap;
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Creating matrix A" << endl;

  // The matrix will be block diagonal, so we can use the row Map as
  // the column Map.
  RCP<crs_matrix_type> A (new crs_matrix_type (rowMap, colMap, (size_t) 0));

  Array<scalar_type> val (3);
  Array<local_ordinal_type> ind (3);

  out << "Filling A_dense" << endl;

  // [1 2 3;
  //  4 5 6;
  //  7 8 10]
  Teuchos::SerialDenseMatrix<int, scalar_type> A_dense (3, 3);
  A_dense(0, 0) = as<scalar_type> (1);
  A_dense(0, 1) = as<scalar_type> (2);
  A_dense(0, 2) = as<scalar_type> (3);
  A_dense(1, 0) = as<scalar_type> (4);
  A_dense(1, 1) = as<scalar_type> (5);
  A_dense(1, 2) = as<scalar_type> (6);
  A_dense(2, 0) = as<scalar_type> (7);
  A_dense(2, 1) = as<scalar_type> (8);
  A_dense(2, 2) = as<scalar_type> (10); // not 9, else singular

  // A_dense(0, 0) = as<scalar_type> (10); // not 9, else singular
  // A_dense(0, 1) = as<scalar_type> (8);
  // A_dense(0, 2) = as<scalar_type> (7);
  // A_dense(1, 0) = as<scalar_type> (6);
  // A_dense(1, 1) = as<scalar_type> (5);
  // A_dense(1, 2) = as<scalar_type> (4);
  // A_dense(2, 0) = as<scalar_type> (3);
  // A_dense(2, 1) = as<scalar_type> (2);
  // A_dense(2, 2) = as<scalar_type> (1);

  out << "Filling A (sparse)" << endl;

  val[0] = A_dense(0, 0);
  val[1] = A_dense(0, 1);
  val[2] = A_dense(0, 2);
  ind[0] = 0;
  ind[1] = 1;
  ind[2] = 2;
  A->insertLocalValues (0, ind (), val ());

  val[0] = A_dense(1, 0);
  val[1] = A_dense(1, 1);
  val[2] = A_dense(1, 2);
  ind[0] = 0;
  ind[1] = 1;
  ind[2] = 2;
  A->insertLocalValues (1, ind (), val ());

  val[0] = A_dense(2, 0);
  val[1] = A_dense(2, 1);
  val[2] = A_dense(2, 2);
  ind[0] = 0;
  ind[1] = 1;
  ind[2] = 2;
  A->insertLocalValues (2, ind (), val ());

  out << "Calling fillComplete on A" << endl;

  A->fillComplete (domMap, ranMap);

  out << "Creating x_exact and b for the test problem" << endl;
  // Randomizing makes results not reproducible across platforms.
  vec_type x_exact (domMap);
  {
    ArrayRCP<scalar_type> x_exact_view = x_exact.get1dViewNonConst ();
    x_exact_view[0] = 1;
    x_exact_view[1] = 4;
    x_exact_view[2] = 9;
  }
  out << "x_exact: " << toString ((x_exact.get1dView ()) ()) << endl;

  vec_type b (ranMap);
  A->apply (x_exact, b); // b := A*x_exact
  out << "b: " << toString ((b.get1dView ()) ()) << endl;

  out << "Creating solver" << endl;
  solver_type solver (A);

  out << "Calling initialize" << endl;
  solver.initialize ();

  out << "Calling compute" << endl;
  solver.compute ();

  out << "Calling solver's describe() method" << endl;
  solver.describe (out, Teuchos::VERB_EXTREME);

  out << "Calling apply" << endl;
  vec_type x_computed (ranMap);
  solver.apply (b, x_computed); // solve A*x_computed=b for x_computed

  out << "x_computed: " << toString ((x_computed.get1dView ()) ()) << endl;

  // Compute the residual.
  vec_type r (ranMap);
  A->apply (x_computed, r); // r := A*x_computed
  r.update (STS::one (), b, -STS::one ()); // r := b - A*x_computed
  const magnitude_type absResNorm = r.norm2 ();
  out << "\\| b - A*x_computed \\|_2 = " << absResNorm << endl;
  const magnitude_type normB = b.norm2 ();
  const magnitude_type relResNorm = absResNorm / normB;
  out << "\\| b - A*x_computed \\|_2 / \\|b\\|_2 = " << relResNorm << endl;

  out << "A->describe() result:" << endl;
  out.setOutputToRootOnly(-1);
  A->describe (out, Teuchos::VERB_EXTREME);

  out.setOutputToRootOnly(0);
  out << "x_exact.describe() result:" << endl;
  out.setOutputToRootOnly(-1);
  x_exact.describe (out, Teuchos::VERB_EXTREME);
  out.setOutputToRootOnly(0);

  out << "b.describe() result:" << endl;
  out.setOutputToRootOnly(-1);
  b.describe (out, Teuchos::VERB_EXTREME);
  out.setOutputToRootOnly(0);

  out << "x_computed.describe() result:" << endl;
  out.setOutputToRootOnly(-1);
  x_computed.describe (out, Teuchos::VERB_EXTREME);
  out.setOutputToRootOnly(0);

  out << "r.describe() result:" << endl;
  out.setOutputToRootOnly(-1);
  r.describe (out, Teuchos::VERB_EXTREME);
  out.setOutputToRootOnly(0);

  out << "solver.describe() result:" << endl;
  out.setOutputToRootOnly(-1);
  solver.describe (out, Teuchos::VERB_EXTREME);
  out.setOutputToRootOnly(0);

  TEUCHOS_TEST_FOR_EXCEPTION(
    relResNorm > 10*STS::eps (), std::logic_error,
    "DenseSolver failed to solve the problem to within a small tolerance "
    << 10*STS::eps () << ".  Relative residual norm: " << relResNorm << ".");

  // Permutation array for LAPACK's LU factorization.
  Array<int> ipiv (A_dense.numCols ());
  // Fill the LU permutation array with zeros.
  std::fill (ipiv.begin (), ipiv.end (), 0);
  // Compute the LU factorization.
  Teuchos::LAPACK<int, scalar_type> lapack;
  int INFO = 0;

  out << "A_dense before GETRF:" << endl;
  {
    OSTab tab1 (rcpFromRef (out));
    out << "[";
    for (int i = 0; i < A_dense.numRows (); ++i) {
      for (int j = 0; j < A_dense.numCols (); ++j) {
        out << A_dense(i,j);
        if (j + 1 < A_dense.numCols ()) {
          out << ", ";
        }
      }
      if (i + 1 < A_dense.numRows ()) {
        out << ";" << endl;
      }
    }
    out << "]" << endl;
  }

  out << "Calling GETRF (" << A_dense.numRows () << ", " << A_dense.numCols ()
      << ", A, " << A_dense.stride () << ", IPIV, INFO)" << endl;

  // Remember that LAPACK's LU factorization overwrites its input.
  lapack.GETRF (A_dense.numRows (), A_dense.numCols (),
                A_dense.values (), A_dense.stride (),
                ipiv.getRawPtr (), &INFO);

  out << "A_dense after GETRF:" << endl;
  {
    OSTab tab1 (rcpFromRef (out));
    out << "[";
    for (int i = 0; i < A_dense.numRows (); ++i) {
      for (int j = 0; j < A_dense.numCols (); ++j) {
        out << A_dense(i,j);
        if (j + 1 < A_dense.numCols ()) {
          out << ", ";
        }
      }
      if (i + 1 < A_dense.numRows ()) {
        out << ";" << endl;
      }
    }
    out << "]" << endl;
  }
  out << "ipiv after GETRF: " << toString (ipiv) << endl
      << "INFO: " << INFO << endl;

  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO < 0, std::logic_error, "Bug in test after calling GETRF: "
    "INFO = " << INFO << " < 0.  "
    "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO > 0, std::logic_error, "Bug in test after calling GETRF: "
    "INFO = " << INFO << " > 0.  "
    "This means that the U factor of the test matrix is exactly singular, "
    "and therefore that the test's authors need to come up with a different "
    "test matrix.  Please report this bug to the Ifpack2 developers.");

  Array<scalar_type> x_lapack (A_dense.numCols ());
  // LAPACK overwrites its input.
  ArrayRCP<const scalar_type> b_view = b.get1dView ();
  std::copy (b_view.begin (), b_view.end (), x_lapack.begin ());

  const int numRhs = 1;
  out << "Calling GETRS ('N', " << A_dense.numRows () << ", " << numRhs
      << ", A, " << A_dense.stride () << ", IPIV, X, "
      << static_cast<int> (x_lapack.size ()) << ", INFO)" << endl;
  lapack.GETRS ('N', A_dense.numRows (), numRhs,
                A_dense.values (), A_dense.stride (),
                ipiv.getRawPtr (), x_lapack.getRawPtr (),
                static_cast<int> (x_lapack.size ()), &INFO);
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO != 0, std::logic_error, "Bug in test after calling GETRS: "
    "INFO = " << INFO << " != 0.  "
    "This means that the U factor of the test matrix is exactly singular, "
    "and therefore that the test's authors need to come up with a different "
    "test matrix.  Please report this bug to the Ifpack2 developers.");
  out << "Exact solution: " << toString (x_lapack) << endl;

  // Compare LAPACK solution against DenseSolver solution.  They
  // should be identical, or nearly so (if LAPACK was nondeterministic
  // -- unlikely with such a small problem).

  // Make a copy of x_lapack as a MultiVector.  (This is NOT a view).
  vec_type x_lapack_mv (x_exact.getMap (), x_lapack ());

  // Measure the difference between x_lapack and x_computed.
  r = x_lapack_mv;
  r.update (STS::one (), x_computed, -STS::one ()); // r := x_computed - x_lapack
  const magnitude_type lapackAbsResNorm = r.norm2 ();
  const magnitude_type lapackRelResNorm = lapackAbsResNorm / x_lapack_mv.norm2 ();

  out << "\\|x_lapack - x_computed\\|_2 / \\|x_lapack\\|_2 = "
      << lapackRelResNorm << endl;

  TEUCHOS_TEST_FOR_EXCEPTION(
    lapackRelResNorm > 10*STS::eps (), std::logic_error,
    "DenseSolver failed to reproduce LAPACK's solution to within a small tolerance "
    << 10*STS::eps () << ".  \\|x_lapack - x_computed\\|_2 / \\|x_lapack\\|_2 = "
    << lapackRelResNorm << ".");
}

// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP_SC_LO_GO(Scalar,LocalOrdinal,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( DenseSolver, LapackComparison, Scalar, LocalOrdinal, GlobalOrdinal)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

