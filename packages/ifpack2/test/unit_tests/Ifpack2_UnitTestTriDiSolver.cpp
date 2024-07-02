// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_UnitTestDenseSolver.cpp
// \brief Unit test for Ifpack2::Details::TriDiSolver.

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Details_TriDiSolver.hpp>
#include <Ifpack2_Details_DenseSolver.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialTriDiMatrix.hpp>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(TriDiSolver, LapackComparison, ScalarType, LocalOrdinalType, GlobalOrdinalType)
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
  //  typedef Teuchos::SerialTriDiMatrix<local_ordinal_type, scalar_type> crs_matrix_type;
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

  typedef Ifpack2::Details::TriDiSolver<row_matrix_type> solver_type;
    //  typedef Ifpack2::Details::DenseSolver<row_matrix_type> solver_type;

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

  RCP<const map_type> rowMap = rcp (new map_type (globalNumRows, localNumRows, indexBase, comm));
  RCP<const map_type> colMap = rowMap;
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Creating matrix A" << endl;

  // The matrix will be block diagonal, so we can use the row Map as
  // the column Map.
  RCP<crs_matrix_type> A = rcp(new crs_matrix_type (rowMap, colMap, 3));

  Array<scalar_type> val (3);
  Array<local_ordinal_type> ind (3);

  out << "Filling A_tridi" << endl;

  Teuchos::SerialTriDiMatrix<int, scalar_type> A_tridi (3, 3);
  A_tridi(0, 0) = as<scalar_type> (2);
  A_tridi(0, 1) = as<scalar_type> (-1);
  //  A_tridi(0, 2) = as<scalar_type> (0);
  A_tridi(1, 0) = as<scalar_type> (-1);
  A_tridi(1, 1) = as<scalar_type> (3);
  A_tridi(1, 2) = as<scalar_type> (-1);
  // A_tridi(2, 0) = as<scalar_type> (7);
  A_tridi(2, 1) = as<scalar_type> (-2);
  A_tridi(2, 2) = as<scalar_type> (1); //

  A_tridi.print(out);

  out << "Filling A CRS" << endl;

  val[0] = A_tridi(0, 0);
  val[1] = A_tridi(0, 1);
  //  val[2] = A_tridi(0, 2);
  val[2] = 0;
  ind[0] = 0;
  ind[1] = 1;
  ind[2] = 2;
  A->insertLocalValues (0, ind (), val ());

  val[0] = A_tridi(1, 0);
  val[1] = A_tridi(1, 1);
  val[2] = A_tridi(1, 2);
  ind[0] = 0;
  ind[1] = 1;
  ind[2] = 2;
  A->insertLocalValues (1, ind (), val ());

  //  val[0] = A_tridi(2, 0);
  val[0] = 0;
  val[1] = A_tridi(2, 1);
  val[2] = A_tridi(2, 2);
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
    "Dense failed to solve the problem to within a small tolerance "
    << 10*STS::eps () << ".  Relative residual norm: " << relResNorm << ".");

  // Permutation array for LAPACK's LU factorization.
  Array<int> ipiv (A_tridi.numRowsCols ());
  // Fill the LU permutation array with zeros.
  std::fill (ipiv.begin (), ipiv.end (), 0);
  // Compute the LU factorization.
  Teuchos::LAPACK<int, scalar_type> lapack;
  int INFO = 0;

  out << "A_tridi before GETRF:" << endl;
  A_tridi.print(out);


  out << "Calling GTTRF (" << A_tridi.numRowsCols () << ", A.DL, A.D, A.DU, A.DU2, IPIV, INFO)" << endl;

  // Remember that LAPACK's LU factorization overwrites its input.
  lapack.GTTRF (A_tridi.numRowsCols(),
                A_tridi.DL(),
                A_tridi.D(),
                A_tridi.DU(),
                A_tridi.DU2(),
                ipiv.getRawPtr (), &INFO);

  out << "A_tridi after GTTRF:" << endl;

  A_tridi.print(out);
  out<<std::endl;

  out << "ipiv after GETRF: " << toString (ipiv) << endl
      << "INFO: " << INFO << endl;

  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO < 0, std::logic_error, "Bug in test after calling GTTRF: "
    "INFO = " << INFO << " < 0.  "
    "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO > 0, std::logic_error, "Bug in test after calling GTTRF: "
    "INFO = " << INFO << " > 0.  "
    "This means that the U factor of the test matrix is exactly singular, "
    "and therefore that the test's authors need to come up with a different "
    "test matrix.  Please report this bug to the Ifpack2 developers.");

  Array<scalar_type> x_lapack (A_tridi.numRowsCols());
  // LAPACK overwrites its input.
  ArrayRCP<const scalar_type> b_view = b.get1dView ();
  std::copy (b_view.begin (), b_view.end (), x_lapack.begin ());

  const int numRhs = 1;
  out << "Calling GTTRS ('N', " << A_tridi.numRowsCols() << ", " << numRhs
      << ", A.DL, A.D, A.DU, A.DU2, IPIV, X, "
      << static_cast<int> (x_lapack.size ()) << ", INFO)" << endl;
  lapack.GTTRS ('N', A_tridi.numRowsCols (), numRhs,
                A_tridi.DL(),
                A_tridi.D(),
                A_tridi.DU(),
                A_tridi.DU2(),
                ipiv.getRawPtr (), x_lapack.getRawPtr (),
                static_cast<int> (x_lapack.size ()), &INFO);
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO != 0, std::logic_error, "Bug in test after calling GTTRS: "
    "INFO = " << INFO << " != 0.  "
    "This means that the U factor of the test matrix is exactly singular, "
    "and therefore that the test's authors need to come up with a different "
    "test matrix.  Please report this bug to the Ifpack2 developers.");
  out << "Exact solution: " << toString (x_lapack) << endl;

  // Compare LAPACK solution against TriDiSolver solution.  They
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
    "TriDiSolver failed to reproduce LAPACK's solution to within a small tolerance "
    << 10*STS::eps () << ".  \\|x_lapack - x_computed\\|_2 / \\|x_lapack\\|_2 = "
    << lapackRelResNorm << ".");
}

// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( TriDiSolver, LapackComparison, SC, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

