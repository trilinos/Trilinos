// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Ifpack2_Factory.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
// Define typedefs and macros for testing over all template parameter
// combinations.
#include "Ifpack2_ETIHelperMacros.h"
#include <type_traits>

// FIXME (mfh 21 Aug 2015) Temporary work-around for Bug 6392.
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)
namespace Ifpack2 {
namespace Details {
  // FIXME (mfh 21 Aug 2015) NONE of the commented-out things work.
  //
  // extern void __attribute__((weak)) registerLinearSolverFactory ();
  // void __attribute__((weak)) registerLinearSolverFactory ();
  // #pragma weak registerLinearSolverLibrary

  extern void registerLinearSolverFactory ();

} // namespace Details
} // namespace Ifpack2
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)

namespace {
// FIXME (mfh 21 Aug 2015) Temporary work-around for Bug 6392.
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)
  TEUCHOS_STATIC_SETUP()
  {
    // if (Ifpack2::Details::registerLinearSolverFactory == NULL) {
    //   std::cout << "-- Ifpack2::Details::registerLinearSolverFactory is NULL" << std::endl;
    // } else {
    //   Ifpack2::Details::registerLinearSolverFactory ();
    // }

    Ifpack2::Details::registerLinearSolverFactory ();
  }
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)

  // Create a very simple square test matrix.  We use the identity
  // matrix here.  The point of this test is NOT to exercise the
  // preconditioner; it's just to check that its LinearSolverFactory
  // can create working preconditioners.  Ifpack2 has more rigorous
  // tests for each of its preconditioners.
  template<class SC, class LO, class GO, class NT>
  Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >
  createTestMatrix (Teuchos::FancyOStream& out,
                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                    const Tpetra::global_size_t gblNumRows)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NT> MAT;
    typedef Tpetra::Map<LO,GO,NT> map_type;
    typedef Teuchos::ScalarTraits<SC> STS;

    Teuchos::OSTab tab0 (out);
    out << "Create test matrix with " << gblNumRows << " row(s)" << endl;
    Teuchos::OSTab tab1 (out);

    TEUCHOS_TEST_FOR_EXCEPTION
      ( gblNumRows == 0, std::invalid_argument, "gblNumRows = 0");

    const GO indexBase = 0;
    RCP<const map_type> rowMap (new map_type (gblNumRows, indexBase, comm));
    // For this particular row matrix, the row Map and the column Map
    // are the same.  Giving a column Map to CrsMatrix's constructor
    // lets us use local indices.
    RCP<const map_type> colMap = rowMap;
    const size_t maxNumEntPerRow = 1;
    RCP<MAT> A (new MAT (rowMap, colMap, maxNumEntPerRow));

    if (rowMap->getLocalNumElements () != 0) {
      Teuchos::Array<SC> vals (1);
      Teuchos::Array<LO> inds (1);
      for (LO lclRow = rowMap->getMinLocalIndex ();
           lclRow <= rowMap->getMaxLocalIndex (); ++lclRow) {
        inds[0] = lclRow;
        vals[0] = STS::one ();
        A->insertLocalValues (lclRow, inds (), vals ());
      }
    }

    RCP<const map_type> domMap = rowMap;
    RCP<const map_type> ranMap = rowMap;
    A->fillComplete (domMap, ranMap);
    return A;
  }

  // Create a very simple square test linear system (matrix,
  // right-hand side(s), and exact solution(s).  We use the identity
  // matrix here.  The point of this test is NOT to exercise the
  // preconditioner; it's just to check that its LinearSolverFactory
  // can create working preconditioners.  Ifpack2 has more rigorous
  // tests for each of its preconditioners.
  template<class SC, class LO, class GO, class NT>
  void
  createTestProblem (Teuchos::FancyOStream& out,
                     Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NT> >& X,
                     Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >& A,
                     Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NT> >& B,
                     const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const Tpetra::global_size_t gblNumRows,
                     const size_t numVecs)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
    typedef Teuchos::ScalarTraits<SC> STS;

    A = createTestMatrix<SC, LO, GO, NT> (out, comm, gblNumRows);
    X = rcp (new MV (A->getDomainMap (), numVecs));
    B = rcp (new MV (A->getRangeMap (), numVecs));

    X->putScalar (STS::one ());
    A->apply (*X, *B);
  }

  template<class SC, class LO, class GO, class NT>
  void
  testSolver (Teuchos::FancyOStream& out,
              bool& success,
              Tpetra::MultiVector<SC, LO, GO, NT>& X,
              const Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >& A,
              const Tpetra::MultiVector<SC, LO, GO, NT>& B,
              const Tpetra::MultiVector<SC, LO, GO, NT>& X_exact,
              const std::string& solverName)
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    typedef Tpetra::Operator<SC,LO,GO,NT> OP;
    typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
    typedef Teuchos::ScalarTraits<SC> STS;
    typedef typename MV::mag_type mag_type;

    Teuchos::OSTab tab0 (out);
    out << "Test solver \"" << solverName << "\" from Ifpack2 package" << endl;
    Teuchos::OSTab tab1 (out);

    RCP<Trilinos::Details::LinearSolver<MV, OP, mag_type> > solver;
    try {
      solver = Trilinos::Details::getLinearSolver<MV, OP, mag_type> ("Ifpack2", solverName);
    } catch (std::exception& e) {
      out << "*** FAILED: getLinearSolver threw an exception: " << e.what () << endl;
      success = false;
      return;
    }
    TEST_EQUALITY_CONST( solver.is_null (), false );
    if (solver.is_null ()) {
      out << "*** FAILED to create solver \"" << solverName
          << "\" from Ifpack2 package" << endl;
      success = false;
      return;
    }

    out << "Set matrix" << endl;
    // Set all matrix entries to nan.
    A->setAllToScalar(STS::nan());
    solver->setMatrix (A);

    out << "Compute symbolic factorization" << endl;
    solver->symbolic ();
    // Set all matrix entries to one, so that A is the identity matrix.
    A->setAllToScalar(STS::one());
    out << "Compute numeric factorization" << endl;
    solver->numeric ();

    out << "Apply solver to \"solve\" AX=B for X.  Ifpack2 only promises "
      "preconditioning, so we don't expect to get X right.  Mainly we check "
      "that solve() doesn't throw." << endl;
    X.putScalar (STS::zero ());
    solver->solve (X, B);

    // Set all matrix entries to 2, so that A is a scaled identity matrix.
    RCP<MV> X2 = rcp (new MV (X.getMap (), X.getNumVectors()));
    A->setAllToScalar(STS::one() + STS::one());
    out << "Recompute numeric factorization" << endl;
    solver->numeric ();

    out << "Apply solver to \"solve\" (2*A)X=B for X.  Ifpack2 only promises "
      "preconditioning, so we don't expect to get X right.  Mainly we check "
      "that solve() doesn't throw." << endl;
    X2->putScalar (STS::zero ());
    solver->solve (*X2, B);

    // We aasume that the preconditioner is linear in A and compare the two resulting vectors.
    Kokkos::View<typename MV::mag_type *, Kokkos::HostSpace> normsX   ("normsX", X.getNumVectors());
    Kokkos::View<typename MV::mag_type *, Kokkos::HostSpace> normsX2  ("normsX2", X.getNumVectors());
    Kokkos::View<typename MV::mag_type *, Kokkos::HostSpace> normsDiff("normsDiff", X.getNumVectors());
    X.norm2(normsX);
    X2->norm2(normsX2);
    // X2 := 2*X2 - X
    X2->update(-STS::one(), X, STS::one()+STS::one());
    X2->norm2(normsDiff);
    for (size_t k = 0; k<X.getNumVectors(); k++) {
      out << "normsX[" << k << "] = " << normsX(k) << ", normsX2[" << k << "] = " << normsX2(k) << ", normsDiff[" << k << "] = " << normsDiff(k) << std::endl;
      TEUCHOS_ASSERT_INEQUALITY(normsDiff(k), <, 100 * STS::eps());
    }
  }

  //
  // The actual unit test.
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SolverFactory, Solve, SC, LO, GO, NT )
  {
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::TypeNameTraits;
    using std::endl;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NT> MAT;
    typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
    typedef Tpetra::RowMatrix<SC,LO,GO,NT> row_matrix_type;

#if ! defined(TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION)
    out << "LinearSolverFactory run-time registration disabled; "
      "not running test" << endl;
    return;
#endif // NOT TRILINOS_HAVE_LINEAR_SOLVER_FACTORY_REGISTRATION

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    const Tpetra::global_size_t gblNumRows = comm->getSize () * 10;
    const size_t numVecs = 3;

    out << "Create test problem" << endl;
    RCP<MV> X_exact, B;
    RCP<MAT> A;
    createTestProblem (out, X_exact, A, B, comm, gblNumRows, numVecs);

    RCP<MV> X = rcp (new MV (X_exact->getMap (), numVecs));

    // FIXME (mfh 26 Jul 2015) Need to test more solvers.  In
    // particular, it's important to test AdditiveSchwarz.
    const std::vector<std::string> solverNames = {
      "CHEBYSHEV",
      "DENSE",
      "AMESOS2"
      "DIAGONAL",
      "ILUT",
      "RELAXATION",
      "RILUK",
      "MDF",
      "RBILUK",
      "FAST_IC",
      "FAST_ILU",
      "FAST_ILU_B",
      "FAST_ILDL",
      "BLOCK RELAXATION",
      // "DATABASE SCHWARZ",  // Skipping because it fails
      "SPARSE_BLOCK_RELAXATION"
      "TRIDI_RELAXATION",
      "BANDED_RELAXATION",
      // "IDENTITY", // Skipping because it does not depend at all on A.
      "LOCAL SPARSE TRIANGULAR SOLVER",
      // "HIPTMAIR", // Skipping because it needs auxiliary matrix
      "HYPRE",
    };
    int numSolversTested = 0;
    for (size_t k = 0; k < solverNames.size(); ++k) {
      const std::string solverName (solverNames[k]);

      // Use Ifpack2's factory to tell us whether the factory supports
      // the given combination of template parameters.  If create()
      // throws, the combination is not supported.
      bool skip = false;
      try {
        (void) Ifpack2::Factory::create<row_matrix_type> (solverName, A);
      }
      catch (...) {
        skip = true;
      }
      skip |= (solverName == "MDF" && !std::is_arithmetic_v<SC>); // skip complex types in mdf for now
      if (! skip) {
        testSolver<SC, LO, GO, NT> (out, success, *X, A, *B, *X_exact, solverName);
        ++numSolversTested;
      }
    }

    out << "Tested " << numSolversTested << " solver(s) of " << solverNames.size()
        << endl;
    if (numSolversTested == 0) {
      out << "*** ERROR: Tested no solvers for template parameters"
          << "SC = " << TypeNameTraits<SC>::name ()
          << ", LO = " << TypeNameTraits<LO>::name ()
          << ", GO = " << TypeNameTraits<GO>::name ()
          << ", NT = " << TypeNameTraits<NT>::name () << endl;
      success = false;
    }
  }

  // Define typedefs that make the Tpetra macros work.
  IFPACK2_ETI_MANGLING_TYPEDEFS()

// Macro that instantiates the unit test
#define LCLINST( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( SolverFactory, Solve, SC, LO, GO, NT )

// Ifpack2's ETI will instantiate the unit test for all enabled type
// combinations.
IFPACK2_INSTANTIATE_SLGN( LCLINST )

} // namespace (anonymous)
