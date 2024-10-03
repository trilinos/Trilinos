// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Tpetra_Core.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "BelosSolverFactory.hpp"
#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include "TpetraCore_ETIHelperMacros.h"

// FIXME (mfh 21 Aug 2015) Temporary work-around for Bug 6392.
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)
namespace Belos {
namespace Details {
//extern void registerLinearSolverFactory ();

// FIXME (mfh 24 Aug 2015) The weak symbol stuff in the main
// registration function doesn't work.  I have to reference the
// Tpetra registration function for now instead.
namespace Tpetra {
  extern void registerLinearSolverFactory ();
} // namespace Tpetra

} // namespace Details
} // namespace Belos
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)

namespace { // (anonymous)

// FIXME (mfh 21 Aug 2015) Temporary work-around for Bug 6392.
#if ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)
TEUCHOS_STATIC_SETUP()
{
  // FIXME (mfh 24 Aug 2015) The weak symbol stuff in the main
  // registration function doesn't work.  I have to reference the
  // Tpetra registration function for now instead.

  // Belos::Details::registerLinearSolverFactory ();
  Belos::Details::Tpetra::registerLinearSolverFactory ();
}
#endif // ! defined(HAVE_TEUCHOS_DYNAMIC_LIBS)

// Create a very simple square test matrix.  We use the identity
// matrix here.  The point of this test is NOT to exercise the solver;
// it's just to check that Belos' LinearSolverFactory can create
// working solvers.  Belos has more rigorous tests for its solvers.
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

// Create a very simple square test linear system (matrix, right-hand
// side(s), and exact solution(s).  We use the identity matrix here.
// The point of this test is NOT to exercise the preconditioner; it's
// just to check that its LinearSolverFactory can create working
// preconditioners.  Belos has more rigorous tests for each of its
// preconditioners.
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
            const Teuchos::RCP<const Tpetra::CrsMatrix<SC, LO, GO, NT> >& A,
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
  out << "Test solver \"" << solverName << "\" from Belos package" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<Trilinos::Details::LinearSolver<MV, OP, mag_type> > solver;
  try {
    solver = Trilinos::Details::getLinearSolver<MV, OP, mag_type> ("Belos", solverName);
  } catch (std::exception& e) {
    out << "*** FAILED: getLinearSolver threw an exception: " << e.what () << endl;
    success = false;
    return;
  }
  TEST_EQUALITY_CONST( solver.is_null (), false );
  if (solver.is_null ()) {
    out << "*** FAILED to create solver \"" << solverName
        << "\" from Belos package" << endl;
    success = false;
    return;
  }

  out << "Set matrix" << endl;
  solver->setMatrix (A);

  out << "Compute symbolic and numeric factorization" << endl;
  solver->symbolic ();
  solver->numeric ();

  out << "Apply solver to \"solve\" AX=B for X.  Belos already has solver "
    "tests; the point is to check that solve() doesn't throw." << endl;
  X.putScalar (STS::zero ());
  solver->solve (X, B);
}

//
// The actual unit test.
//
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( LinearSolverFactory, Solve, SC, LO, GO, NT )
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TypeNameTraits;
  using std::endl;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NT> MAT;
  typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
  typedef Tpetra::Operator<SC,LO,GO,NT> OP;

  return;

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

  Belos::SolverFactory<SC, MV, OP> factory;
  // FIXME (mfh 23 Aug 2015) Not all Belos solvers can handle solves
  // with the identity matrix.  BiCGSTAB might need a bit of work, for
  // example.  I'm not so worried about that for now but we should go
  // back and revisit this at some point.
  //
  // Teuchos::Array<std::string> solverNames = factory.supportedSolverNames ();
  // const int numSolvers = static_cast<int> (solverNames.size ());
  const char* solverNames[2] = {"CG", "GMRES"};
  const int numSolvers = 2;

  int numSolversTested = 0;
  for (int k = 0; k < numSolvers; ++k) {
    const std::string solverName (solverNames[k]);

    // Use Belos' factory to tell us whether the factory supports the
    // given combination of template parameters.  If create() throws,
    // the combination is not supported.
    bool skip = false;
    try {
      (void) factory.create (solverName, Teuchos::null);
    }
    catch (...) {
      skip = true;
    }
    if (! skip) {
      testSolver<SC, LO, GO, NT> (out, success, *X, A, *B, *X_exact, solverName);
      ++numSolversTested;
    }
  }

  out << "Tested " << numSolversTested << " solver(s) of " << numSolvers
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
TPETRA_ETI_MANGLING_TYPEDEFS()

// Macro that instantiates the unit test
#define LCLINST( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( LinearSolverFactory, Solve, SC, LO, GO, NT )

// Tpetra's ETI will instantiate the unit test for all enabled type
// combinations.
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCLINST )

} // namespace (anonymous)
