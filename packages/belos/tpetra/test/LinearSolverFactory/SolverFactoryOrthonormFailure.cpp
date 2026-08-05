// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_MultiVector.hpp"
#include "TpetraCore_ETIHelperMacros.h"
#include <vector>
#include <string>

namespace { // (anonymous)

// Create a very simple square test matrix.
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
      if (lclRow < rowMap->getMaxLocalIndex()) {
        vals[0] = 1.;
      }
      else {
        vals[0] = 1.e-8;
      }
      A->insertLocalValues (lclRow, inds (), vals ());
    }
  }

  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;
  A->fillComplete (domMap, ranMap);
  return A;
}

// Create a very simple square test linear system (matrix, right-hand
// side(s), and exact solution(s).
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

  B->putScalar (STS::zero());
  auto map = B->getMap();
  // Set Vector 0 and Vector 1 at Global Row 0
  if (map->isNodeGlobalElement(0)) {
    // Get the local index for global row 0 on this processor
    size_t localRow = map->getLocalElement(0);
        
    // replaceLocalValue(localRowIndex, vectorIndex, value)
    B->replaceLocalValue(localRow, 0, 1.0); // First vector
    B->replaceLocalValue(localRow, 1, 2.0); // Second vector
  }

  // Set Vector 2 at Global Row 1
  if (map->isNodeGlobalElement(1)) {
    size_t localRow = map->getLocalElement(1);
    B->replaceLocalValue(localRow, 2, 1.0); // Third vector
  }

  X->putScalar (STS::zero());
}

template<class SC, class LO, class GO, class NT>
void
testSolver (Teuchos::FancyOStream& out,
            bool& success,
            const Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NT> >& X,
            const Teuchos::RCP<const Tpetra::CrsMatrix<SC, LO, GO, NT> >& A,
            const Teuchos::RCP<const Tpetra::MultiVector<SC, LO, GO, NT> >& B,
            Tpetra::MultiVector<SC, LO, GO, NT>& /* X_exact */,
            const std::string& solverName)
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::Operator<SC,LO,GO,NT> OP;
  typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
  typedef Teuchos::ScalarTraits<SC> STS;

  Teuchos::OSTab tab0 (out);
  out << "Test solver \"" << solverName << "\" from Belos package" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<Belos::SolverManager<SC, MV, OP> > solver;
  Belos::SolverFactory<SC, MV, OP> factory;

  // Set up Belos solver parameters.
  Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList (solverName);
  belosList->set ("Verbosity", Belos::Errors + Belos::Warnings);
  belosList->set("Maximum Iterations", 10);
  if (solverName == "BLOCK GMRES") {
    belosList->set("Flexible Gmres", false);
    belosList->set("Num Blocks", 1);
    belosList->set("Block Size", 3);
    belosList->set("Adaptive Block Size", false);
  }
  if (solverName == "GCRODR") {
    //belosList->set("Num Blocks", 2);
    //belosList->set("Num Recycled Blocks", 1);
  }
  belosList->set("Convergence Tolerance", 1.e-8);
  belosList->set("Orthogonalization", "ICGS");
  belosList->set("Orthogonalization Constant", 1.e-16);

  try {
    solver = factory.create (solverName, belosList);
  } catch (std::exception& e) {
    out << "*** FAILED: Belos::SolverFactory::create threw an exception: "
        << e.what () << endl;
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

  out << "Create the Belos::LinearProblem to solve" << endl;
  typedef Belos::LinearProblem<SC, MV, OP> linear_problem_type;
  X->putScalar (STS::zero ());

  RCP<linear_problem_type> problem (new linear_problem_type (A, X, B));
  problem->setProblem ();

  out << "Set up the solver" << endl;
  solver->setProblem (problem);

  out << "Apply solver to \"solve\" AX=B for X, and check if it fails to converge with 'OrthonormFailure'." << endl;
  Belos::ReturnType ret;

  try {
    ret = solver->solve ();
  } catch (std::exception& e) {
    out << "*** FAILED: Belos::SolverFactory::solve threw an unexpected exception: "
        << e.what () << endl;
    success = false;
    return;
  }
  out << "ret = " << convertReturnTypeToString(ret)
      << endl;

  if ( (ret != Belos::OrthonormFailure) ) {
    success = false;
    return;
  }

  success = true;
}

//
// The actual unit tests start here.
//

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( SolverFactory, CreateAndSolve, SC, LO, GO, NT )
{
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TypeNameTraits;
  using std::endl;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NT> MAT;
  typedef Tpetra::MultiVector<SC,LO,GO,NT> MV;
  typedef Tpetra::Operator<SC,LO,GO,NT> OP;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const Tpetra::global_size_t gblNumRows = comm->getSize () * 10;
  const size_t numVecs = 3;

  out << "Create test problem" << endl;
  RCP<MV> X_exact, B;
  RCP<MAT> A;
  createTestProblem (out, X_exact, A, B, comm, gblNumRows, numVecs);

  RCP<MV> X = rcp (new MV (X_exact->getMap (), numVecs));

  Belos::SolverFactory<SC, MV, OP> factory;
  const char* solverNames[2] = {
    "BLOCK GMRES",
    "GCRODR"};
  const int numSolvers = 2;

  int numSolversTested = 0;
  for (int k = 0; k < numSolvers; ++k) {
    const std::string solverName (solverNames[k]);
    out << "Testing k = " << k << ", solverName = " << solverName << std::endl;

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
      testSolver<SC, LO, GO, NT> (out, success, X, A, B, *X_exact, solverName);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( SolverFactory, CreateAndSolve, SC, LO, GO, NT )

// Tpetra's ETI will instantiate the unit test for all enabled type
// combinations.
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCLINST )

} // namespace (anonymous)
