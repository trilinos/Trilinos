// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Solve.hpp"

#include <Ifpack2_AdditiveSchwarz.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <TpetraExt_MatrixMatrix.hpp>

namespace Ifpack2 {
namespace Test {

bool solveTwice(const sparse_mat_type& A,
                const sparse_mat_type& A2,
                const multivector_type& b,
                const multivector_type& b2,
                const int overlapLevel,
                const int numBlocks,
                const int maxIters,
                const double tol,
                const bool reorder,
                const std::string& innerPrecondType,
                const bool expectOneIter) {
  using std::cout;
  using std::endl;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  typedef Tpetra::RowMatrix<scalar_type, LO, GO, node_type> row_matrix_type;
  typedef Ifpack2::AdditiveSchwarz<row_matrix_type> outer_prec_type;
  typedef Tpetra::Operator<scalar_type, LO, GO, node_type> op_type;

  bool success                        = true;
  RCP<const Teuchos::Comm<int> > comm = A.getRowMap()->getComm();
  const int myRank                    = comm->getRank();

  RCP<sparse_mat_type> A_ptr  = Teuchos::rcp(new sparse_mat_type(A, Teuchos::Copy));
  RCP<multivector_type> b_ptr = Teuchos::rcp(new multivector_type(b, Teuchos::Copy));
  RCP<multivector_type> x_ptr = Teuchos::rcp(new multivector_type(A.getRowMap(), b.getNumVectors()));

  // Set parameters
  ParameterList innerPlist;
  if (innerPrecondType == "AMESOS2") {
    // Create the nested parameter list specifically for Amesos2
    innerPlist.setName("Amesos2");
    Teuchos::ParameterList& amesos2List = innerPlist.sublist("KLU2");
  }

  // Create (outer) Additive Schwarz preconditioner
  RCP<outer_prec_type> additiveSchwarz(new outer_prec_type(A_ptr));

  // Set outer preconditioner parameters
  ParameterList ASlist;
  ASlist.set("inner preconditioner name", innerPrecondType);
  ASlist.set("inner preconditioner parameters", innerPlist);
  ASlist.set("schwarz: combine mode", "ADD");
  ASlist.set("schwarz: overlap level", overlapLevel);
  ASlist.set("schwarz: use reordering", reorder);

  additiveSchwarz->setParameters(ASlist);

  // Compute (set up) the (outer) preconditioner
  additiveSchwarz->initialize();
  additiveSchwarz->compute();

  int verbosity = Belos::FinalSummary;
  // Set GMRES (iterative linear solver) parameters
  RCP<ParameterList> solverParams = parameterList();
  solverParams->set("Num Blocks", numBlocks);
  solverParams->set("Block Size", 1);
  solverParams->set("Maximum Iterations", maxIters);
  solverParams->set("Maximum Restarts", 10);
  solverParams->set("Convergence Tolerance", tol);
  solverParams->set("Implicit Residual Scaling", "None");
  solverParams->set("Explicit Residual Scaling", "None");
  solverParams->set("Verbosity", verbosity);

  // Create the GMRES solver using a "factory" and
  // the list of solver parameters created above.
  typedef Belos::SolverFactory<scalar_type, multivector_type, op_type> belos_factory_type;
  typedef Belos::SolverManager<scalar_type, multivector_type, op_type> solver_type;
  RCP<solver_type> solver = belos_factory_type().create("GMRES", solverParams);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  typedef Belos::LinearProblem<scalar_type, multivector_type, op_type> problem_type;
  RCP<problem_type> problem(new problem_type(A_ptr, x_ptr, b_ptr));
  problem->setRightPrec(additiveSchwarz);

  // Tell the solver what problem you want to solve.
  if (myRank == 0)
    cout << "Perform first solve with new Ifpack2 AdditiveSchwarz Amesos2 preconditioner" << endl;
  solver->setProblem(problem);
  solver->reset(Belos::Problem);
  Belos::ReturnType result = solver->solve();

  // Ask the solver how many iterations the last solve() took.
  int numIters = solver->getNumIters();
  if ((numIters != 1) && expectOneIter)
    success = false;

  if (myRank == 0) {
    if (solver->isLOADetected()) {
      cout << "Detected a loss of accuracy!" << endl;
    }
    cout << "The Belos solve took " << numIters << " iteration"
         << (numIters != 1 ? "s" : "");
    if (result == Belos::Converged) {
      cout << " to converge." << endl;
    } else {
      cout << ", but did not converge." << endl;
      success = false;
    }
    cout << "It achieved a tolerance of: " << solver->achievedTol() << endl
         << endl;
  }

  // -----------------------------------------------------------------------
  // Copy the values from A2 to A and b2 to b and perform second solve reusing preconditioner object
  // -----------------------------------------------------------------------
  Tpetra::MatrixMatrix::Add(A2, false, STS::one(), A2, false, STS::zero(), A_ptr);
  b_ptr->update(STS::one(), b2, STS::zero());
  x_ptr->putScalar(0.0);

  additiveSchwarz->compute();

  // Reset the iterative solver.
  if (myRank == 0)
    cout << "Perform second solve reusing Ifpack2 AdditiveSchwarz Amesos2 preconditioner" << endl;
  solver->reset(Belos::Problem);
  result = solver->solve();

  // Ask the solver how many iterations the last solve() took.
  int numIters2 = solver->getNumIters();
  if ((numIters2 != 1) && expectOneIter)
    success = false;

  if (myRank == 0) {
    if (solver->isLOADetected()) {
      cout << "Detected a loss of accuracy!" << endl;
    }
    cout << "The Belos solve took " << numIters2 << " iteration"
         << (numIters2 != 1 ? "s" : "");
    if (result == Belos::Converged) {
      cout << " to converge." << endl;
    } else {
      cout << ", but did not converge." << endl;
      success = false;
    }
    cout << "It achieved a tolerance of: " << solver->achievedTol() << endl
         << endl;
  }

  // -----------------------------------------------------------------------
  // Now recreate the preconditioner to check if the number of iterations is the same
  // -----------------------------------------------------------------------
  additiveSchwarz = Teuchos::rcp(new outer_prec_type(A_ptr));
  additiveSchwarz->setParameters(ASlist);

  // Compute (set up) the (outer) preconditioner
  additiveSchwarz->initialize();
  additiveSchwarz->compute();
  problem->setRightPrec(additiveSchwarz);

  // Reset the iterative solver.
  if (myRank == 0)
    cout << "Perform second solve with new Ifpack2 AdditiveSchwarz Amesos2 preconditioner" << endl;
  x_ptr->putScalar(0.0);
  solver->reset(Belos::Problem);
  result = solver->solve();

  // Ask the solver how many iterations the last solve() took.
  int numIters3 = solver->getNumIters();
  if ((numIters3 != 1) && expectOneIter)
    success = false;

  if (myRank == 0) {
    if (solver->isLOADetected()) {
      cout << "Detected a loss of accuracy!" << endl;
    }
    cout << "The Belos solve took " << numIters3 << " iteration"
         << (numIters3 != 1 ? "s" : "");
    if (result == Belos::Converged) {
      cout << " to converge." << endl;
    } else {
      cout << ", but did not converge." << endl;
      success = false;
    }
    cout << "It achieved a tolerance of: " << solver->achievedTol() << endl
         << endl;
    if (numIters2 != numIters3) {
      cout << "Ifpack2::AdditiveSchwarz::Amesos2 preconditioner resulted in different iteration count: "
           << numIters2 << " vs. " << numIters3 << endl;
    }
  }

  return success;
}

}  // namespace Test
}  // namespace Ifpack2
