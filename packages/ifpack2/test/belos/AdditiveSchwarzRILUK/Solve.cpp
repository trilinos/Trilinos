// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Solve.hpp"

//#include <Ifpack2_Factory.hpp>
#include <Ifpack2_AdditiveSchwarz.hpp>
#include <Ifpack2_RILUK.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

namespace Ifpack2 {
namespace Test {

void
solve (const sparse_mat_type& A,
       const multivector_type& b,
       multivector_type& x,
       const int ilukFillLevel,
       const int overlapLevel,
       const int numBlocks,
       const int maxIters,
       const double tol,
       const bool reorder,
       const std::string& innerPrecondType)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::cout;
  using std::endl;

  typedef Tpetra::RowMatrix<scalar_type, LO, GO, node_type> row_matrix_type;
  typedef Ifpack2::AdditiveSchwarz<row_matrix_type> outer_prec_type;
  typedef Tpetra::Operator<scalar_type, LO, GO, node_type> op_type;

  RCP<const Teuchos::Comm<int> > comm = A.getRowMap ()->getComm ();
  const int myRank = comm->getRank ();

  RCP<const sparse_mat_type> A_ptr = rcpFromRef (A);
  RCP<const multivector_type> b_ptr = rcpFromRef (b);
  RCP<multivector_type> x_ptr = rcpFromRef (x);

  // Set RILUK parameters
  ParameterList innerPlist;
  innerPlist.set ("fact: drop tolerance", 0.0);
  innerPlist.set ("fact: iluk level-of-fill", ilukFillLevel);
  innerPlist.set ("fact: relax value", 0.0);

  //Teuchos::RCP<prec_type> innerPrec;
  //Ifpack2::Factory factory;
  //innerPrec = factory.create (innerPrecondType, A_ptr);
  //innerPrec->setParameters (innerPlist);

  // Create (outer) Additive Schwarz preconditioner
  RCP<outer_prec_type> additiveSchwarz (new outer_prec_type (A_ptr));

  // Set outer preconditioner parameters
  ParameterList ASlist;
  ASlist.set ("inner preconditioner name", innerPrecondType);
  ASlist.set ("inner preconditioner parameters", innerPlist);
  ASlist.set ("schwarz: combine mode", "ZERO");
  ASlist.set ("schwarz: overlap level", overlapLevel);
  ASlist.set ("schwarz: use reordering", reorder);

  additiveSchwarz->setParameters (ASlist);

  // Compute (set up) the (outer) preconditioner
  additiveSchwarz->initialize ();
  additiveSchwarz->compute ();

  // Set GMRES (iterative linear solver) parameters
  RCP<ParameterList> solverParams = parameterList ();
  solverParams->set ("Num Blocks", numBlocks);
  solverParams->set ("Block Size", 1);
  solverParams->set ("Maximum Iterations", maxIters);
  solverParams->set ("Maximum Restarts", 10);
  solverParams->set ("Convergence Tolerance", tol);
  solverParams->set ("Implicit Residual Scaling", "None");
  solverParams->set ("Explicit Residual Scaling", "None");

  // Create the GMRES solver using a "factory" and
  // the list of solver parameters created above.
  typedef Belos::SolverFactory<scalar_type, multivector_type, op_type> belos_factory_type;
  typedef Belos::SolverManager<scalar_type, multivector_type, op_type> solver_type;
  RCP<solver_type> solver = belos_factory_type ().create ("GMRES", solverParams);

  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  typedef Belos::LinearProblem<scalar_type, multivector_type, op_type> problem_type;
  RCP<problem_type> problem (new problem_type (A_ptr, x_ptr, b_ptr));
  problem->setRightPrec (additiveSchwarz);

  // Tell the solver what problem you want to solve.
  solver->setProblem (problem);
  solver->reset (Belos::Problem);
  Belos::ReturnType result = solver->solve ();

  // Ask the solver how many iterations the last solve() took.
  const int numIters = solver->getNumIters ();

  if (myRank == 0) {
    if (solver->isLOADetected ()) {
      cout << "Detected a loss of accuracy!" << endl;
    }
    cout << "The Belos solve took " << numIters << " iteration"
         << (numIters != 1 ? "s" : "");
    if (result == Belos::Converged) {
      cout << " to converge." << endl;
    } else {
      cout << ", but did not converge." << endl;
    }
    cout << "It achieved a tolerance of: " << solver->achievedTol () << endl;
  }
}

} // namespace Test
} // namespace Ifpack2
