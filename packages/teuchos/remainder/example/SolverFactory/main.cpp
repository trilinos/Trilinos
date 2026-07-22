// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Common.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

template<class ScalarType>
void testSolvers () {
  typedef Common::MultiVector<ScalarType> MV;
  typedef Common::Operator<ScalarType> OP;
  // To make the example simpler, we assume that ScalarType =
  // NormType.  For this to be correct, this would imply that
  // ScalarType is real.
  typedef ScalarType NormType;

  std::vector<std::pair<std::string, std::string> > solvers;
  solvers.push_back (std::make_pair ("A", "1"));
  solvers.push_back (std::make_pair ("A", "2"));
  solvers.push_back (std::make_pair ("B", "3"));
  solvers.push_back (std::make_pair ("B", "4"));
  solvers.push_back (std::make_pair ("C", "5"));
  solvers.push_back (std::make_pair ("C", "6"));

  for (size_t k = 0; k < solvers.size (); ++k) {
    const std::string packageName = solvers[k].first;
    const std::string solverName = solvers[k].second;
    cout << "Package \"" << packageName << "\", solver \"" << solverName
         << "\":" << endl;
    Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> > solver =
      Trilinos::Details::getLinearSolver<MV, OP, NormType> (packageName, solverName);
    if (solver.get () == NULL) {
      std::ostringstream err;
      err << "Solver \"" << solvers[k].second << "\" from package \""
          << solvers[k].first << "\" does not exist!";
      throw std::logic_error (err.str ());
    }

    Teuchos::RCP<OP> A = Teuchos::rcp (new OP ());
    // your code for filling in the matrix A would go here

    solver->setMatrix (A);
    solver->symbolic ();
    solver->numeric ();

    MV X, B;
    // your code for filling in X and B would go here

    solver->solve (X, B);
    cout << "Finished solver->solve(X, B)" << endl << endl;

    // This is a proxy for a residual calculation.  Some solvers
    // compute the residual on their own, but you have to ask them.
    A->apply (X, B);
  }
}


int main () {
  int err = EXIT_SUCCESS;

  cout << "Test ScalarType=float" << endl;
  try {
    testSolvers<float> ();
  } catch (std::exception& e) {
    cout << "testSolvers<float>() threw an exception: " << e.what () << endl;
    return EXIT_FAILURE;
  }

  cout << endl << "Test ScalarType=double" << endl;
  try {
    testSolvers<double> ();
  } catch (std::exception& e) {
    cout << "testSolvers<double>() threw an exception: " << e.what () << endl;
    return EXIT_FAILURE;
  }

  cout << endl << "Test ScalarType=int (should not work)" << endl;
  try {
    testSolvers<int> ();
    cout << "testSolvers<int>() should not have worked!" << endl;
    err = EXIT_FAILURE;
  } catch (std::exception&) {
    cout << "Of course testSolvers<int>() threw an exception: "
      "no packages registered themselves for ScalarType=int.  "
      "This is correct behavior in that case." << endl;
  }

  return err;
}

