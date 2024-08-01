// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace { // (anonymous)

// 18x18 symmetric positive definite matrix from the University of
// Florida sparse matrix collection: Oberwolfach/LF10 (the main (K)
// matrix).
const char symPosDefMatrixString[] =
"%%MatrixMarket matrix coordinate real symmetric\n"
"%-------------------------------------------------------------------------------\n"
"% UF Sparse Matrix Collection, Tim Davis\n"
"% http://www.cise.ufl.edu/research/sparse/matrices/Oberwolfach/LF10\n"
"% name: Oberwolfach/LF10\n"
"% [Oberwolfach: linear 1D beam]\n"
"% id: 1438\n"
"% date: 2004\n"
"% author: J. Lienemann, A. Greiner, J. Korvink\n"
"% ed: E. Rudnyi\n"
"% fields: name title A id notes aux date author ed kind\n"
"% aux: M E B C\n"
"% kind: model reduction problem\n"
"%-------------------------------------------------------------------------------\n"
"% notes:\n"
"% Primary matrix in this model reduction problem is the Oberwolfach K matrix\n"
"%-------------------------------------------------------------------------------\n"
"18 18 50\n"
"1 1 3.5344800000000003\n"
"2 1 -477.1548\n"
"3 1 1.7672400000000001\n"
"2 2 171775.728\n"
"4 2 -85887.864\n"
"5 2 477.1548\n"
"3 3 7.068960000000001\n"
"4 3 -477.1548\n"
"5 3 1.7672400000000001\n"
"4 4 171775.728\n"
"6 4 -85887.864\n"
"7 4 477.1548\n"
"5 5 7.068960000000001\n"
"6 5 -477.1548\n"
"7 5 1.7672400000000001\n"
"6 6 171775.728\n"
"8 6 -85887.864\n"
"9 6 477.1548\n"
"7 7 7.068960000000001\n"
"8 7 -477.1548\n"
"9 7 1.7672400000000001\n"
"8 8 171775.728\n"
"10 8 -85887.864\n"
"11 8 477.1548\n"
"9 9 7.068960000000001\n"
"10 9 -477.1548\n"
"11 9 1.7672400000000001\n"
"10 10 171775.728\n"
"12 10 -85887.864\n"
"13 10 477.1548\n"
"11 11 7.068960000000001\n"
"12 11 -477.1548\n"
"13 11 1.7672400000000001\n"
"12 12 171775.728\n"
"14 12 -85887.864\n"
"15 12 477.1548\n"
"13 13 7.068960000000001\n"
"14 13 -477.1548\n"
"15 13 1.7672400000000001\n"
"14 14 171775.728\n"
"16 14 -85887.864\n"
"17 14 477.1548\n"
"15 15 7.068960000000001\n"
"16 15 -477.1548\n"
"17 15 1.7672400000000001\n"
"16 16 171775.728\n"
"18 16 477.1548\n"
"17 17 7.068960000000001\n"
"18 17 1.7672400000000001\n"
"18 18 3.5344800000000003\n";

#if 0
// Test matrix: nonsymmetric real.  It's small enough to store as a
// string, so we don't have to worry about adding a file to the
// repository.
const char matrixString[] =
"%%MatrixMarket matrix coordinate real general\n"
"%-------------------------------------------------------------------------------\n"
"% UF Sparse Matrix Collection, Tim Davis\n"
"    % http://www.cise.ufl.edu/research/sparse/matrices/vanHeukelum/cage4\n"
"    % name: vanHeukelum/cage4\n"
"    % [DNA electrophoresis, 4 monomers in polymer. A. van Heukelum, Utrecht U.]\n"
"    % id: 905\n"
"    % date: 2003\n"
"    % author: A. van Heukelum\n"
"    % ed: T. Davis\n"
"    % fields: title A name id date author ed kind\n"
"    % kind: directed weighted graph\n"
"%-------------------------------------------------------------------------------\n"
"9 9 49\n"
"1 1 .75\n"
"2 1 .075027667114587\n"
"4 1 .0916389995520797\n"
"5 1 .0375138335572935\n"
"8 1 .0458194997760398\n"
"1 2 .137458499328119\n"
"2 2 .687569167786467\n"
"3 2 .0916389995520797\n"
"5 2 .0375138335572935\n"
"6 2 .0458194997760398\n"
"2 3 .112541500671881\n"
"3 3 .666666666666667\n"
"4 3 .13745849932812\n"
"6 3 .0458194997760398\n"
"7 3 .0375138335572935\n"
"1 4 .112541500671881\n"
"3 4 .075027667114587\n"
"4 4 .729097498880199\n"
"7 4 .0375138335572935\n"
"8 4 .0458194997760398\n"
"1 5 .137458499328119\n"
"2 5 .075027667114587\n"
"5 5 .537513833557293\n"
"6 5 .075027667114587\n"
"7 5 .0916389995520797\n"
"9 5 .0833333333333333\n"
"2 6 .112541500671881\n"
"3 6 .0916389995520797\n"
"5 6 .13745849932812\n"
"6 6 .445874834005214\n"
"8 6 .13745849932812\n"
"9 6 .075027667114587\n"
"3 7 .075027667114587\n"
"4 7 .13745849932812\n"
"5 7 .112541500671881\n"
"7 7 .470791832661453\n"
"8 7 .112541500671881\n"
"9 7 .0916389995520797\n"
"1 8 .112541500671881\n"
"4 8 .0916389995520797\n"
"6 8 .075027667114587\n"
"7 8 .0916389995520797\n"
"8 8 .54581949977604\n"
"9 8 .0833333333333333\n"
"5 9 .25\n"
"6 9 .150055334229174\n"
"7 9 .183277999104159\n"
"8 9 .25\n"
"9 9 .166666666666667\n";
#endif // 0

using Teuchos::Comm;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using std::endl;
typedef Tpetra::MultiVector<>     MV;
typedef MV::scalar_type           ST;
typedef MV::mag_type              MT;
typedef Teuchos::ScalarTraits<ST> STS;
typedef Teuchos::ScalarTraits<MT> STM;
typedef Tpetra::Operator<>        OP;
typedef Tpetra::CrsMatrix<> crs_matrix_type;

// Report both computed-solution and exact-solution residuals.
void
testResiduals (bool& success,
               Teuchos::FancyOStream& out,
               const Tpetra::Operator<>& A,
               const Tpetra::MultiVector<>& X,
               const Tpetra::MultiVector<>& B,
               const Tpetra::MultiVector<>& X_exact,
               const MT& tol,
               const Belos::ReturnType& ret)
{
  const size_t numRhs = X.getNumVectors ();

  // R := B - A*X
  auto R = rcp (new MV (B, Teuchos::Copy));
  A.apply (X, *R, Teuchos::NO_TRANS, -STS::one (), STS::one ());
  Teuchos::Array<MT> norms (numRhs);
  R->norm2 (norms ());

  // R_exact := B - A*X_exact
  auto R_exact = rcp (new MV (B, Teuchos::Copy));
  A.apply (X_exact, *R_exact, Teuchos::NO_TRANS, -STS::one (), STS::one ());
  Teuchos::Array<MT> normsExact (numRhs);
  R_exact->norm2 (normsExact ());

  Teuchos::Array<MT> B_norms (numRhs);
  B.norm2 (B_norms ());

  bool badRes = false;
  out << "Results:" << endl;
  {
    Teuchos::OSTab tab1 (out);

    out << "Exact solution relative residuals: [";
    for (size_t j = 0; j < numRhs; ++j) {
      const MT exactRelRes = normsExact[j] / B_norms[j];
      out << exactRelRes;
      if (j + 1 < numRhs) {
        out << ", ";
      }
    }
    out << "]" << endl;

    out << "Computed solution relative residuals: [";
    for (size_t j = 0; j < numRhs; ++j) {
      const MT computedRelRes = norms[j] / B_norms[j];
      out << computedRelRes;
      if (j + 1 < numRhs) {
        out << ", ";
      }
      if (computedRelRes > tol) {
        badRes = true;
      }
    }
    out << "]" << endl;
  }

  if (ret != Belos::Converged || badRes) {
    success = false;
    out << "Belos did NOT converge!" << endl;
  } else {
    out << "Belos converged!" << endl;
  }
}


TEUCHOS_UNIT_TEST( MultipleSolves, GMRES )
{
  Teuchos::OSTab tab0 (out);
  out << "Test whether Belos correctly overwrites the initial guess with the "
      << "computed solution, or incorrectly does a +=, as some users have "
      << "reported" << endl;
  Teuchos::OSTab tab1 (out);

  // Get the default communicator for tests.
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  std::istringstream matrixFile (symPosDefMatrixString);
  RCP<crs_matrix_type> A = reader_type::readSparse (matrixFile, comm);

  size_t numRhs = 1; // number of right-hand sides in the linear system
  int freq = 1; // frequency of status test output

  // Maximum number of iterations per solve.  We use twice the problem
  // size, since the problem to solve is pretty hard for being 18x18.
  const int maxIters =
    2 * static_cast<int> (A->getRangeMap ()->getGlobalNumElements ());
  // Relative residual tolerance
  const MT tol = static_cast<MT> (10) * STM::squareroot (STM::eps ());

  const std::string solverName ("GMRES");
  RCP<ParameterList> belosList (new ParameterList ("Belos"));
  belosList->set ("Maximum Iterations", maxIters);
  belosList->set ("Convergence Tolerance", tol);
  belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                  /* Belos::StatusTestDetails + */ Belos::FinalSummary);
  belosList->set ("Output Frequency", freq);

  auto X = rcp (new MV (A->getDomainMap (), numRhs));
  auto X_exact = rcp (new MV (A->getDomainMap (), numRhs));
  auto B = rcp (new MV (A->getRangeMap (), numRhs));

  // Fill initial guess with zeros.
  X->putScalar (STS::one ());

  //
  // First solve: Exact solution is all ones.
  //
  X_exact->putScalar (STS::one ());
  A->apply (*X_exact, *B);

  typedef Belos::LinearProblem<ST, MV, OP> problem_type;
  auto problem = rcp (new problem_type (A, X, B));
  bool set = problem->setProblem ();
  TEUCHOS_TEST_FOR_EXCEPTION( ! set, std::logic_error, "Failed to set Belos::LinearProblem!");

  // Create a Belos solver.
  RCP<Belos::SolverManager<ST, MV, OP> > solver;
  {
    Belos::SolverFactory<ST, MV, OP> factory;
    solver = factory.create (solverName, belosList);
    solver->setProblem (problem);
  }

  // Ask Belos to solve the linear system.
  Belos::ReturnType ret = solver->solve ();

  // Evaluate the result.
  testResiduals (success, out, *A, *X, *B, *X_exact, tol, ret);

  //
  // Second solve: Exact solution is all 1.5.  We pick this because
  // it's not too far off from the previous exact solution -- e.g.,
  // not the opposite sign -- so it shouldn't take too long to
  // converge.
  //
  const ST one = STS::one ();
  const ST two = one + one;
  X_exact->putScalar (one + one / two);
  A->apply (*X_exact, *B);

  // Don't change the initial guess this time.  Keep the solution from
  // last time.  The point is to make sure that the LinearProblem
  // doesn't do something weird to it.

  set = problem->setProblem ();
  TEUCHOS_TEST_FOR_EXCEPTION( ! set, std::logic_error, "Failed to set Belos::LinearProblem!");

  // Ask Belos to solve the linear system.
  ret = solver->solve ();

  // Evaluate the result.
  testResiduals (success, out, *A, *X, *B, *X_exact, tol, ret);
}


TEUCHOS_UNIT_TEST( MultipleSolves, CG )
{
  Teuchos::OSTab tab0 (out);
  out << "Test whether Belos correctly overwrites the initial guess with the "
      << "computed solution, or incorrectly does a +=, as some users have "
      << "reported" << endl;
  Teuchos::OSTab tab1 (out);

  // Get the default communicator for tests.
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  std::istringstream matrixFile (symPosDefMatrixString);
  RCP<crs_matrix_type> A = reader_type::readSparse (matrixFile, comm);

  size_t numRhs = 1; // number of right-hand sides in the linear system
  int freq = 1; // frequency of status test output

  // Maximum number of iterations per solve.  We use twice the problem
  // size, since the problem to solve is pretty hard for being 18x18.
  const int maxIters =
    2 * static_cast<int> (A->getRangeMap ()->getGlobalNumElements ());
  // Relative residual tolerance
  const MT tol = static_cast<MT> (10) * STM::squareroot (STM::eps ());

  const std::string solverName ("CG");
  RCP<ParameterList> belosList (new ParameterList ("Belos"));
  belosList->set ("Maximum Iterations", maxIters);
  belosList->set ("Convergence Tolerance", tol);
  belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                  /* Belos::StatusTestDetails + */ Belos::FinalSummary);
  belosList->set ("Output Frequency", freq);

  auto X = rcp (new MV (A->getDomainMap (), numRhs));
  auto X_exact = rcp (new MV (A->getDomainMap (), numRhs));
  auto B = rcp (new MV (A->getRangeMap (), numRhs));

  // Fill initial guess with zeros.
  X->putScalar (STS::one ());

  //
  // First solve: Exact solution is all ones.
  //
  X_exact->putScalar (STS::one ());
  A->apply (*X_exact, *B);

  typedef Belos::LinearProblem<ST, MV, OP> problem_type;
  auto problem = rcp (new problem_type (A, X, B));
  bool set = problem->setProblem ();
  TEUCHOS_TEST_FOR_EXCEPTION( ! set, std::logic_error, "Failed to set Belos::LinearProblem!");

  // Create a Belos solver.
  RCP<Belos::SolverManager<ST, MV, OP> > solver;
  {
    Belos::SolverFactory<ST, MV, OP> factory;
    solver = factory.create (solverName, belosList);
    solver->setProblem (problem);
  }

  // Ask Belos to solve the linear system.
  Belos::ReturnType ret = solver->solve ();

  // Evaluate the result.
  testResiduals (success, out, *A, *X, *B, *X_exact, tol, ret);

  //
  // Second solve: Exact solution is all 1.5.  We pick this because
  // it's not too far off from the previous exact solution -- e.g.,
  // not the opposite sign -- so it shouldn't take too long to
  // converge.
  //
  const ST one = STS::one ();
  const ST two = one + one;
  X_exact->putScalar (one + one / two);
  A->apply (*X_exact, *B);

  // Don't change the initial guess this time.  Keep the solution from
  // last time.  The point is to make sure that the LinearProblem
  // doesn't do something weird to it.

  set = problem->setProblem ();
  TEUCHOS_TEST_FOR_EXCEPTION( ! set, std::logic_error, "Failed to set Belos::LinearProblem!");

  // Ask Belos to solve the linear system.
  ret = solver->solve ();

  // Evaluate the result.
  testResiduals (success, out, *A, *X, *B, *X_exact, tol, ret);
}

} // namespace (anonymous)
