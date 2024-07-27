// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ExampleTridiagSpmdLinearOp.hpp"
#include "sillyCgSolve.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


//
// This example program is meant to show how easy it is to create MPI Thyra
// objects and use them with an ANA (CG in this case).
//
// This example uses a silly concrete tridiagonal matrix class called
// ExampleTridiagSpmdLinearOp that demonstrates how to write such subclasses.
//
template<class Scalar>
bool runCgSolveExample(
  const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm,
  const int procRank,
  const int numProc,
  const int localDim,
  const Scalar diagScale,
  const bool showAllTests,
  const bool dumpAll,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance,
  const int maxNumIters
  )
{

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  bool success = true;
  bool result;

  // ToDo: Get VerboseObjectBase to automatically setup for parallel
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "\n***\n*** Running silly CG solver using scalar type = \'" << ST::name() << "\' ...\n***\n";
  Teuchos::Time timer("");
  timer.start(true);

  //
  // (A) Setup a simple linear system with tridiagonal operator:
  //
  //       [  a*2   -1                ]
  //       [ -1    a*2  -1            ]
  //  A =  [         .   .    .       ]
  //       [            -1  a*2    -1 ]
  //       [                 -1   a*2 ]
  //

  // (A.1) Create the tridiagonal matrix operator
  *out << "\nConstructing tridiagonal matrix A of local dimension = " << localDim
       << " and diagonal multiplier = " << diagScale << " ...\n";
  const Thyra::Ordinal
    lowerDim = ( procRank == 0         ? localDim - 1 : localDim ),
    upperDim = ( procRank == numProc-1 ? localDim - 1 : localDim );
  Teuchos::Array<Scalar> lower(lowerDim), diag(localDim), upper(upperDim);
  const Scalar one = ST::one(), diagTerm = as<Scalar>(2.0)* diagScale * ST::one();
  int k = 0, kl = 0;
  // First local row
  if(procRank > 0) { lower[kl] = -one; ++kl; };  diag[k] = diagTerm; upper[k] = -one;
  // Middle local rows
  for( k = 1; k < localDim - 1; ++k, ++kl ) {
    lower[kl] = -one; diag[k] = diagTerm; upper[k] = -one;
  }
  // Last local row
  lower[kl] = -one; diag[k] = diagTerm; if(procRank < numProc-1) upper[k] = -one;
  RCP<const Thyra::LinearOpBase<Scalar> > A =
    rcp(new ExampleTridiagSpmdLinearOp<Scalar>(comm, localDim, lower, diag, upper));
  *out << "\nGlobal dimension of A = " << A->domain()->dim() << std::endl;

  // (A.2) Testing the linear operator constructed linear operator
  *out << "\nTesting the constructed linear operator A ...\n";
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.dump_all(dumpAll);
  linearOpTester.set_all_error_tol(tolerance);
  linearOpTester.set_all_warning_tol(1e-2*tolerance);
  linearOpTester.show_all_tests(true);
  linearOpTester.check_adjoint(false);
  linearOpTester.check_for_symmetry(true);
  linearOpTester.show_all_tests(showAllTests);
  result = linearOpTester.check(*A, out.ptr());
  if(!result) success = false;

  // (A.3) Create RHS vector b and set to a random value
  RCP<Thyra::VectorBase<Scalar> > b = createMember(A->range());
  Thyra::seed_randomize<Scalar>(0);
  Thyra::randomize( -ST::one(), +ST::one(), b.ptr() );

  // (A.4) Create LHS vector x and set to zero
  RCP<Thyra::VectorBase<Scalar> > x = createMember(A->domain());
  Thyra::V_S( x.ptr(), ST::zero() );

  //
  // (B) Solve the linear system with the silly CG solver
  //
  *out << "\nSolving the linear system with sillyCgSolve(...) ...\n";
  {
    OSTab tab(out);
    result = sillyCgSolve(*A, *b, maxNumIters, tolerance, x.ptr(), *out);
    if(!result) success = false;
  }

  //
  // (C) Check that the linear system was solved to the specified tolerance
  //
  RCP<Thyra::VectorBase<Scalar> > r = createMember(A->range());
  // r = b
  Thyra::V_V(r.ptr(), *b);
  // r = -A*x + r
  Thyra::apply<Scalar>(*A, Thyra::NOTRANS, *x, r.ptr(), -ST::one(), ST::one());
  const ScalarMag r_nrm = Thyra::norm(*r), b_nrm = Thyra::norm(*b);
  const ScalarMag rel_err = r_nrm/b_nrm, relaxTol = 10.0*tolerance;
  result = rel_err <= relaxTol;
  if(!result) success = false;
  *out << "\nChecking the residual ourselves ...\n";
  {
    OSTab tab(out);
    *out
      << "\n||b-A*x||/||b|| = "<<r_nrm<<"/"<<b_nrm<<" = "<<rel_err<<(result?" <= ":" > ")
      <<"10.0*tolerance = "<<relaxTol<<": "<<(result?"passed":"failed")<<std::endl;
  }
  timer.stop();
  *out << "\nTotal time = " << timer.totalElapsedTime() << " sec\n";

  return success;

} // end runCgSolveExample()


//
// Actual main driver program
//
int main(int argc, char *argv[])
{

  using Teuchos::CommandLineProcessor;

  bool success = false;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {
    bool result;
    const int procRank = Teuchos::GlobalMPISession::getRank();
    const int numProc = Teuchos::GlobalMPISession::getNProc();

    const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >
      comm = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();

    //
    // Read in command-line options
    //

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    int localDim = 500;
    clp.setOption( "local-dim", &localDim,
      "Local dimension of the linear system." );

    double diagScale = 1.001;
    clp.setOption( "diag-scale", &diagScale,
      "Scaling of the diagonal to improve conditioning." );

    bool showAllTests = false;
    clp.setOption( "show-all-tests", "show-summary-only", &showAllTests,
      "Show all LinearOpTester tests or not" );

    double tolerance = 1e-4;
    clp.setOption( "tol", &tolerance,
      "Relative tolerance for linear system solve." );

    int maxNumIters = 300;
    clp.setOption( "max-num-iters", &maxNumIters,
      "Maximum of CG iterations." );

    bool dumpAll = false;
    clp.setOption( "dump-all", "no-dump-all", &dumpAll,
      "Determines if vectors are printed or not." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc, argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEUCHOS_TEST_FOR_EXCEPTION( localDim < 2, std::logic_error,
      "Error, localDim=" << localDim << " < 2 is not allowed!" );

#ifdef HAVE_TEUCHOS_INST_FLOAT
    result = runCgSolveExample<float>(comm, procRank, numProc, localDim, diagScale,
      showAllTests, dumpAll, tolerance, maxNumIters);
    success = result;
#endif

    result = runCgSolveExample<double>(comm, procRank, numProc, localDim, diagScale,
      showAllTests, dumpAll, tolerance, maxNumIters);
    success = result;

#if defined(HAVE_TEUCHOS_INST_COMPLEX_FLOAT) && defined(HAVE_TEUCHOS_INST_FLOAT)
    result = runCgSolveExample<std::complex<float> >(comm, procRank, numProc, localDim,
      diagScale, showAllTests, dumpAll, tolerance, maxNumIters);
    success = result;
#endif

#if defined(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE)
    result = runCgSolveExample<std::complex<double> >(comm, procRank, numProc, localDim,
      diagScale, showAllTests, dumpAll, tolerance, maxNumIters);
    success = result;
#endif

    if (procRank==0) {
      if (success)
        *out << "\nAll of the tests seem to have run successfully!\n";
      else
        *out << "\nOh no! at least one of the tests failed!\n";
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success)

  return success ? EXIT_SUCCESS : EXIT_FAILURE;

} // end main()
