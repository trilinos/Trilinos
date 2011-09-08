// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "ExampleTridiagSerialLinearOp.hpp"
#include "sillyCgSolve.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


//
// This example program is meant to show how easy it is to create
// serial Thyra objects and use them with an ANA (CG in this case).
//
// This example uses a silly concrete tridiagonal matrix class called
// ExampleSpmdTridiagLinearOp that demonstrates how to write and use such
// subclasses.
//
template<class Scalar>
bool runCgSolveExample(
  const int dim,
  const Scalar diagScale,
  const bool symOp,
  const bool showAllTests,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance,
  const int maxNumIters
  )
{

  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::multiply;
  using Thyra::scale;
  typedef typename ST::magnitudeType  ScalarMag;
  bool success = true;
  bool result;
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "\n***\n*** Running silly CG solver using scalar type = \'"
       << ST::name() << "\' ...\n***\n";
  Teuchos::Time timer("");
  timer.start(true);

  //
  // (A) Setup a simple linear system with tridiagonal operator:
  //
  //       [   a*2   -1                         ]
  //       [ -r(1)  a*2       -1                ]
  //  A =  [          .        .        .       ]
  //       [             -r(n-2)      a*2    -1 ]
  //       [                      -r(n-1)   a*2 ]
  //

  // (A.1) Create the tridiagonal matrix operator
  *out << "\nConstructing tridiagonal matrix A of dimension = " << dim
       << " and diagonal multiplier = " << diagScale << " ...\n";
  Teuchos::Array<Scalar> lower(dim-1), diag(dim), upper(dim-1);
  const Scalar
    up = -ST::one(),
    diagTerm = as<Scalar>(2.0) * diagScale * ST::one(),
    low = -(symOp ? ST::one() : ST::random());
  int k = 0;
  // First row
  diag[k] = diagTerm; upper[k] = up;
  // Middle rows
  for( k = 1; k < dim - 1; ++k ) {
    lower[k-1] = low; diag[k] = diagTerm; upper[k] = up;
  }
  // Last row
  lower[k-1] = low; diag[k] = diagTerm;
  RCP<const Thyra::LinearOpBase<Scalar> > A =
    rcp(new ExampleTridiagSerialLinearOp<Scalar>(dim, lower, diag, upper));

  // (A.2) Testing the linear operator constructed linear operator
  *out << "\nTesting the constructed linear operator A ...\n";
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.enable_all_tests(false);
  linearOpTester.check_linear_properties(true);
  linearOpTester.set_all_error_tol(tolerance);
  linearOpTester.set_all_warning_tol(1e-2*tolerance);
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

  // (A.5) Create the final linear system
  if(!symOp) {
    *out << "\nSetting up normal equations for unsymmetric system A^H*(A*x-b) => new A*x = b ...\n";
    // A^H*A
    RCP<const Thyra::LinearOpBase<Scalar> > AtA = multiply(adjoint(A), A);
    // A^H*b
    RCP<Thyra::VectorBase<Scalar> > nb = createMember(AtA->range());
    Thyra::apply<Scalar>(*A, Thyra::CONJTRANS, *b, nb.ptr());
    A = AtA;
    b = nb;
  }

  // (A.6) Testing the linear operator used with the solve
  *out << "\nTesting the linear operator used with the solve ...\n";
  linearOpTester.check_for_symmetry(true);
  result = linearOpTester.check(*A, out.ptr());
  if(!result) success = false;

  //
  // (B) Solve the linear system with the silly CG solver
  //
  *out << "\nSolving the linear system with sillyCgSolve(...) ...\n";
  {
    OSTab tab2(out);
    result = sillyCgSolve(*A, *b, maxNumIters, tolerance, x.ptr(), *out);
  }
  if(!result) success = false;

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

  bool success = true;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  // Above is needed to run in an MPI build with some MPI implementations

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read in command-line options
    //


    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    int dim = 500;
    clp.setOption( "dim", &dim,
      "Dimension of the linear system." );

    double diagScale = 1.001;
    clp.setOption( "diag-scale", &diagScale,
      "Scaling of the diagonal to improve conditioning." );

    bool symOp = true;
    clp.setOption( "sym-op", "unsym-op", &symOp,
      "Determines if the operator is symmetric or not." );

    bool showAllTests = false;
    clp.setOption( "show-all-tests", "show-summary-only", &showAllTests,
      "Show all LinearOpTester tests or not" );

    double tolerance = 1e-4;
    clp.setOption( "tol", &tolerance,
      "Relative tolerance for linear system solve." );

    int maxNumIters = 300;
    clp.setOption( "max-num-iters", &maxNumIters,
      "Maximum of CG iterations." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPTION( dim < 2, std::logic_error, "Error, dim=" << dim << " < 2 is not allowed!" );

#if defined(HAVE_THYRA_FLOAT)
    result = runCgSolveExample<float>(dim, diagScale, symOp, showAllTests,
      tolerance, maxNumIters);
    if(!result) success = false;
#endif

    result = runCgSolveExample<double>(dim, diagScale, symOp, showAllTests,
      tolerance, maxNumIters);
    if(!result) success = false;

#ifdef HAVE_THYRA_COMPLEX

#if defined(HAVE_THYRA_FLOAT)
    result = runCgSolveExample<std::complex<float> >(dim, diagScale, symOp, showAllTests,
      tolerance, maxNumIters);
    if(!result) success = false;
#endif

    result = runCgSolveExample<std::complex<double> >(dim, diagScale, symOp, showAllTests,
      tolerance, maxNumIters);
    if(!result) success = false;

#endif // HAVE_THYRA_COMPLEX

#ifdef HAVE_TEUCHOS_GNU_MP

    result = runCgSolveExample<mpf_class>(dim, diagScale, symOp, showAllTests, 
      tolerance, maxNumIters);
    if(!result) success = false;

#ifdef HAVE_THYRA_COMPLEX

    //result = runCgSolveExample<std::complex<mpf_class> >(dim, mpf_class(diagScale), symOp,
    showAllTests, mpf_class(tolerance), maxNumIters);
    //if(!result) success = false;
    //The above commented-out code throws a floating-point exception?

#endif // HAVE_THYRA_COMPLEX

#endif // HAVE_TEUCHOS_GNU_MP

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

    if(success)
      *out << "\nCongratulations! All of the tests checked out!\n";
    else
      *out << "\nOh no! At least one of the tests failed!\n";
  
  return success ? 0 : 1;

} // end main()
