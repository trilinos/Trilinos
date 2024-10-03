// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "sillyPowerMethod.hpp"
#include "ExampleTridiagSerialLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"

//
// This example function is meant to show how easy it is to create
// serial Thyra objects and use them with an ANA.
//
// This example uses a silly concrete tridiagonal matrix class
// called SillyTridiagSerialLinearOp that demonstrates how
// to write such subclasses.
//
template<class Scalar>
bool runPowerMethodExample(
  const int dim,
  const int maxNumIters,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance,
  const bool  dumpAll
  )
{

  using Teuchos::OSTab;
  using Teuchos::outArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  //typedef typename ST::magnitudeType    ScalarMag; //unused

  bool success = true;
  bool result;

  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();

  *out << "\n***\n*** Running power method example using scalar type = \'"
    << ST::name() << "\' ...\n***\n" << std::scientific;

  //
  // (1) Setup the initial tridiagonal operator
  //
  //       [  2  -1             ]
  //       [ -1   2  -1         ]
  //  A =  [      .   .   .     ]
  //       [          -1  2  -1 ]
  //       [             -1   2 ]
  //
  *out << "\n(1) Constructing tridiagonal matrix A of dimension = " << dim << " ...\n";
  Teuchos::Array<Scalar> lower(dim-1), diag(dim), upper(dim-1);
  const Scalar one = ST::one(), two = Scalar(2)*one;
  int k = 0;
  diag[k] = two; upper[k] = -one;                        //  First row
  for( k = 1; k < dim - 1; ++k ) {
    lower[k-1] = -one; diag[k] = two; upper[k] = -one;   //  Middle rows
  }
  lower[k-1] = -one; diag[k] = two;                      //  Last row
  Teuchos::RCP<ExampleTridiagSerialLinearOp<Scalar> >
    A = Teuchos::rcp(new ExampleTridiagSerialLinearOp<Scalar>(dim, lower, diag, upper));
  if (dumpAll) *out << "\nA =\n" << *A;

  //
  // (2) Run the power method ANA
  //
  *out << "\n(2) Running the power method on matrix A ...\n";
  Scalar     lambda      = ST::zero();
  {
    OSTab tab(out);
    result = sillyPowerMethod(*A, maxNumIters, tolerance, outArg(lambda), *out);
    if(!result) success = false;
    *out << "\nEstimate of dominate eigenvalue lambda = " << lambda << std::endl;
  }

  //
  // (3) Increase dominance of first eigenvalue
  //
  *out << "\n(3) Increasing first diagonal entry by factor of 10 ...\n";
  diag[0] *= 10;
  A->initialize(dim, lower, diag, upper);
  if (dumpAll) *out << "A =\n" << *A;

  //
  // (4) Run the power method ANA
  //
  *out << "\n(4) Running the power method again on matrix A ...\n";
  {
    OSTab tab(out);
    result = sillyPowerMethod(*A, maxNumIters, tolerance, outArg(lambda), *out);
    if(!result) success = false;
    *out << "\nEstimate of dominate eigenvalue lambda = " << lambda << std::endl;
  }

  return success;

} // end runPowerMethodExample()


//
// Main driver program for serial implementation of the power method.
//
// Note that several different scalar types are used here (float,
// double, std::complex<float>, std::complex<double> etc.).
//
int main(int argc, char *argv[])
{

  using Teuchos::as;
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

    int dim = 4;
    clp.setOption( "dim", &dim, "Dimension of the linear system." );

    bool dumpAll = false;
    clp.setOption( "dump-all", "no-dump", &dumpAll,
      "Determines if quantities are dumped or not." );

    double tolerance = 1e-3;
    clp.setOption( "tol", &tolerance, "Final tolerance of eigen system." );

    double maxItersDimFactor = 10.0;
    clp.setOption( "max-iters-dim-factor", &maxItersDimFactor,
      "Factor to multiple dim by to get maxIters." );

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);
    if (parse_return != CommandLineProcessor::PARSE_SUCCESSFUL)
      return parse_return;

    TEUCHOS_TEST_FOR_EXCEPTION( dim < 2, std::logic_error, "Error, dim=" << dim << " < 2 is not allowed!" );

    int maxNumIters = static_cast<int>(maxItersDimFactor*dim);

#if defined(HAVE_TEUCHOS_INST_FLOAT)
    // Run using float
    result = runPowerMethodExample<float>(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;
#endif

    // Run using double
    result = runPowerMethodExample<double>(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;

#if defined(HAVE_TEUCHOS_INST_COMPLEX_FLOAT) && defined(HAVE_TEUCHOS_INST_FLOAT)
    // Run using std::complex<float>
    result = runPowerMethodExample<std::complex<float> >(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;
#endif

#if defined(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE)
    // Run using std::complex<double>
    result = runPowerMethodExample<std::complex<double> >(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;
#endif

#ifdef HAVE_TEUCHOS_GNU_MP

    // Run using mpf_class
    result = runPowerMethodExample<mpf_class >(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;

    // Run using std::complex<mpf_class>
    //result = runPowerMethodExample<std::complex<mpf_class> >(
    //  dim, maxNumIters, tolerance, dumpAll);
    //if(!result) success = false;
    //The above commented-out code throws a floating-point exception?

#endif // HAVE_TEUCHOS_GNU_MP

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success)

  if (success)
    *out << "\nCongratulations! All of the tests checked out!\n";
  else
    *out << "\nOh no! At least one of the tests failed!\n";

  return success ? 0 : 1;

} // end main()
