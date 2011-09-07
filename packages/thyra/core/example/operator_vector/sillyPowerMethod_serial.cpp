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
  typedef typename ST::magnitudeType    ScalarMag;
  
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

    TEST_FOR_EXCEPTION( dim < 2, std::logic_error, "Error, dim=" << dim << " < 2 is not allowed!" );

    int maxNumIters = static_cast<int>(maxItersDimFactor*dim);
    
#if defined(HAVE_THYRA_FLOAT)
    // Run using float
    result = runPowerMethodExample<float>(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;
#endif

    // Run using double
    result = runPowerMethodExample<double>(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;

#ifdef HAVE_THYRA_COMPLEX
    
#if defined(HAVE_THYRA_FLOAT)
    // Run using std::complex<float>
    result = runPowerMethodExample<std::complex<float> >(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;
#endif
    
    // Run using std::complex<double>
    result = runPowerMethodExample<std::complex<double> >(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;

#endif // HAVE_THYRA_COMPLEX

#ifdef HAVE_TEUCHOS_GNU_MP
    
    // Run using mpf_class
    result = runPowerMethodExample<mpf_class >(
      dim, maxNumIters, tolerance, dumpAll);
    if(!result) success = false;

#ifdef HAVE_THYRA_COMPLEX
    
    // Run using std::complex<mpf_class>
    //result = runPowerMethodExample<std::complex<mpf_class> >(
    //  dim, maxNumIters, tolerance, dumpAll);
    //if(!result) success = false;
    //The above commented-out code throws a floating-point exception?
 
#endif // HAVE_THYRA_COMPLEX


#endif // HAVE_TEUCHOS_GNU_MP
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, *out, success)
  
  if (success)
    *out << "\nCongratulations! All of the tests checked out!\n";
  else
    *out << "\nOh no! At least one of the tests failed!\n";
  
  return success ? 0 : 1;

} // end main()
