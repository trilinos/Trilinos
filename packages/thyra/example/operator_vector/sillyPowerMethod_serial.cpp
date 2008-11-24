// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "sillyPowerMethod.hpp"
#include "ExampleTridiagSerialLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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
  const int    dim
  ,const int   maxNumIters
  ,const bool  verbose
  ,const bool  dumpAll
  )
{

  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;
  
  bool success = true;
  bool result;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = (verbose ? Teuchos::VerboseObjectBase::getDefaultOStream() : Teuchos::null);

  if(verbose)
    *out << "\n***\n*** Running power method example using scalar type = \'" << ST::name() << "\' ...\n***\n" << std::scientific;

  //
  // (1) Setup the initial tridiagonal operator
  //
  //       [  2  -1             ]
  //       [ -1   2  -1         ]
  //  A =  [      .   .   .     ]
  //       [          -1  2  -1 ]
  //       [             -1   2 ]
  //
  if(verbose) *out << "\n(1) Constructing tridiagonal matrix A of dimension = " << dim << " ...\n";
  std::vector<Scalar> lower(dim-1), diag(dim), upper(dim-1);
  const Scalar one = ST::one(), two = Scalar(2)*one;
  int k = 0;
  diag[k] = two; upper[k] = -one;                        //  First row
  for( k = 1; k < dim - 1; ++k ) {
    lower[k-1] = -one; diag[k] = two; upper[k] = -one;   //  Middle rows
  }
  lower[k-1] = -one; diag[k] = two;                      //  Last row
  Teuchos::RCP<ExampleTridiagSerialLinearOp<Scalar> >
    A = Teuchos::rcp( new ExampleTridiagSerialLinearOp<Scalar>(dim,&lower[0],&diag[0],&upper[0]) );
  if( verbose && dumpAll ) *out << "\nA =\n" << *A;

  //
  // (2) Run the power method ANA
  //
  if(verbose) *out << "\n(2) Running the power method on matrix A ...\n";
  Scalar     lambda      = ST::zero();
  ScalarMag  tolerance   = ScalarMag(1e-3)*Teuchos::ScalarTraits<ScalarMag>::one();
  {
    OSTab tab(out);
    result = sillyPowerMethod(*A,maxNumIters,tolerance,&lambda,out.get());
    if(!result) success = false;
    if(verbose) *out << "\nEstimate of dominate eigenvalue lambda = " << lambda << std::endl;
  }

  //
  // (3) Increase dominance of first eigenvalue
  //
  if(verbose) *out << "\n(3) Increasing first diagonal entry by factor of 10 ...\n";
  diag[0] *= 10;
  A->initialize(dim,&lower[0],&diag[0],&upper[0]);
  if( verbose && dumpAll ) *out << "A =\n" << *A;

  //
  // (4) Run the power method ANA
  //
  if(verbose) *out << "\n(4) Running the power method again on matrix A ...\n";
  {
    OSTab tab(out);
    result = sillyPowerMethod(*A,maxNumIters,tolerance,&lambda,out.get());
    if(!result) success = false;
    if(verbose) *out << "\nEstimate of dominate eigenvalue lambda = " << lambda << std::endl;
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

  using Teuchos::CommandLineProcessor;
 
  bool success = true;
  bool verbose = true;
  bool result;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read in command-line options
    //
    
    int    dim          = 4;
    bool   dumpAll      = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "dim", &dim, "Dimension of the linear system." );
    clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPTION( dim < 2, std::logic_error, "Error, dim=" << dim << " < 2 is not allowed!" );

    int    maxNumIters  = 10*dim;
    
    // Run using float
    result = runPowerMethodExample<float>(dim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;

    // Run using double
    result = runPowerMethodExample<double>(dim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;

#ifdef HAVE_TEUCHOS_COMPLEX
    
    // Run using std::complex<float>
    result = runPowerMethodExample<std::complex<float> >(dim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;
    
    // Run using std::complex<double>
    result = runPowerMethodExample<std::complex<double> >(dim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;

#endif // HAVE_TEUCHOS_COMPLEX

#ifdef HAVE_TEUCHOS_GNU_MP
    
    // Run using mpf_class
    result = runPowerMethodExample<mpf_class >(dim,maxNumIters,verbose,dumpAll);
    if(!result) success = false;

#ifdef HAVE_TEUCHOS_COMPLEX
    
    // Run using std::complex<mpf_class>
    //result = runPowerMethodExample<std::complex<mpf_class> >(dim,maxNumIters,verbose,dumpAll);
    //if(!result) success = false;
    //The above commented-out code throws a floating-point exception?
 
#endif // HAVE_TEUCHOS_COMPLEX


#endif
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)
  
  if (verbose) {
    if(success)  *out << "\nCongratulations! All of the tests checked out!\n";
    else         *out << "\nOh no! At least one of the tests failed!\n";
  }
  
  return success ? 0 : 1;

} // end main()
