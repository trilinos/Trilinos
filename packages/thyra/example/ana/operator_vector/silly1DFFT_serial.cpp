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

#include "ComplexFFTLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_ScalarTraits.hpp"

//
// This example program does some stuff with FFT.
//
template<class RealScalar>
bool run1DFFTExample(
  const int                                                      N
  ,const bool                                                    verbose
  ,const RealScalar                                              tolerance
  )
{
  using Teuchos::RefCountPtr; using Teuchos::rcp;
  typedef Teuchos::ScalarTraits<RealScalar> ST;
  bool success = true;
  bool result;
  if(verbose)
    std::cout << "\n***\n*** Running 1D FFT example using real scalar type = \'" << ST::name() << "\' ...\n***\n";
  Teuchos::Time timer("");
  timer.start(true);
  if(verbose) std::cout << "\nConstructing a 1D FFT linear operator A ...\n";
  Teuchos::RefCountPtr< const Thyra::LinearOpWithSolveBase< std::complex<RealScalar> > >
    A = Teuchos::rcp( new ComplexFFTLinearOp<RealScalar>(N) );
  if(verbose) std::cout << "\nTesting the constructed linear operator A ...\n";
  Thyra::LinearOpTester< std::complex<RealScalar> > linearOpTester;
  linearOpTester.set_all_error_tol(tolerance);
  linearOpTester.set_all_warning_tol(RealScalar(RealScalar(1e-2)*tolerance));
  linearOpTester.show_all_tests(true);
  result = linearOpTester.check(*A,verbose?&std::cout:0);
  if(!result) success = false;
  timer.stop();
  if(verbose) std::cout << "\nTotal time = " << timer.totalElapsedTime() << " sec\n";

  return success;

} // end run1DFFTExample()

//
// Actual main driver program
//
int main(int argc, char *argv[])
{

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;
  bool result;

  try {

    //
    // Read in command-line options
    //

    int    N           = 4;

    CommandLineProcessor  clp(false); // Don't throw exceptions

    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "N", &N, "Power of 2 for size of the FFT." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPTION( N < 1, std::logic_error, "Error, N=" << N << " < 1 is not allowed!" );

    // Run using float
    result = run1DFFTExample<float>(N,verbose,1e-5);
    if(!result) success = false;

    // Run using double
    result = run1DFFTExample<double>(N,verbose,1e-13);
    if(!result) success = false;

#ifdef HAVE_TEUCHOS_GNU_MP

    // Run using mpf_class
    //result = run1DFFTExample<mpf_class>(N,verbose,1e-13);
    //if(!result) success = false;

#endif // HAVE_TEUCHOS_GNU_MP

  }
  catch( const std::exception &excpt ) {
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    std::cerr << "*** Caught an unknown exception\n";
    success = false;
  }

  if (verbose) {
    if(success)   std::cout << "\nCongratulations! All of the tests checked out!\n";
    else          std::cout << "\nOh no! At least one of the tests failed!\n";
  }
  
  return success ? 0 : 1;

} // end main()
