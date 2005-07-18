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
#include "RealComplexFFTLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_ScalarTraits.hpp"

//
// Creates random complex multi-vectors for FFT with symmetic entries
//

template<class RealScalar>
class SymmetricComplexMultiVectorRandomizer : public Thyra::MultiVectorRandomizerBase< std::complex<RealScalar> > {
public:

  typedef std::complex<RealScalar> Scalar;

  bool isCompatible( const Thyra::VectorSpaceBase<Scalar> &space ) const
    {
      return space.isInCore();
    }

  void randomize( Thyra::MultiVectorBase<Scalar> *mv )
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
#     ifdef _DEBUG
      TEST_FOR_EXCEPT( mv == NULL );
      TEST_FOR_EXCEPT( mv->range()->dim() % 2 != 0 );
#     endif
      Thyra::ExplicitMutableMultiVectorView<Scalar> ev_mv(*mv);
      const Thyra::Index n = ev_mv.subDim();
      for( Thyra::Index j = 1; j <= ev_mv.numSubCols(); ++j ) {
        for( Thyra::Index i = 1; i <= n/2; ++i ) {
          const Scalar val = ST::random();
          ev_mv(i,j)     = val;
          ev_mv(n-i+1,j) = ST::conjugate(val);
        }
      }
    }
  
};

//
// This example program does some stuff with FFT.
//
template<class RealScalar>
bool run1DFFTExample(
  const int                                                      N
  ,const bool                                                    verbose
  ,const bool                                                    dumpAll
  ,const RealScalar                                              tolerance
  ,const int                                                     outputPrec
  )
{
  using Teuchos::RefCountPtr; using Teuchos::rcp;
  typedef std::complex<RealScalar> ComplexScalar;
  typedef Teuchos::ScalarTraits<RealScalar> ST;
  const std::string indentSpacer = "  ";
  bool success = true;
  bool result;

  if(outputPrec > 0) std::cout.precision(outputPrec);

  if(verbose)
    std::cout << "\n***\n*** Running 1D FFT example using real scalar type = \'" << ST::name() << "\' ...\n***\n";

  Teuchos::Time timer("");
  timer.start(true);

  if(verbose) std::cout << "\nConstructing a 1D complex-to-complex FFT linear operator C ...\n";

  Teuchos::RefCountPtr< const Thyra::LinearOpWithSolveBase<ComplexScalar> >
    C = Teuchos::rcp( new ComplexFFTLinearOp<RealScalar>(N) );
  if(verbose) std::cout << "\nTesting the LinearOpBase interface of the constructed linear operator C ...\n";
  Thyra::LinearOpTester<ComplexScalar> linearOpTester;
  linearOpTester.set_all_error_tol(tolerance);
  linearOpTester.set_all_warning_tol(RealScalar(RealScalar(1e-2)*tolerance));
  linearOpTester.show_all_tests(true);
  linearOpTester.dump_all(dumpAll);
  result = linearOpTester.check(*C,verbose?&std::cout:0,indentSpacer,indentSpacer);
  if(!result) success = false;

  if(verbose) std::cout << "\nTesting the LinearOpWithSolveBase interface of the constructed linear operator C ...\n";

  Thyra::LinearOpWithSolveTester<ComplexScalar> linearOpWithSolveTester;
  linearOpWithSolveTester.set_all_solve_tol(tolerance);
  linearOpWithSolveTester.set_all_slack_error_tol(RealScalar(RealScalar(1e+1)*tolerance));
  linearOpWithSolveTester.set_all_slack_warning_tol(tolerance);
  linearOpWithSolveTester.show_all_tests(true);
  linearOpWithSolveTester.dump_all(dumpAll);
  result = linearOpWithSolveTester.check(*C,verbose?&std::cout:0,indentSpacer,indentSpacer);
  if(!result) success = false;

  if(verbose) std::cout << "\nConstructing a 1D real-to-complex FFT linear operator R ...\n";

  Teuchos::RefCountPtr< const Thyra::LinearOpWithSolveBase< ComplexScalar, RealScalar > >
    R = Teuchos::rcp( new RealComplexFFTLinearOp<RealScalar>(N) );
  if(verbose) std::cout << "\nTesting the LinearOpBase interface of the constructed linear operator R ...\n";
  SymmetricComplexMultiVectorRandomizer<RealScalar> symmetricComplexMultiVectorRandomizer;
  Thyra::LinearOpTester<ComplexScalar,RealScalar> RlinearOpTester;
  RlinearOpTester.set_all_error_tol(tolerance);
  RlinearOpTester.set_all_warning_tol(RealScalar(RealScalar(1e-2)*tolerance));
  RlinearOpTester.show_all_tests(true);
  RlinearOpTester.dump_all(dumpAll);
  result = RlinearOpTester.check(*R,&symmetricComplexMultiVectorRandomizer,NULL,verbose?&std::cout:0,indentSpacer,indentSpacer);
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

    int    N             = 4;
    bool   dumpAll       = false;
    int    outputPrec    = -1;

    CommandLineProcessor  clp(false); // Don't throw exceptions

    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "N", &N, "Power of 2 for size of the FFT." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Print all objects or not." );
    clp.setOption( "output-prec", &outputPrec, "Precision for outputting floating point numbers." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    TEST_FOR_EXCEPTION( N < 0, std::logic_error, "Error, N=" << N << " < 1 is not allowed!" );

    // Run using float
    result = run1DFFTExample<float>(N,verbose,dumpAll,1e-5,outputPrec);
    if(!result) success = false;

    // Run using double
    result = run1DFFTExample<double>(N,verbose,dumpAll,1e-13,outputPrec);
    if(!result) success = false;

#ifdef HAVE_TEUCHOS_GNU_MP

    // Run using mpf_class
    //result = run1DFFTExample<mpf_class>(N,verbose,dumpAll,1e-13,outputPrec);
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
