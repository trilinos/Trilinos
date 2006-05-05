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

#include "Teuchos_CommandLineProcessor.hpp"

#ifndef __sun

#include "ComplexFFTLinearOp.hpp"
#include "RealComplexFFTLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_ListedMultiVectorRandomizer.hpp"
#include "Thyra_DefaultSerialVectorSpaceConverter.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_arrayArg.hpp"
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
      return space.hasInCoreView();
    }

  void randomize( Thyra::MultiVectorBase<Scalar> *mv )
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
#     ifdef _DEBUG
      TEST_FOR_EXCEPT( mv == NULL );
      TEST_FOR_EXCEPT( mv->range()->dim() % 2 != 0 );
#     endif
      Thyra::DetachedMultiVectorView<Scalar> ev_mv(*mv);
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
  using Teuchos::OSTab;
  typedef std::complex<RealScalar> ComplexScalar;
  typedef Teuchos::ScalarTraits<RealScalar> RST;
  bool success = true;
  bool result;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = ( verbose ? Teuchos::VerboseObjectBase::getDefaultOStream() : Teuchos::null );

  if(outputPrec > 0) out->precision(outputPrec);

  if(verbose)
    *out << "\n***\n*** Running 1D FFT example using real scalar type = \'" << RST::name() << "\' ...\n***\n";

  Teuchos::Time timer("");
  timer.start(true);

  if(verbose) *out << "\nConstructing a 1D complex-to-complex FFT linear operator C ...\n";

  Teuchos::RefCountPtr< const Thyra::LinearOpWithSolveBase<ComplexScalar> >
    C = Teuchos::rcp( new ComplexFFTLinearOp<RealScalar>(N) );

  if(verbose) *out << "\nConstructing as set of simple known vectors to be used as random domain and range vectors ...\n";
  Thyra::DefaultSerialVectorSpaceConverter<RealScalar,ComplexScalar>
    realToComplexConverter;
  RefCountPtr<const Thyra::VectorSpaceBase<RealScalar> >
    realDomainVecSpc = realToComplexConverter.createVectorSpaceFrom(*C->domain());
  RefCountPtr<Thyra::MultiVectorBase<RealScalar> >
    realDomainVec = Thyra::createMember(realDomainVecSpc);
  Thyra::seed_randomize<RealScalar>(0);
  Thyra::randomize( RealScalar(-RST::one()), RST::one(), &*realDomainVec );
  if(verbose && dumpAll)
    *out << "\nrealDomainVec:\n" << *realDomainVec;
  RefCountPtr<Thyra::MultiVectorBase<ComplexScalar> >
    complexDomainVec = Thyra::createMember(C->domain()),
    complexRangeVec = Thyra::createMember(C->range());
  realToComplexConverter.convert(*realDomainVec,&*complexDomainVec);
  Thyra::apply( *C, Thyra::NOTRANS, *complexDomainVec, &*complexRangeVec );
  if(verbose && dumpAll)
    *out << "\ncomplexDomainVec:\n" << *complexDomainVec << "\ncomplexRangeVec:\n" << *complexRangeVec;
  Thyra::ListedMultiVectorRandomizer<RealScalar>
    realDomainRand( Teuchos::arrayArg<RefCountPtr<const Thyra::MultiVectorBase<RealScalar> > >(realDomainVec)(), 1 );
  Thyra::ListedMultiVectorRandomizer<ComplexScalar>
    complexDomainRand( Teuchos::arrayArg<RefCountPtr<const Thyra::MultiVectorBase<ComplexScalar> > >(complexDomainVec)(), 1 ),
    complexRangeRand( Teuchos::arrayArg<RefCountPtr<const Thyra::MultiVectorBase<ComplexScalar> > >(complexRangeVec)(), 1 );

  if(verbose) *out << "\nTesting the LinearOpBase interface of the constructed linear operator C ...\n";
  Thyra::LinearOpTester<ComplexScalar> linearOpTester;
  linearOpTester.set_all_error_tol(tolerance);
  linearOpTester.set_all_warning_tol(RealScalar(RealScalar(1e-2)*tolerance));
  linearOpTester.show_all_tests(true);
  linearOpTester.dump_all(dumpAll);
  result = linearOpTester.check(*C,&complexRangeRand,&complexDomainRand,out.get());
  if(!result) success = false;

  if(verbose) *out << "\nTesting the LinearOpWithSolveBase interface of the constructed linear operator C ...\n";

  Thyra::LinearOpWithSolveTester<ComplexScalar> linearOpWithSolveTester;
  linearOpWithSolveTester.set_all_solve_tol(tolerance);
  linearOpWithSolveTester.set_all_slack_error_tol(RealScalar(RealScalar(1e+1)*tolerance));
  linearOpWithSolveTester.set_all_slack_warning_tol(tolerance);
  linearOpWithSolveTester.show_all_tests(true);
  linearOpWithSolveTester.dump_all(dumpAll);
  result = linearOpWithSolveTester.check(*C,out.get());
  if(!result) success = false;

  if(verbose) *out << "\nConstructing a 1D real-to-complex FFT linear operator R ...\n";

  Teuchos::RefCountPtr< const Thyra::LinearOpWithSolveBase< ComplexScalar, RealScalar > >
    R = Teuchos::rcp( new RealComplexFFTLinearOp<RealScalar>(N) );

  if(verbose) *out << "\nTesting the LinearOpBase interface of the constructed linear operator R ...\n";
  SymmetricComplexMultiVectorRandomizer<RealScalar> symmetricComplexMultiVectorRandomizer;
  Thyra::LinearOpTester<ComplexScalar,RealScalar> RlinearOpTester;
  RlinearOpTester.set_all_error_tol(tolerance);
  RlinearOpTester.set_all_warning_tol(RealScalar(RealScalar(1e-2)*tolerance));
  RlinearOpTester.show_all_tests(true);
  RlinearOpTester.dump_all(dumpAll);
  result = RlinearOpTester.check(*R,&complexRangeRand,&realDomainRand,out.get());
  if(!result) success = false;

  timer.stop();

  if(verbose) *out << "\nTotal time = " << timer.totalElapsedTime() << " sec\n";

  return success;

} // end run1DFFTExample()

#endif // __sun

//
// Actual main driver program
//
int main(int argc, char *argv[])
{

  bool success = true;

  using Teuchos::CommandLineProcessor;

  bool verbose = true;
  bool result;

  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read in command-line options
    //

    int    N             = 4;
    bool   dumpAll       = false;
    int    outputPrec    = -1;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    clp.setOption( "verbose", "quiet", &verbose, "Determines if any output is printed or not." );
    clp.setOption( "N", &N, "Power of 2 for size of the FFT." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Print all objects or not." );
    clp.setOption( "output-prec", &outputPrec, "Precision for outputting floating point numbers." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

#ifndef __sun

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

#endif // ifndef __sun

  }
  catch( const std::exception &excpt ) {
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    std::cerr << "*** Caught an unknown exception\n";
    success = false;
  }

#ifndef __sun

  if (verbose) {
    if(success)   *out << "\nCongratulations! All of the tests checked out!\n";
    else          *out << "\nOh no! At least one of the tests failed!\n";
  }
  
  return success ? 0 : 1;

#else // ifndef __sun

  if (verbose) {
    *out << "\nError, the test was never run since __sun was defined and this test does not build on the Sun compiler!\n";
  }
  
  return 1;

#endif // __sun

} // end main()
