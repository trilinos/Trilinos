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

#include "Thyra_ProductVectorSpace.hpp"
#include "Thyra_SerialVectorSpaceStd.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

/** \brief Main test driver function for composite product spaces
 */
template <class Scalar>
bool run_product_space_tests(
  const int                                                     n
  ,const int                                                    numBlocks
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &tol
  ,const bool                                                   showAllTests
  ,const bool                                                   dumpAll
  ,std::ostream                                                 *out
  )
{

  using Thyra::relErr;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;

  const std::string is = " ";

  if(out) *out << "\n*** Entering run_product_space_tests<"<<ST::name()<<">(...) ...\n";

  bool success = true, result;
  Scalar sresult1, sresult2;

  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.warning_tol(ScalarMag(0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);

  std::vector<Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> > >
    vecSpaces(numBlocks);
  Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<Scalar> >
    spaceBlock = Teuchos::rcp(new Thyra::SerialVectorSpaceStd<Scalar>(n));
  for( int i = 0; i < numBlocks; ++i )
    vecSpaces[i] = spaceBlock;

  if(out) *out << "\nA)Performing basic tests on product vectors with serial constituent vectors ...\n";

  if(out) *out << "\nCreating a product space ps with numBlocks="<<numBlocks<<" and n="<<n<<"vector elements per block ...\n";

  Thyra::ProductVectorSpace<Scalar> ps(numBlocks,&vecSpaces[0]);

  if(out) *out << "\nps.numBlocks()=";
  result = ps.numBlocks() == numBlocks;
  if(!result) success = false;
  if(out) *out
    << ps.numBlocks() << " == numBlocks=" << numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;

  if(out) *out << "\nTesting the product space ps ...\n";

  if(out) *out << "\nps.dim()=";
  result = ps.dim() == n*numBlocks;
  if(!result) success = false;
  if(out) *out
    << ps.dim() << " == n*numBlocks=" << n*numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;

  if(out) *out << "\nTesting the VectorSpaceBase interface of ps ...\n";
  result = vectorSpaceTester.check(ps,out,is,is);
  if(!result) success = false;

  if(out) *out
    << "\nB) Test the compatibility of serial vectors and product vectors with serial blocks."
    << "\n   These tests demonstrate the principle of how all in-core vectors are compatible ...\n";

  const Scalar
    one   = ST::one(),
    two   = Scalar(2)*one,
    three = Scalar(3)*one;

  if(out) *out << "\nCreating a serial vector space ss with numBlocks*n=" << numBlocks*n << " vector elements ...\n";

  Thyra::SerialVectorSpaceStd<Scalar> ss(numBlocks*n);

  if(out) *out << "\nTesting the serial space ss ...\n";
  result = vectorSpaceTester.check(ss,out,is,is);
  if(!result) success = false;

  if(out) *out << "\nCreating product vectors; pv1, pv2 ...\n";
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    pv1 = createMember(ps),
    pv2 = createMember(ps);

  if(out) *out << "\nassign(&pv1,2.0) ...\n";
  Thyra::assign( &*pv1, two );

  if(out) *out << "\nassign(&pv1,3.0) ...\n";
  Thyra::assign( &*pv2, three );

  if(out) *out << "\nCreating serial vectors; sv1, sv2 ...\n";
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    sv1 = createMember(ss),
    sv2 = createMember(ss);

  if(out) *out << "\nassign(&sv1,*pv1) ...\n";
  Thyra::assign( &*sv1, *pv1 );

  if(out) *out << "\nsum(sv1)=";
  sresult1 = Thyra::sum(*sv1);
  sresult2 = two*Scalar(ps.dim());
  result = ( ST::magnitude( Thyra::relErr( sresult1, sresult2 ) )
             < ST::magnitude( tol ) );
  if(!result) success = false;
  if(out) *out
    << sresult1 << " == 2*ps.dim()=" << sresult2
    << " : " << ( result ? "passed" : "failed" ) << std::endl;
  
  if(out && dumpAll) *out
    << "\nsv1 =\n" << *sv1;

  if(out) *out << "\nassign(&pv2,*sv1) ...\n";
  Thyra::assign( &*pv2, *sv1 );

  if(out) *out << "\nsum(pv2)=";
  sresult1 = Thyra::sum(*pv2);
  sresult2 = two*Scalar(ps.dim());
  result = ( ST::magnitude( Thyra::relErr( sresult1, sresult2 ) )
             < ST::magnitude( tol ) );
  if(!result) success = false;
  if(out) *out
    << sresult1 << " == 2*ps.dim()=" << sresult2
    << " : " << ( result ? "passed" : "failed" ) << std::endl;
  
  if(out && dumpAll) *out
    << "\npv2 =\n" << *pv2;

  // ToDo: Finish tests!

  if(out) *out
    << "\n*** Leaving run_product_space_tests<"<<ST::name()<<">(...) ...\n";

  return success;

} // end run_product_space_tests() [Doxygen looks for this!]

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;

  std::ostream &out = std::cout;

  try {

    //
    // Read options from command-line
    //

    int n              = 4;
    int numBlocks      = 4;
    bool showAllTests  = true;
    bool dumpAll       = false;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "n", &n, "Number of elements in each constituent vector." );
    clp.setOption( "num-blocks", &numBlocks, "blocks to create." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    clp.setOption( "show-all-tests", "no-show-all-tests", &showAllTests, "Determines if all tests are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Run the tests
    //

    if( !run_product_space_tests<float>(n,numBlocks,float(1e-5),showAllTests,dumpAll,verbose?&out:NULL) ) success = false;
    if( !run_product_space_tests<double>(n,numBlocks,double(1e-13),showAllTests,dumpAll,verbose?&out:NULL) ) success = false;
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
    if( !run_product_space_tests<std::complex<float> >(n,numBlocks,float(1e-5),showAllTests,dumpAll,verbose?&out:NULL) ) success = false;
    if( !run_product_space_tests<std::complex<double> >(n,numBlocks,double(1e-13),showAllTests,dumpAll,verbose?&out:NULL) ) success = false;
#endif
#ifdef HAVE_TEUCHOS_GNU_MP
    //if( !run_product_space_tests<mpf_class>(n,numBlocks,mpf_class(1e-13),showAllTests,dumpAll,verbose?&out:NULL) ) success = false;
    // Above commented out code will not compile because its ScalarTraits specialization does not support eps()
#endif

  } // end try
  catch( const std::exception &excpt ) {
    if(verbose)
      std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    if(verbose)
      std::cerr << "*** Caught an unknown exception!\n";
    success = false;
  }

  if(verbose) {
    if(success)
      out << "\nAll of the tests seem to have run successfully!\n";
    else
      out << "\nOh no! at least one of the test failed!\n";	
  }
  
  return success ? 0 : 1;

} // end main() [Doxygen looks for this!]
