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

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_VerboseObject.hpp"

/** \brief Main test driver function for composite product spaces
 */
template <class Scalar>
bool run_product_space_tests(
  const int                                                     n
  ,const int                                                    numBlocks
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &tol
  ,const bool                                                   showAllTests
  ,const bool                                                   dumpAll
  ,Teuchos::FancyOStream                                        *out_arg
  )
{

  using Thyra::relErr;
  using Teuchos::OSTab;
  using Teuchos::rcp;
  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::fancyOStream(rcp(out_arg,false));

  if(out.get()) *out << "\n*** Entering run_product_space_tests<"<<ST::name()<<">(...) ...\n";

  bool success = true, result;
  Scalar sresult1, sresult2;

  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.warning_tol(ScalarMag(0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);

  Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > >
    vecSpaces(numBlocks);
  const RCP<const Teuchos::Comm<Thyra::Index> >
    comm = Teuchos::DefaultComm<Thyra::Index>::getComm();
  const int numProcs = size(*comm);
  RCP<const Thyra::VectorSpaceBase<Scalar> >
    spaceBlock = Thyra::defaultSpmdVectorSpace<Scalar>(comm,n,-1);
  for( int i = 0; i < numBlocks; ++i )
    vecSpaces[i] = spaceBlock;
  
  if(out.get()) *out << "\nA) Performing basic tests on product vectors with SPMD constituent vectors ...\n";

  if(out.get()) *out << "\nCreating a product space ps with numBlocks="<<numBlocks<<" and n="<<n<<"vector elements per block ...\n";

  Thyra::DefaultProductVectorSpace<Scalar> ps(numBlocks,&vecSpaces[0]);

  if(out.get()) *out << "\nps.numBlocks()=";
  result = ps.numBlocks() == numBlocks;
  if(!result) success = false;
  if(out.get()) *out
    << ps.numBlocks() << " == numBlocks=" << numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;

  if(out.get()) *out << "\nTesting the product space ps ...\n";

  if(out.get()) *out << "\nps.dim()=";
  result = ps.dim() == numProcs*n*numBlocks;
  if(!result) success = false;
  if(out.get()) *out
    << ps.dim() << " == numProcs*n*numBlocks=" << numProcs*n*numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;
  
  if(out.get()) *out << "\nTesting the VectorSpaceBase interface of ps ...\n";
  result = vectorSpaceTester.check(ps,out.get());
  if(!result) success = false;
  
  if(out.get()) *out << "\nB) Testing a nested product space of product vector spaces called pps ...\n";

  Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > >
    blockVecSpaces(numBlocks);
  for( int i = 0; i < numBlocks; ++i )
    blockVecSpaces[i] = Teuchos::rcp(&ps,false);

  Thyra::DefaultProductVectorSpace<Scalar> pps(numBlocks,&blockVecSpaces[0]);
  
  if(out.get()) *out << "\nTesting the VectorSpaceBase interface of pps ...\n";
  result = vectorSpaceTester.check(pps,out.get());
  if(!result) success = false;

  if(numProcs==1) {

    if(out.get()) *out
      << "\nC) Test the compatibility of serial vectors and product vectors with serial blocks."
      << "\n   These tests demonstrate the principle of how all in-core vectors are compatible ...\n";
    
    const Scalar
      one   = ST::one(),
      two   = Scalar(2)*one,
      three = Scalar(3)*one;
    
    if(out.get()) *out << "\nCreating a serial vector space ss with numBlocks*n=" << numBlocks*n << " vector elements ...\n";
    
    const RCP<const Thyra::VectorSpaceBase<Scalar> > ss
      = Thyra::defaultSpmdVectorSpace<Scalar>(numBlocks*n);
    
    if(out.get()) *out << "\nTesting the serial space ss ...\n";
    result = vectorSpaceTester.check(*ss, out.get());
    if(!result) success = false;
    
    if(out.get()) *out << "\nCreating product vectors; pv1, pv2 ...\n";
    RCP<Thyra::VectorBase<Scalar> >
      pv1 = createMember(ps),
      pv2 = createMember(ps);
    
    if(out.get()) *out << "\nassign(&pv1,2.0) ...\n";
    Thyra::assign( &*pv1, two );
    
    if(out.get()) *out << "\nassign(&pv1,3.0) ...\n";
    Thyra::assign( &*pv2, three );
    
    if(out.get()) *out << "\nCreating serial vectors; sv1, sv2 ...\n";
    RCP<Thyra::VectorBase<Scalar> >
      sv1 = createMember(ss),
      sv2 = createMember(ss);
    
    if(out.get()) *out << "\nassign(&sv1,*pv1) ...\n";
    Thyra::assign( &*sv1, *pv1 );
    
    if(out.get()) *out << "\nsum(sv1)=";
    sresult1 = Thyra::sum(*sv1);
    sresult2 = two*Scalar(ps.dim());
    result = ( ST::magnitude( Thyra::relErr( sresult1, sresult2 ) )
               < ST::magnitude( tol ) );
    if(!result) success = false;
    if(out.get()) *out
      << sresult1 << " == 2*ps.dim()=" << sresult2
      << " : " << ( result ? "passed" : "failed" ) << std::endl;
    
    if(out.get() && dumpAll) *out
      << "\nsv1 =\n" << *sv1;
    
    if(out.get()) *out << "\nassign(&pv2,*sv1) ...\n";
    Thyra::assign( &*pv2, *sv1 );
    
    if(out.get()) *out << "\nsum(pv2)=";
    sresult1 = Thyra::sum(*pv2);
    sresult2 = two*Scalar(ps.dim());
    result = ( ST::magnitude( Thyra::relErr( sresult1, sresult2 ) )
               < ST::magnitude( tol ) );
    if(!result) success = false;
    if(out.get()) *out
      << sresult1 << " == 2*ps.dim()=" << sresult2
      << " : " << ( result ? "passed" : "failed" ) << std::endl;
    
    if(out.get() && dumpAll) *out
      << "\npv2 =\n" << *pv2;
    
    // ToDo: Finish tests!

  }
    
  if(out.get()) *out
    << "\n*** Leaving run_product_space_tests<"<<ST::name()<<">(...) ...\n";
  
  return success;

} // end run_product_space_tests() [Doxygen looks for this!]

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
  using Teuchos::RCP;

  bool success = true;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //

    int n              = 4;
    int numBlocks      = 4;
    bool showAllTests  = true;
    bool dumpAll       = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
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

#ifdef HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<float>(n,numBlocks,float(1e-5),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<double>(n,numBlocks,double(1e-13),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_TEUCHOS_COMPLEX)
#ifdef THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<std::complex<float> >(n,numBlocks,float(1e-5),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_product_space_tests<std::complex<double> >(n,numBlocks,double(1e-13),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif // defined(HAVE_TEUCHOS_COMPLEX)
#ifdef HAVE_TEUCHOS_GNU_MP
    //if( !run_product_space_tests<mpf_class>(n,numBlocks,mpf_class(1e-13),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
    // Above commented out code will not compile because its ScalarTraits specialization does not support eps()
#endif // HAVE_TEUCHOS_GNU_MP

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,success)

  if(verbose) {
    if(success)
      *out << "\nAll of the tests seem to have run successfully!\n";
    else
      *out << "\nOh no! at least one of the test failed!\n";	
  }
  
  return success ? 0 : 1;

} // end main() [Doxygen looks for this!]
