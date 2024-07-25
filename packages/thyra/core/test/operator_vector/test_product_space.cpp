// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_VectorStdOpsTester.hpp"
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

  Thyra::VectorSpaceTester<Scalar> vectorSpaceTester;
  vectorSpaceTester.warning_tol(ScalarMag(0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(showAllTests);
  vectorSpaceTester.dump_all(dumpAll);

  Thyra::VectorStdOpsTester<Scalar> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(ScalarMag(0.1)*tol);
  vectorStdOpsTester.error_tol(tol);

  Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > >
    vecSpaces(numBlocks);
  const RCP<const Teuchos::Comm<Thyra::Ordinal> >
    comm = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();
  const int numProcs = size(*comm);
  RCP<const Thyra::VectorSpaceBase<Scalar> >
    spaceBlock = Thyra::defaultSpmdVectorSpace<Scalar>(comm,n,-1);
  for( int i = 0; i < numBlocks; ++i )
    vecSpaces[i] = spaceBlock;
  
  if(out.get()) *out << "\nA) Performing basic tests on product vectors with SPMD constituent vectors ...\n";

  if(out.get()) *out << "\nCreating a product space ps with numBlocks="<<numBlocks<<" and n="<<n<<"vector elements per block ...\n";

  RCP<Thyra::DefaultProductVectorSpace<Scalar> > ps =
    Thyra::productVectorSpace<Scalar>(vecSpaces());

  if(out.get()) *out << "\nps->numBlocks()=";
  result = ps->numBlocks() == numBlocks;
  if(!result) success = false;
  if(out.get()) *out
    << ps->numBlocks() << " == numBlocks=" << numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;

  if(out.get()) *out << "\nTesting the product space ps ...\n";

  if(out.get()) *out << "\nps->dim()=";
  result = ps->dim() == numProcs*n*numBlocks;
  if(!result) success = false;
  if(out.get()) *out
    << ps->dim() << " == numProcs*n*numBlocks=" << numProcs*n*numBlocks
    << " : " << ( result ? "passed" : "failed" ) << std::endl;
  
  if(out.get()) *out << "\nTesting the VectorSpaceBase interface of ps ...\n";
  TEUCHOS_TEST_ASSERT(vectorSpaceTester.check(*ps, out.get()), *out, success);
  
  if(out.get()) *out << "\nTesting standard vector ops for ps ...\n";
  TEUCHOS_TEST_ASSERT(vectorStdOpsTester.checkStdOps(*ps, out.get()), *out, success);
  
  if(out.get()) *out << "\nB) Testing a nested product space of product vector spaces called pps ...\n";

  Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > >
    blockVecSpaces(numBlocks);
  for( int i = 0; i < numBlocks; ++i )
    blockVecSpaces[i] = ps;

  RCP<Thyra::DefaultProductVectorSpace<Scalar> > pps =
    Thyra::productVectorSpace<Scalar>(blockVecSpaces());
  
  if(out.get()) *out << "\nTesting the VectorSpaceBase interface of pps ...\n";
  TEUCHOS_TEST_ASSERT(vectorSpaceTester.check(*pps, out.get()), *out, success);
  
  if(out.get()) *out << "\nTesting standard vector ops for pps ...\n";
  TEUCHOS_TEST_ASSERT(vectorStdOpsTester.checkStdOps(*pps, out.get()), *out, success);
    
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

#if defined(HAVE_TEUCHOS_INST_FLOAT) && defined(HAVE_TEUCHOS_BLASFLOAT)
    if( !run_product_space_tests<float>(n,numBlocks,float(1e-4),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif
    if( !run_product_space_tests<double>(n,numBlocks,double(1e-13),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_TEUCHOS_INST_COMPLEX_FLOAT) && defined(HAVE_TEUCHOS_INST_FLOAT) && defined(HAVE_TEUCHOS_BLASFLOAT)
    if( !run_product_space_tests<std::complex<float> >(n,numBlocks,float(1e-4),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#if defined(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE)
    if( !run_product_space_tests<std::complex<double> >(n,numBlocks,double(1e-12),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#ifdef HAVE_TEUCHOS_GNU_MP
    //if( !run_product_space_tests<mpf_class>(n,numBlocks,mpf_class(1e-13),showAllTests,dumpAll,verbose?&*out:NULL) ) success = false;
    // Above commented out code will not compile because its ScalarTraits specialization does not support eps()
#endif

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
