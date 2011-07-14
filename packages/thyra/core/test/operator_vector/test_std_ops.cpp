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

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_VectorStdOpsTester.hpp"
#include "Thyra_MultiVectorStdOpsTester.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_DefaultComm.hpp"

namespace Thyra {

/** \brief Main test driver that runs tests on all standard operators
 */
template <class Scalar>
bool run_std_ops_tests(
  const int                                                       n
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType    max_rel_err
  ,const bool                                                     dumpAll
  ,Teuchos::FancyOStream                                          *out_arg
  )
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Thyra::Ordinal Index;

  RCP<Teuchos::FancyOStream>
    out = rcp(new Teuchos::FancyOStream(rcp(out_arg,false)));

  VectorStdOpsTester<Scalar> vectorStdOpsTester;
  vectorStdOpsTester.warning_tol(ScalarMag(0.1)*max_rel_err);
  vectorStdOpsTester.error_tol(max_rel_err);
  MultiVectorStdOpsTester<Scalar> multiVectorStdOpsTester;
  multiVectorStdOpsTester.warning_tol(ScalarMag(0.1)*max_rel_err);
  multiVectorStdOpsTester.error_tol(max_rel_err);

  if(out.get()) *out << "\n*** Entering run_std_ops_tests<"<<ST::name()<<">(...) ...\n";

  bool success = true;

  if(out.get()) *out << "\nCreating a serial vector space svs with n="<<n<<" vector elements ...\n";

  const RCP<const Teuchos::Comm<Ordinal> >
    comm = Teuchos::DefaultComm<Ordinal>::getComm();

  const RCP<Thyra::VectorSpaceBase<Scalar> > svs =
    Thyra::defaultSpmdVectorSpace<Scalar>(comm,n,-1);

  if(out.get()) *out << "\nTesting standard vector ops with svs ...\n";
  if(!vectorStdOpsTester.checkStdOps(*svs, OSTab(out).get(), dumpAll)) success = false;

  if(out.get()) *out << "\nTesting standard multi-vector ops with svs ...\n";
  if(!multiVectorStdOpsTester.checkStdOps(*svs, OSTab(out).get(), dumpAll)) success = false;

  const int numBlocks = 3;

  if(out.get()) *out << "\nCreating a product space pvs with numBlocks="<<numBlocks<<" and n="<<n<<"vector elements per block ...\n";

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > >
    vecSpaces(numBlocks);
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
    spaceBlock = Thyra::defaultSpmdVectorSpace<Scalar>(comm,n,-1);
  for( int i = 0; i < numBlocks; ++i )
    vecSpaces[i] = spaceBlock;

  RCP<Thyra::DefaultProductVectorSpace<Scalar> > pvs =
    rcp(new Thyra::DefaultProductVectorSpace<Scalar>(numBlocks,&vecSpaces[0]));

  if(out.get()) *out << "\nTesting standard vector ops with pvs ...\n";
  if(!vectorStdOpsTester.checkStdOps(*pvs, OSTab(out).get(), dumpAll)) success = false;

  if(out.get()) *out << "\nTesting standard multi-vector ops with pvs ...\n";
  if(!multiVectorStdOpsTester.checkStdOps(*pvs, OSTab(out).get(), dumpAll)) success = false;

  if(out.get()) *out << "\nCreating a nested product space ppvs with numBlocks="<<numBlocks<<" product spaces as components ...\n";

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > >
    blockVecSpaces(numBlocks);
  for( int i = 0; i < numBlocks; ++i )
    blockVecSpaces[i] = pvs;

  RCP<Thyra::DefaultProductVectorSpace<Scalar> > ppvs = 
    rcp(new Thyra::DefaultProductVectorSpace<Scalar>(numBlocks, &blockVecSpaces[0]));

  if(out.get()) *out << "\nTesting standard vector ops with ppvs ...\n";
  if(!vectorStdOpsTester.checkStdOps(*ppvs, OSTab(out).get(), dumpAll)) success = false;

  if(out.get()) *out << "\nTesting standard multi-vector ops with ppvs ...\n";
  if(!multiVectorStdOpsTester.checkStdOps(*ppvs, OSTab(out).get(), dumpAll)) success = false;

  return success;

}

} // namespace Thyra

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool verbose = true;
  bool dumpAll = false;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  // KL 27 Jul 2006 -- Commenting out unused variables 
  // const int procRank = Teuchos::GlobalMPISession::getRank();
  // const int numProc = Teuchos::GlobalMPISession::getNProc();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from the command-line
    //


    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);

    int local_dim = 4;
    clp.setOption( "local-dim", &local_dim, "Number of vector elements per process." );

    double eps_scale = 200.0;
    clp.setOption( "eps-scale", &eps_scale, "Constant (greater than 1) to scale eps by in error tests." );

    clp.setOption( "verbose", "quiet", &verbose,
      "Determines if any output is printed or not." );

    clp.setOption( "dump-all", "no-dump", &dumpAll, "Determines if quantities are dumped or not." );

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    //
    // Run the tests
    //

#if defined(HAVE_THYRA_FLOAT)
    if( !Thyra::run_std_ops_tests<float>(local_dim,float(eps_scale*Teuchos::ScalarTraits<float>::eps()),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
    if( !Thyra::run_std_ops_tests<double>(local_dim,double(eps_scale*Teuchos::ScalarTraits<double>::eps()),dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_THYRA_COMPLEX) && defined(HAVE_THYRA_FLOAT)
    if( !Thyra::run_std_ops_tests<std::complex<float> >(local_dim,float(eps_scale*Teuchos::ScalarTraits<float>::eps()),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#if defined(HAVE_THYRA_COMPLEX)
    if( !Thyra::run_std_ops_tests<std::complex<double> >(local_dim,double(eps_scale*Teuchos::ScalarTraits<double>::eps()),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#ifdef HAVE_TEUCHOS_GNU_MP
    //if( !Thyra::run_std_ops_tests<mpf_class>(local_dim,mpf_class(max_rel_err),dumpAll,verbose?&*out:NULL) ) success = false;
    // RAB: 4/16/2005: We can not instantiate the above since rmax() is not supported by this types ScalarTraits class
    // and it is needed by the class RTOpPack::ROpMaxIndexLessThanBound.  This can be fixed using a template
    // conditional but I have not done this yet.
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
      *out << "\nAll of the tests seem to have run successfully!\n";
    else
      *out << "\nOh no! at least one of the test failed!\n";	
  }
  
  return success ? 0 : 1;

}
