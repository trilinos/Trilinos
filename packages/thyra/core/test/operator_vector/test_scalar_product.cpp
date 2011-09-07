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

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"

#ifndef SUN_CXX

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"
#include "Thyra_LinearOpScalarProd.hpp"
#include "Thyra_EuclideanScalarProd.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultSerialVectorSpaceConverter.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

/** \brief Main test driver function for scalar products
 */
template <class Scalar>
bool run_scalar_product_tests(
  const int                                                     n
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &tol
  ,const bool                                                   dumpAll
  ,Teuchos::FancyOStream                                        *out_arg
  )
{

  using Thyra::relErr;
  using Thyra::testBoolExpr;
  typedef Teuchos::ScalarTraits<Scalar>    ST;
  typedef typename ST::magnitudeType       ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  using Teuchos::Ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::as;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::OSTab;

  RCP<Teuchos::FancyOStream>
    out = rcp(new Teuchos::FancyOStream(rcp(out_arg,false)));

  if(out.get()) *out << "\n*** Entering run_scalar_product_tests<"<<ST::name()<<">(...) ...\n" << std::boolalpha;

  bool success = true, result;

  RCP<Thyra::DefaultSpmdVectorSpace<Scalar> >
    domain = Thyra::defaultSpmdVectorSpace<Scalar>(n/2),
    range  = Thyra::defaultSpmdVectorSpace<Scalar>(n);

  RCP<Thyra::DefaultSpmdMultiVector<Scalar> >
    op_coeff = rcp(new Thyra::DefaultSpmdMultiVector<Scalar>(range,domain)),
    op       = rcp(new Thyra::DefaultSpmdMultiVector<Scalar>(range,domain));

  RCP<Thyra::DefaultDiagonalLinearOp<Scalar> >
    domainScalarProdOp = rcp(
      new Thyra::DefaultDiagonalLinearOp<Scalar>(
        rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(domain)
        )
      ),
    rangeScalarProdOp = rcp(
      new Thyra::DefaultDiagonalLinearOp<Scalar>(
        rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(range)
        )
      );

  Thyra::seed_randomize<Scalar>(0);
  Thyra::randomize<Scalar>( as<Scalar>(-ST::one()), ST::one(),
    Ptr<Thyra::MultiVectorBase<Scalar> >(op_coeff.ptr()) );
  if(out.get() && dumpAll) *out << "\nop_coeff =\n" << *op_coeff;
  RCP<const Thyra::VectorSpaceConverterBase<ScalarMag,Scalar> >
    vecSpcConverterFromMag = rcp(new Thyra::DefaultSerialVectorSpaceConverter<ScalarMag,Scalar>());
  RCP<const Thyra::VectorSpaceBase<ScalarMag> >
    magDomain = vecSpcConverterFromMag->createVectorSpaceFrom(*rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(domain)),
    magRange  = vecSpcConverterFromMag->createVectorSpaceFrom(*rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(range));
  RCP<Thyra::VectorBase<ScalarMag> >
    _domainScalarProdDiag = createMember(magDomain),
    _rangeScalarProdDiag  = createMember(*magRange);
  Thyra::randomize( ScalarMag(ScalarMag(+1)*SMT::one()), ScalarMag(ScalarMag(+2)*SMT::one()), _domainScalarProdDiag.ptr() );
  Thyra::randomize( ScalarMag(ScalarMag(+1)*SMT::one()), ScalarMag(ScalarMag(+2)*SMT::one()), _rangeScalarProdDiag.ptr() );
  vecSpcConverterFromMag->convert( *_domainScalarProdDiag, &*domainScalarProdOp->getNonconstDiag() );
  vecSpcConverterFromMag->convert( *_rangeScalarProdDiag, &*rangeScalarProdOp->getNonconstDiag() );

  RCP<const Thyra::EuclideanScalarProd<Scalar> >
    euclideanScalarProd = rcp(new Thyra::EuclideanScalarProd<Scalar>());
  RCP<const Thyra::LinearOpScalarProd<Scalar> >
    domainScalarProd = rcp(new Thyra::LinearOpScalarProd<Scalar>(domainScalarProdOp.create_weak())),
    rangeScalarProd = rcp(new Thyra::LinearOpScalarProd<Scalar>(rangeScalarProdOp.create_weak()));

  const ScalarMag warning_tol = ScalarMag(1e-2)*tol, error_tol = tol;
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.linear_properties_warning_tol(warning_tol);
  linearOpTester.linear_properties_error_tol(error_tol);
  linearOpTester.adjoint_warning_tol(warning_tol);
  linearOpTester.adjoint_error_tol(error_tol);
  linearOpTester.show_all_tests(true);
  linearOpTester.dump_all(dumpAll);

  if(out.get()) *out << "\nTesting LinearOpBase with Euclidean domain and range scalar products ...\n";
  {
    OSTab tab(out);
    Thyra::assign( &*op, *op_coeff );
    if(out.get() && dumpAll) *out << "\nop =\n" << *op;
    if(out.get() && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",op->domain()->isEuclidean(),true,out.get());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",op->range()->isEuclidean(),true,out.get());
    if(!result) success = false;
    result = linearOpTester.check(*op,out.get());
    if(!result) success = false;
  }
 
  if(out.get()) *out << "\nTesting LinearOpBase with non-Euclidean domain and Euclidean range scalar products ...\n";
  {
    OSTab tab(out);
    range->setScalarProd(euclideanScalarProd);
    domain->setScalarProd(domainScalarProd);
    op->initialize(range,domain);
    Thyra::assign( &*op, *op_coeff );
    if(out.get() && dumpAll) *out << "\nop =\n" << *op;
    if(out.get() && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",op->domain()->isEuclidean(),false,out.get());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",op->range()->isEuclidean(),true,out.get());
    if(!result) success = false;
    result = linearOpTester.check(*op,out.get());
    if(!result) success = false;
  }
    
  if(out.get()) *out << "\nTesting LinearOpBase with Euclidean domain and non-Euclidean range scalar products ...\n";
  {
    OSTab tab(out);
    range->setScalarProd(rangeScalarProd);
    domain->setScalarProd(euclideanScalarProd);
    op->initialize(range,domain);
    Thyra::assign( &*op, *op_coeff );
    if(out.get() && dumpAll) *out << "\nop =\n" << *op;
    if(out.get() && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",op->domain()->isEuclidean(),true,out.get());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",op->range()->isEuclidean(),false,out.get());
    if(!result) success = false;
    result = linearOpTester.check(*op,out.get());
    if(!result) success = false;
  }
    
  if(out.get()) *out << "\nTesting LinearOpBase with non-Euclidean domain and non-Euclidean range scalar products ...\n";
  {
    OSTab tab(out);
    range->setScalarProd(rangeScalarProd);
    domain->setScalarProd(domainScalarProd);
    op->initialize(range,domain);
    Thyra::assign( &*op, *op_coeff );
    if(out.get() && dumpAll) *out << "\nop =\n" << *op;
    if(out.get() && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",op->domain()->isEuclidean(),false,out.get());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",op->range()->isEuclidean(),false,out.get());
    if(!result) success = false;
    result = linearOpTester.check(*op,out.get());
    if(!result) success = false;
  }
    
  if(out.get()) *out << "\n*** Leaving run_scalar_product_tests<"<<ST::name()<<">(...) ...\n";

  return success;

} // end run_scalar_product_tests() [Doxygen looks for this!]

#endif // SUN_CXX

int main( int argc, char* argv[] ) {

  bool success = true;
  bool verbose = true;

  using Teuchos::CommandLineProcessor;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from command-line
    //

    int n         = 4;
    bool dumpAll  = false;

    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "n", &n, "Number of elements in each constituent vector." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

#ifndef SUN_CXX

    //
    // Run the tests
    //

#ifdef HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_scalar_product_tests<float>(n,float(1e-5),dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_scalar_product_tests<double>(n,double(1e-14),dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_THYRA_COMPLEX)
#ifdef HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_scalar_product_tests<std::complex<float> >(n,float(1e-5),dumpAll,verbose?&*out:NULL) ) success = false;
#endif // HAVE_THYRA_TEUCHOS_BLASFLOAT
    if( !run_scalar_product_tests<std::complex<double> >(n,double(1e-14),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#ifdef HAVE_TEUCHOS_GNU_MP
    if( !run_scalar_product_tests<mpf_class>(n,mpf_class(1e-14),dumpAll,verbose?&*out:NULL) ) success = false;
#endif

#endif // ifndef SUN_CXX

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

#ifndef SUN_CXX

  if(verbose) {
    if(success)
      *out << "\nAll of the tests seem to have run successfully!\n";
    else
      *out << "\nOh no! at least one of the test failed!\n";	
  }
  
  return success ? 0 : 1;

#else // ifndef SUN_CXX

  if (verbose) {
    std::cout << "\nError, the test was never run since SUN_CXX was defined and this test does not build on the Sun compiler!\n";
  }
  
  return 1;

#endif //ifndef SUN_CXX

} // end main() [Doxygen looks for this!]
