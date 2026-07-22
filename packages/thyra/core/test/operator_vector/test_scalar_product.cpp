// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"

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
  const int n,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &tol,
  const bool dumpAll,
  Teuchos::FancyOStream *out_arg
  )
{

  using Thyra::relErr;
  using Thyra::testBoolExpr;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  using Teuchos::Ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::as;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::OSTab;

  RCP<Teuchos::FancyOStream>
    out = rcp(new Teuchos::FancyOStream(rcp(out_arg,false)));

  if (nonnull(out)) *out << "\n*** Entering run_scalar_product_tests<"<<ST::name()<<">(...) ...\n" << std::boolalpha;

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
  if (nonnull(out) && dumpAll) *out << "\nop_coeff =\n" << *op_coeff;
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

  if (nonnull(out)) *out << "\nTesting LinearOpBase with Euclidean domain and range scalar products ...\n";
  {
    OSTab tab(out);
    Thyra::assign<Scalar>( op.ptr(), *op_coeff );
    if (nonnull(out) && dumpAll) *out << "\nop =\n" << *op;
    if (nonnull(out) && dumpAll) *out
      << "\nop' =\n"
      << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",
      op->domain()->isEuclidean(), true, out.ptr());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",
      op->range()->isEuclidean(), true, out.ptr());
    if(!result) success = false;
    result = linearOpTester.check(*op, out.ptr());
    if(!result) success = false;
  }
 
  if (nonnull(out)) *out << "\nTesting LinearOpBase with non-Euclidean domain and Euclidean range scalar products ...\n";
  {
    OSTab tab(out);
    range->setScalarProd(euclideanScalarProd);
    domain->setScalarProd(domainScalarProd);
    op->initialize(range,domain);
    Thyra::assign<Scalar>( op.ptr(), *op_coeff );
    if (nonnull(out) && dumpAll) *out << "\nop =\n" << *op;
    if (nonnull(out) && dumpAll) *out
      << "\nop' =\n"
      << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",
      op->domain()->isEuclidean(), false, out.ptr());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",
      op->range()->isEuclidean(), true, out.ptr());
    if(!result) success = false;
    result = linearOpTester.check(*op,out.ptr());
    if(!result) success = false;
  }
    
  if (nonnull(out)) *out << "\nTesting LinearOpBase with Euclidean domain and non-Euclidean range scalar products ...\n";
  {
    OSTab tab(out);
    range->setScalarProd(rangeScalarProd);
    domain->setScalarProd(euclideanScalarProd);
    op->initialize(range,domain);
    Thyra::assign<Scalar>( op.ptr(), *op_coeff );
    if (nonnull(out) && dumpAll) *out << "\nop =\n" << *op;
    if (nonnull(out) && dumpAll) *out
      << "\nop' =\n"
      << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",
      op->domain()->isEuclidean(), true, out.ptr());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",
      op->range()->isEuclidean(), false, out.ptr());
    if(!result) success = false;
    result = linearOpTester.check(*op,out.ptr());
    if(!result) success = false;
  }
    
  if (nonnull(out)) *out << "\nTesting LinearOpBase with non-Euclidean domain and non-Euclidean range scalar products ...\n";
  {
    OSTab tab(out);
    range->setScalarProd(rangeScalarProd);
    domain->setScalarProd(domainScalarProd);
    op->initialize(range,domain);
    Thyra::assign<Scalar>( op.ptr(), *op_coeff );
    if (nonnull(out) && dumpAll) *out << "\nop =\n" << *op;
    if (nonnull(out) && dumpAll) *out
      << "\nop' =\n"
      << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
    result = testBoolExpr("op->domain()->isEuclidean()",
      op->domain()->isEuclidean(), false, out.ptr());
    if(!result) success = false;
    result = testBoolExpr("op->range()->isEuclidean()",
      op->range()->isEuclidean(), false, out.ptr());
    if(!result) success = false;
    result = linearOpTester.check(*op,out.ptr());
    if(!result) success = false;
  }
    
  if (nonnull(out)) *out << "\n*** Leaving run_scalar_product_tests<"<<ST::name()<<">(...) ...\n";

  return success;

} // end run_scalar_product_tests() [Doxygen looks for this!]

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

    //
    // Run the tests
    //

#if defined(HAVE_TEUCHOS_INST_FLOAT) && defined(HAVE_TEUCHOS_BLASFLOAT)
    if( !run_scalar_product_tests<float>(n,float(1e-5),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
    if( !run_scalar_product_tests<double>(n,double(1e-14),dumpAll,verbose?&*out:NULL) ) success = false;
#if defined(HAVE_TEUCHOS_INST_COMPLEX_FLOAT) && defined(HAVE_TEUCHOS_INST_FLOAT) && defined(HAVE_TEUCHOS_BLASFLOAT)
    if( !run_scalar_product_tests<std::complex<float> >(n,float(1e-5),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#if defined(HAVE_TEUCHOS_INST_COMPLEX_DOUBLE)
    if( !run_scalar_product_tests<std::complex<double> >(n,double(1e-14),dumpAll,verbose?&*out:NULL) ) success = false;
#endif
#ifdef HAVE_TEUCHOS_GNU_MP
    if( !run_scalar_product_tests<mpf_class>(n,mpf_class(1e-14),dumpAll,verbose?&*out:NULL) ) success = false;
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

} // end main() [Doxygen looks for this!]
