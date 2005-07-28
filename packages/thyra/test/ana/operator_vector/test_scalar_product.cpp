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

#include "Thyra_SerialVectorSpaceStd.hpp"
#include "Thyra_SerialMultiVectorStd.hpp"
#include "Thyra_LinearOpScalarProd.hpp"
#include "Thyra_EuclideanScalarProd.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DiagonalLinearOp.hpp"
#include "Thyra_SerialVectorSpaceConverterStd.hpp"
#include "Thyra_ScaledAdjointLinearOp.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

/** \brief Main test driver function for scalar products
 */
template <class Scalar>
bool run_scalar_product_tests(
  const int                                                     n
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &tol
  ,const bool                                                   dumpAll
  ,std::ostream                                                 *out
  )
{

  using Thyra::relErr;
  typedef Teuchos::ScalarTraits<Scalar>    ST;
  typedef typename ST::magnitudeType       ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;

  if(out) *out << "\n*** Entering run_scalar_product_tests<"<<ST::name()<<">(...) ...\n";

  bool success = true, result;

  RefCountPtr<Thyra::ScalarProdVectorSpaceBase<Scalar> >
    domain = rcp(new Thyra::SerialVectorSpaceStd<Scalar>(n/2)),
    range  = rcp(new Thyra::SerialVectorSpaceStd<Scalar>(n));
  RefCountPtr<Thyra::SerialMultiVectorStd<Scalar> >
    op_coeff = rcp(new Thyra::SerialMultiVectorStd<Scalar>(range,domain)),
    op       = rcp(new Thyra::SerialMultiVectorStd<Scalar>(range,domain));
  Thyra::seed_randomize<Scalar>(0);
  Thyra::randomize( Scalar(Scalar(-1)*ST::one()), Scalar(Scalar(+1)*ST::one()), &*op_coeff );
  if(out && dumpAll) *out << "\nop_coeff =\n" << *op_coeff;
  RefCountPtr<const Thyra::VectorSpaceConverterBase<ScalarMag,Scalar> >
    vecSpcConverterFromMag = rcp(new Thyra::SerialVectorSpaceConverterStd<ScalarMag,Scalar>());
  RefCountPtr<const Thyra::VectorSpaceBase<ScalarMag> >
    magDomain = vecSpcConverterFromMag->createVectorSpaceFrom(*rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(domain)),
    magRange  = vecSpcConverterFromMag->createVectorSpaceFrom(*rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(range));
  RefCountPtr<Thyra::VectorBase<ScalarMag> >
    _domainScalarProdDiag = createMember(magDomain),
    _rangeScalarProdDiag  = createMember(*magRange);
  Thyra::randomize( ScalarMag(ScalarMag(+1)*SMT::one()), ScalarMag(ScalarMag(+2)*SMT::one()), &*_domainScalarProdDiag );
  Thyra::randomize( ScalarMag(ScalarMag(+1)*SMT::one()), ScalarMag(ScalarMag(+2)*SMT::one()), &*_rangeScalarProdDiag );
  RefCountPtr<Thyra::VectorBase<Scalar> >
    domainScalarProdDiag  = createMember(rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(domain)),
    rangeScalarProdDiag   = createMember(rcp_implicit_cast<const Thyra::VectorSpaceBase<Scalar> >(range));
  vecSpcConverterFromMag->convert( *_domainScalarProdDiag, &*domainScalarProdDiag );
  vecSpcConverterFromMag->convert( *_rangeScalarProdDiag, &*rangeScalarProdDiag );

  const ScalarMag warning_tol = ScalarMag(1e-2)*tol, error_tol = tol;
  Thyra::LinearOpTester<Scalar> linearOpTester;
  linearOpTester.linear_properties_warning_tol(warning_tol);
  linearOpTester.linear_properties_error_tol(error_tol);
  linearOpTester.adjoint_warning_tol(warning_tol);
  linearOpTester.adjoint_error_tol(error_tol);
  linearOpTester.show_all_tests(true);

  if(out) *out << "\nTesting LinearOpBase with Euclidean domain and range scalar products ...\n";
  Thyra::assign( &*op, *op_coeff );
  if(out && dumpAll) *out << "\nop =\n" << *op;
  if(out && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
  result = linearOpTester.check(*op,out);
  if(!result) success = false;
 
  if(out) *out << "\nTesting LinearOpBase with non-Euclidean domain and Euclidean range scalar products ...\n";
  range->setScalarProd(rcp(new Thyra::EuclideanScalarProd<Scalar>()));
  domain->setScalarProd(
    rcp(
      new Thyra::LinearOpScalarProd<Scalar>(
        rcp(new Thyra::DiagonalLinearOp<Scalar>(domainScalarProdDiag))
        )
      )
    );
  op->initialize(range,domain);
  Thyra::assign( &*op, *op_coeff );
  if(out && dumpAll) *out << "\nop =\n" << *op;
  if(out && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
  result = linearOpTester.check(*op,out);
  if(!result) success = false;
  
  if(out) *out << "\nTesting LinearOpBase with Euclidean domain and non-Euclidean range scalar products ...\n";
  range->setScalarProd(
    rcp(
      new Thyra::LinearOpScalarProd<Scalar>(
        rcp(new Thyra::DiagonalLinearOp<Scalar>(rangeScalarProdDiag))
        )
      )
    );
  domain->setScalarProd(rcp(new Thyra::EuclideanScalarProd<Scalar>()));
  op->initialize(range,domain);
  Thyra::assign( &*op, *op_coeff );
  if(out && dumpAll) *out << "\nop =\n" << *op;
  if(out && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
  result = linearOpTester.check(*op,out);
  if(!result) success = false;
  
  if(out) *out << "\nTesting LinearOpBase with non-Euclidean domain and non-Euclidean range scalar products ...\n";
  range->setScalarProd(
    rcp(
      new Thyra::LinearOpScalarProd<Scalar>(
        rcp(new Thyra::DiagonalLinearOp<Scalar>(rangeScalarProdDiag))
        )
      )
    );
  domain->setScalarProd(
    rcp(
      new Thyra::LinearOpScalarProd<Scalar>(
        rcp(new Thyra::DiagonalLinearOp<Scalar>(domainScalarProdDiag))
        )
      )
    );
  op->initialize(range,domain);
  Thyra::assign( &*op, *op_coeff );
  if(out && dumpAll) *out << "\nop =\n" << *op;
  if(out && dumpAll) *out << "\nop' =\n" << *Thyra::adjoint(Teuchos::rcp_implicit_cast<const Thyra::LinearOpBase<Scalar> >(op));
  result = linearOpTester.check(*op,out);
  if(!result) success = false;

  if(out) *out << "\n*** Leaving run_scalar_product_tests<"<<ST::name()<<">(...) ...\n";

  return success;

} // end run_scalar_product_tests() [Doxygen looks for this!]

#endif // __sun

int main( int argc, char* argv[] ) {

  bool success = true;
  bool verbose = true;

  using Teuchos::CommandLineProcessor;

  std::ostream &out = std::cout;

  try {

    //
    // Read options from command-line
    //

    int n         = 4;
    bool dumpAll  = false;

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not." );
    clp.setOption( "n", &n, "Number of elements in each constituent vector." );
    clp.setOption( "dump-all", "no-dump-all", &dumpAll, "Determines if vectors are printed or not." );
    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

#ifndef __sun

    //
    // Run the tests
    //

    if( !run_scalar_product_tests<float>(n,float(1e-5),dumpAll,verbose?&out:NULL) ) success = false;
    if( !run_scalar_product_tests<double>(n,double(1e-14),dumpAll,verbose?&out:NULL) ) success = false;
#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)
    if( !run_scalar_product_tests<std::complex<float> >(n,float(1e-5),dumpAll,verbose?&out:NULL) ) success = false;
    if( !run_scalar_product_tests<std::complex<double> >(n,double(1e-14),dumpAll,verbose?&out:NULL) ) success = false;
#endif
#ifdef HAVE_TEUCHOS_GNU_MP
    if( !run_scalar_product_tests<mpf_class>(n,mpf_class(1e-14),dumpAll,verbose?&out:NULL) ) success = false;
#endif

#endif // ifndef __sun

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

#ifndef __sun

  if(verbose) {
    if(success)
      out << "\nAll of the tests seem to have run successfully!\n";
    else
      out << "\nOh no! at least one of the test failed!\n";	
  }
  
  return success ? 0 : 1;

#else // ifndef __sun

  if (verbose) {
    std::cout << "\nError, the test was never run since __sun was defined and this test does not build on the Sun compiler!\n";
  }
  
  return 1;

#endif //ifndef __sun

} // end main() [Doxygen looks for this!]
