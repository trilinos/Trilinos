// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_VECTOR_TESTER_HPP
#define THYRA_VECTOR_TESTER_HPP

#include "Thyra_VectorTester_decl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_VerbosityLevel.hpp"


namespace Thyra {


template<class Scalar>
VectorTester<Scalar>::VectorTester(
  const ScalarMag     warning_tol_in
  ,const ScalarMag    error_tol_in
  ,const int          num_random_vectors_in
  ,const bool         show_all_tests_in
  ,const bool         dump_all_in
  )
  :warning_tol_(warning_tol_in)
  ,error_tol_(error_tol_in)
  ,num_random_vectors_(num_random_vectors_in)
  ,show_all_tests_(show_all_tests_in)
  ,dump_all_(dump_all_in)
{}


template<class Scalar>
bool VectorTester<Scalar>::check(
  const VectorBase<Scalar>       &v
  ,Teuchos::FancyOStream         *out_arg
  ) const
{

  using std::endl;
  using Teuchos::describe;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  //typedef typename ST::magnitudeType    ScalarMag;

  Teuchos::RCP<FancyOStream> out = Teuchos::rcp(out_arg,false);
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  OSTab tab(out,1,"THYRA");

  bool result, success = true;

  if(out.get()) *out <<endl<< "*** Entering Thyra::VectorTester<"<<ST::name()<<">::check(v,...) ...\n";

  if(out.get()) *out <<endl<< "Testing a VectorBase object described as:\n" << describe(v,verbLevel);

  if(out.get()) *out <<endl<< "A) Creating temporary vector t1, t2, t3, and t4 from v.space() ...\n";
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
    vs = v.space();
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    t1 = createMember(vs), t2 = createMember(vs), t3 = createMember(vs), t4 = createMember(vs);

  if(out.get()) *out <<endl<< "B) Testing VectorBase::applyOp(...) by calling a few standard RTOp operations ... ";

  const Scalar
    one   = ST::one(),
    two   = Scalar(2)*one,
    three = Scalar(3)*one;

  {
    using Teuchos::inoutArg;

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;
    
    *testOut <<endl<< "assign(t1.ptr(),2.0) ...\n";
    Thyra::assign( t1.ptr(), two );
    if(dump_all()) *testOut <<endl<< "\nt1 =\n" << describe(*t1,verbLevel);
    
    result = Teuchos::testRelErr<Scalar>(
      "sum(t1)", sum(*t1), "2*vs->dim()", two*Scalar(vs->dim()),
      "error_tol()", error_tol(), "warning_tol()", warning_tol(),
      inoutArg(*testOut)
      );
    if(!result) these_results = false;
    
    *testOut <<endl<< "assign(t2.ptr(),3.0) ...\n";
    Thyra::assign( t2.ptr(), three );
    if(dump_all()) *testOut <<endl<< "t2 =\n" << *t1;
    
    result = Teuchos::testRelErr<Scalar>(
      "sum(t2)",sum(*t2),"3*vs->dim()",three*Scalar(vs->dim()),
      "error_tol()",error_tol(),"warning_tol()",warning_tol(),
      inoutArg(*testOut)
      );
    if(!result) these_results = false;
    
    result = Teuchos::testRelErr<Scalar>(
      "vs->scalarProd(*t1,*t2)",vs->scalarProd(*t1,*t2),"2*3*vs->dim()",two*three*Scalar(vs->dim()),
      "error_tol()",error_tol(),"warning_tol()",warning_tol(),
      inoutArg(*testOut)
      );
    if(!result) these_results = false;

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  }
    
  // ToDo: Test the rest of the specific VectorBase interface on v1

  if(out.get()) *out <<endl<< "C) Checking the MultiVectorBase interface of v ...\n";
  result = multiVectorTester_.check(v, out.ptr());
  if(!result) success = false;

  if(out.get()) *out <<endl<< "*** Leaving Thyra::VectorTester<"<<ST::name()<<">::check(v,...) ...\n";
  
  return success;

}


} // namespace Thyra


#endif // THYRA_VECTOR_TESTER_HPP
