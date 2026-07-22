// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_VECTOR_SPACE_TESTER_HPP
#define THYRA_VECTOR_SPACE_TESTER_HPP

#include "Thyra_VectorSpaceTester_decl.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_TestingHelpers.hpp"


namespace Thyra {


template <class Scalar>
VectorSpaceTester<Scalar>::VectorSpaceTester(
  const ScalarMag     warning_tol_in
  ,const ScalarMag    error_tol_in
  ,const int          num_random_vectors_in
  ,const int          num_mv_cols_in
  ,const bool         show_all_tests_in
  ,const bool         dump_all_in
  )
  :num_mv_cols_(num_mv_cols_in) // declared first!
  ,warning_tol_(warning_tol_in)
  ,error_tol_(error_tol_in)
  ,num_random_vectors_(num_random_vectors_in)
  ,show_all_tests_(show_all_tests_in)
  ,dump_all_(dump_all_in)
{}


template <class Scalar>
bool VectorSpaceTester<Scalar>::check(
  const VectorSpaceBase<Scalar>  &vs
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

  if (nonnull(out)) *out <<endl<< "*** Entering Thyra::VectorSpaceTester<"<<ST::name()<<">::check(vs,...) ...\n";

  if (nonnull(out)) *out <<endl<< "Testing a vector space vs described as:\n" << describe(vs,verbLevel);

  if (nonnull(out))
    *out
      <<endl<< "A) Calling basic query functions ...\n"
      <<endl<< "vs.dim() = " << vs.dim()
      <<endl<< "vs.hasInCoreView() = " << vs.hasInCoreView() << std::endl;

  if (nonnull(out)) *out <<endl<< "B) Checking that vs is compatible with itself, cloning, etc. ...\n";

  TEUCHOS_TEST_ASSERT(vs.isCompatible(vs), *out, success);

  const RCP<const VectorSpaceBase<Scalar> > vs_clone = vs.clone();
  if (nonnull(vs_clone)) {
    TEUCHOS_TEST_ASSERT(vs.isCompatible(*vs_clone), *out, success);
  }

  if (nonnull(out)) *out <<endl<< "C) Creating a randomized vector member v ...\n";
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    v = createMember(vs);
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),v.ptr());

  if (nonnull(out)) *out <<endl<< "D) Testing the VectorBase interface of v ...\n";

  result = vectorTester_.check(*v,out.get());
  if(!result) success = false;

  if (nonnull(out)) *out <<endl<< "C) Creating a randomized MultiVector member mv ...\n";
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
    mv = createMembers(vs,num_mv_cols());
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),mv.ptr());

  if (nonnull(out)) *out <<endl<< "D) Testing the MultiVectorBase interface of mv ...\n";

  result = vectorTester_.multiVectorTester().check(*mv, out.ptr());
  if(!result) success = false;

  if (nonnull(out)) {
    if(success)
      *out << endl <<"Congratulations, this VectorSpaceBase object seems to check out!\n";
    else
      *out << endl <<"Oh no, at least one of the tests performed with this VectorSpaceBase object failed (see above failures)!\n";
    *out << endl << "*** Leaving VectorSpaceTester<"<<ST::name()<<">::check(vs,...)\n";
  }
  
  return success;

}


} // namespace Thyra


#endif // THYRA_VECTOR_SPACE_TESTER_HPP
