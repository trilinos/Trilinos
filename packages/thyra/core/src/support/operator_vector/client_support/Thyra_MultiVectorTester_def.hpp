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

#ifndef THYRA_MULTI_VECTOR_TESTER_DEF_HPP
#define THYRA_MULTI_VECTOR_TESTER_DEF_HPP

#include "Thyra_MultiVectorTester_decl.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_TestingHelpers.hpp"


namespace Thyra {


template<class Scalar>
MultiVectorTester<Scalar>::MultiVectorTester(
  const ScalarMag warning_tol_in,
  const ScalarMag error_tol_in,
  const int num_random_vectors_in,
  const bool show_all_tests_in,
  const bool dump_all_in
  )
  :warning_tol_(warning_tol_in),
   error_tol_(error_tol_in),
   num_random_vectors_(num_random_vectors_in),
   show_all_tests_(show_all_tests_in),
   dump_all_(dump_all_in)
{}


template<class Scalar>
bool MultiVectorTester<Scalar>::checkMultiVector(
  const VectorSpaceBase<Scalar> &vs,
  const Ptr<Teuchos::FancyOStream> &out_inout
  ) const
{

  using Teuchos::as;
  using Teuchos::describe;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::tuple;
  using Teuchos::null;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  //typedef typename ST::magnitudeType ScalarMag;

  RCP<FancyOStream> out;
  if (!is_null(out_inout))
    out = Teuchos::rcpFromPtr(out_inout);
  else
    out = Teuchos::fancyOStream(rcp(new Teuchos::oblackholestream));

  const Teuchos::EVerbosityLevel verbLevel =
    (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  OSTab tab(out,1,"THYRA");

  bool success = true;
  
  *out << "\n*** Entering "<<this->description()<<"::checkMultiVector(vs,...) ...\n";

  *out << "\nTesting MultiVectorBase objects created from vs = " << describe(vs, verbLevel);

  const Ordinal dim = vs.dim();
  const Scalar scalarDim = as<Scalar>(dim);

  int tc = 0;

  *out << "\n"<<tc<<") Checking non-contiguous non-const multi-vector views ...\n";
  ++tc;
  {
    OSTab tab2(out);
    const int numCols = 6;
    const RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, numCols);
    assign<Scalar>(mv.ptr(), ST::zero());
    const Scalar
      one = as<Scalar>(1.0),
      three = as<Scalar>(3.0),
      five = as<Scalar>(5.0);
    {
      const RCP<MultiVectorBase<Scalar> > mvView = mv->subView(tuple<int>(1, 3, 5)());
      assign<Scalar>(mvView->col(0).ptr(), one);
      assign<Scalar>(mvView->col(1).ptr(), three);
      assign<Scalar>(mvView->col(2).ptr(), five);
    }
    TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mv->col(0)), ST::zero(), error_tol_,
      *out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mv->col(1)), as<Scalar>(one*scalarDim), error_tol_,
      *out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mv->col(2)), ST::zero(), error_tol_,
      *out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mv->col(3)), as<Scalar>(three*scalarDim), error_tol_,
      *out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mv->col(4)), ST::zero(), error_tol_,
      *out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mv->col(5)), as<Scalar>(five*scalarDim), error_tol_,
      *out, success);
  }

  *out << "\n"<<tc<<") Checking non-contiguous const multi-vector views ...\n";
  ++tc;
  {
    OSTab tab2(out);
    const int numCols = 6;
    const RCP<MultiVectorBase<Scalar> > mv = createMembers(vs, numCols);
    const Scalar
      one = as<Scalar>(1.0),
      three = as<Scalar>(3.0),
      five = as<Scalar>(5.0);
    assign<Scalar>(mv.ptr(), ST::zero());
    assign<Scalar>(mv->col(1).ptr(), one);
    assign<Scalar>(mv->col(3).ptr(), three);
    assign<Scalar>(mv->col(5).ptr(), five);
    {
      const RCP<const MultiVectorBase<Scalar> > mvView =
        mv.getConst()->subView(tuple<int>(1, 3, 4, 5)());
      TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mvView->col(0)), as<Scalar>(one*scalarDim), error_tol_,
        *out, success);
      TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mvView->col(1)), as<Scalar>(three*scalarDim), error_tol_,
        *out, success);
      TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mvView->col(2)), ST::zero(), error_tol_,
        *out, success);
      TEUCHOS_TEST_FLOATING_EQUALITY( sum(*mvView->col(3)), as<Scalar>(five*scalarDim), error_tol_,
        *out, success);
    }
  }

  if(success)
    *out << "\nCongratulations, this MultiVectorBase objects"
         << " created form this vector space seems to check out!\n";
  else
    *out << "\nOh no, at least one of the tests performed failed!\n";

  *out << "\n*** Leaving "<<this->description()<<"::checkMultiVector(vs,...) ...\n";

  return success;

}


template<class Scalar>
bool MultiVectorTester<Scalar>::check(
  const MultiVectorBase<Scalar> &mv,
  const Ptr<Teuchos::FancyOStream> &out_inout
  ) const
{

  using Teuchos::describe;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  //typedef typename ST::magnitudeType ScalarMag;

  RCP<FancyOStream> out;
  if (!is_null(out_inout))
    out = Teuchos::rcpFromPtr(out_inout);
  else
    out = Teuchos::fancyOStream(rcp(new Teuchos::oblackholestream));

  const Teuchos::EVerbosityLevel verbLevel =
    (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  OSTab tab(out,1,"THYRA");

  bool result, success = true;

  *out << "\n*** Entering Thyra::MultiVectorTester<"<<ST::name()<<">::check(mv,...) ...\n";

  *out << "\nTesting a MultiVectorBase object mv described as:\n" << describe(mv,verbLevel);
  
  // ToDo: Test the specific MultiVectorBase interface
  
  *out << "\nChecking the LinearOpBase interface of mv ...\n";
  result =linearOpTester_.check(mv, out.ptr());
  if(!result) success = false;

  if(success)
    *out << "\nCongratulations, this MultiVectorBase object seems to check out!\n";
  else
    *out << "\nOh no, at least one of the tests performed with this MultiVectorBase object failed (see above failures)!\n";
  
  *out << "\n*** Leaving MultiVectorTester<"<<ST::name()<<">::check(mv,...)\n";

  return success;

}


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_TESTER_DEF_HPP
