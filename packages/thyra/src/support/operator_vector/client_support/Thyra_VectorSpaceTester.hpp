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

#ifndef THYRA_VECTOR_SPACE_TESTER_HPP
#define THYRA_VECTOR_SPACE_TESTER_HPP

#include "Thyra_VectorSpaceTesterDecl.hpp"
#include "Thyra_TestingTools.hpp"

namespace Thyra {

template <class Scalar>
VectorSpaceTester<Scalar>::VectorSpaceTester(
  const ScalarMag     warning_tol
  ,const ScalarMag    error_tol
  ,const int          num_random_vectors
  ,const int          num_mv_cols
  ,const bool         show_all_tests
  ,const bool         dump_all
  )
  :num_mv_cols_(num_mv_cols) // declared first!
  ,warning_tol_(warning_tol)
  ,error_tol_(error_tol)
  ,num_random_vectors_(num_random_vectors)
  ,show_all_tests_(show_all_tests)
  ,dump_all_(dump_all)
{
  vectorTester_.warning_tol(warning_tol);
  vectorTester_.error_tol(error_tol);
  vectorTester_.num_random_vectors(num_random_vectors);
  vectorTester_.show_all_tests(show_all_tests);
  vectorTester_.dump_all(dump_all);
}

template <class Scalar>
bool VectorSpaceTester<Scalar>::check(
  const VectorSpaceBase<Scalar>  &vs
  ,std::ostream                  *out
  ,const std::string             &leadingIndent
  ,const std::string             &indentSpacer
  ) const
{
  using std::endl;
  using Teuchos::describe;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;

  const std::string &li = leadingIndent, &is = indentSpacer;
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  bool result, success = true;

  if(out) *out <<endl<<li<< "*** Entering Thyra::VectorSpaceTester<"<<ST::name()<<">::check(vs,...) ...\n";

  if(out) *out <<endl<<li<< "Testing a vector space vs described as:\n" << describe(vs,verbLevel,li,is);

  if(out)
    *out
      <<endl<<li<< "A) Calling basic query functions ...\n"
      <<endl<<li<< "vs.dim() = " << vs.dim()
      <<endl<<li<< "vs.isInCore() = " << vs.isInCore() << std::endl;

  if(out) *out <<endl<<li<< "B) Checking that vs is compatible with itself ...\n";

  if(out) *out <<endl<<li<< "vs.isCompatible(vs)=";
  result = vs.isCompatible(vs);
  if(!result) success = false;
  if(out) *out << result << " == true : " << passfail(result) << std::endl;

  if(out) *out <<endl<<li<< "C) Creating a randomized vector member v ...\n";
  Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> >
    v = createMember(vs);
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),&*v);

  if(out) *out <<endl<<li<< "D) Testing the VectorBase interface of v ...\n";

  result = vectorTester_.check(*v,out,li+is,is);
  if(!result) success = false;

  if(out) *out <<endl<<li<< "C) Creating a randomized MultiVector member mv ...\n";
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<Scalar> >
    mv = createMembers(vs,num_mv_cols());
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),&*mv);

  if(out) *out <<endl<<li<< "D) Testing the MultiVectorBase interface of mv ...\n";

  result = vectorTester_.multiVectorTester().check(*mv,out,li+is,is);
  if(!result) success = false;

  if(out) *out <<endl<<li<< "*** Leaving Thyra::VectorSpaceTester<"<<ST::name()<<">::check(vs,...) ...\n";
  
  return success;

}

} // namespace Thyra

#endif // THYRA_VECTOR_SPACE_TESTER_HPP
