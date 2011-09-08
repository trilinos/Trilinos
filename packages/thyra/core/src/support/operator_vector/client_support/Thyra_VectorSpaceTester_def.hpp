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

#ifndef THYRA_VECTOR_SPACE_TESTER_HPP
#define THYRA_VECTOR_SPACE_TESTER_HPP

#include "Thyra_VectorSpaceTester_decl.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"


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

  if(out.get()) *out <<endl<< "*** Entering Thyra::VectorSpaceTester<"<<ST::name()<<">::check(vs,...) ...\n";

  if(out.get()) *out <<endl<< "Testing a vector space vs described as:\n" << describe(vs,verbLevel);

  if(out.get())
    *out
      <<endl<< "A) Calling basic query functions ...\n"
      <<endl<< "vs.dim() = " << vs.dim()
      <<endl<< "vs.hasInCoreView() = " << vs.hasInCoreView() << std::endl;

  if(out.get()) *out <<endl<< "B) Checking that vs is compatible with itself ...\n";

  if(out.get()) *out <<endl<< "vs.isCompatible(vs)=";
  result = vs.isCompatible(vs);
  if(!result) success = false;
  if(out.get()) *out << result << " == true : " << passfail(result) << std::endl;

  if(out.get()) *out <<endl<< "C) Creating a randomized vector member v ...\n";
  Teuchos::RCP<Thyra::VectorBase<Scalar> >
    v = createMember(vs);
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),v.ptr());

  if(out.get()) *out <<endl<< "D) Testing the VectorBase interface of v ...\n";

  result = vectorTester_.check(*v,out.get());
  if(!result) success = false;

  if(out.get()) *out <<endl<< "C) Creating a randomized MultiVector member mv ...\n";
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >
    mv = createMembers(vs,num_mv_cols());
  randomize(Scalar(-ST::one()),Scalar(+ST::one()),mv.ptr());

  if(out.get()) *out <<endl<< "D) Testing the MultiVectorBase interface of mv ...\n";

  result = vectorTester_.multiVectorTester().check(*mv, out.ptr());
  if(!result) success = false;

  if(out.get()) {
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
