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


namespace Thyra {


template<class Scalar>
MultiVectorTester<Scalar>::MultiVectorTester(
  const ScalarMag     warning_tol
  ,const ScalarMag    error_tol
  ,const int          num_random_vectors
  ,const bool         show_all_tests
  ,const bool         dump_all
  )
  :warning_tol_(warning_tol)
  ,error_tol_(error_tol)
  ,num_random_vectors_(num_random_vectors)
  ,show_all_tests_(show_all_tests)
  ,dump_all_(dump_all)
{
  linearOpTester_.set_all_warning_tol(warning_tol);
  linearOpTester_.set_all_error_tol(error_tol);
  linearOpTester_.num_random_vectors(num_random_vectors);
  linearOpTester_.show_all_tests(show_all_tests);
  linearOpTester_.dump_all(dump_all);
}


template<class Scalar>
bool MultiVectorTester<Scalar>::check(
  const MultiVectorBase<Scalar>  &mv
  ,Teuchos::FancyOStream         *out_arg
  ) const
{

  using std::endl;
  using Teuchos::describe;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType    ScalarMag;

  Teuchos::RCP<FancyOStream> out = Teuchos::rcp(out_arg,false);
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  OSTab tab(out,1,"THYRA");

  bool result, success = true;

  if(out.get()) *out <<endl<< "*** Entering Thyra::MultiVectorTester<"<<ST::name()<<">::check(mv,...) ...\n";

  if(out.get()) *out <<endl<< "Testing a MultiVectorBase object mv described as:\n" << describe(mv,verbLevel);
  
  // ToDo: Test the specific MultiVectorBase interface
  
  if(out.get()) *out <<endl<< "Checking the LinearOpBase interface of mv ...\n";
  result =linearOpTester_.check(mv,out.get());
  if(!result) success = false;

  if(out.get()) {
    if(success)
      *out << endl <<"Congratulations, this MultiVectorBase object seems to check out!\n";
    else
      *out << endl <<"Oh no, at least one of the tests performed with this MultiVectorBase object failed (see above failures)!\n";
    *out << endl << "*** Leaving MultiVectorTester<"<<ST::name()<<">::check(mv,...)\n";
  }
  
  return success;

}


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_TESTER_DEF_HPP
