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

#ifndef THYRA_LINEAR_OP_TESTER_HPP
#define THYRA_LINEAR_OP_TESTER_HPP

#include "Thyra_LinearOpTesterDecl.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Teuchos_arrayArg.hpp"

namespace Thyra {

template<class Scalar>
LinearOpTester<Scalar>::LinearOpTester(
  const bool          check_linear_properties
  ,const ScalarMag    linear_properties_warning_tol
  ,const ScalarMag    linear_properties_error_tol
  ,const bool         check_adjoint
  ,const ScalarMag    adjoint_warning_tol
  ,const ScalarMag    adjoint_error_tol
  ,const bool         check_for_symmetry
  ,const ScalarMag    symmetry_warning_tol
  ,const ScalarMag    symmetry_error_tol
  ,const int          num_random_vectors
  ,const bool         show_all_tests
  ,const bool         dump_all
  )
  :check_linear_properties_(check_linear_properties)
  ,linear_properties_warning_tol_(linear_properties_warning_tol)
  ,linear_properties_error_tol_(linear_properties_error_tol)
  ,check_adjoint_(check_adjoint)
  ,adjoint_warning_tol_(adjoint_warning_tol)
  ,adjoint_error_tol_(adjoint_error_tol)
  ,check_for_symmetry_(check_for_symmetry)
  ,symmetry_warning_tol_(symmetry_warning_tol)
  ,symmetry_error_tol_(symmetry_error_tol)
  ,num_random_vectors_(num_random_vectors)
  ,show_all_tests_(show_all_tests)
  ,dump_all_(dump_all)
{}

template<class Scalar>
void LinearOpTester<Scalar>::set_all_warning_tol( const ScalarMag warning_tol )
{
  linear_properties_warning_tol_  = warning_tol;
  adjoint_warning_tol_            = warning_tol;
  symmetry_warning_tol_           = warning_tol;
}

template<class Scalar>
void LinearOpTester<Scalar>::set_all_error_tol( const ScalarMag error_tol )
{
  linear_properties_error_tol_  = error_tol;
  adjoint_error_tol_            = error_tol;
  symmetry_error_tol_           = error_tol;
}

template<class Scalar>
bool LinearOpTester<Scalar>::check(
  const LinearOpBase<Scalar>  &op
  ,std::ostream               *out
  ,const std::string          &leadingIndent
  ,const std::string          &indentSpacer
  ) const
{

  using std::endl;
  using Teuchos::arrayArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  bool success = true, result;
  const Scalar zero = ST::zero();
  const Scalar one = ST::one();
  const Scalar half = Scalar(0.5)*one;
  const std::string &li = leadingIndent, &is = indentSpacer;
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  // ToDo 04/28/2005:
  // * Test the MultiVectorBase apply() function and output to the VectorBase apply() function!

  if(out) {
    *out <<endl<<li<< "*** Entering LinearOpTester<"<<ST::name()<<">::check(op,...) ...\n";
    if(show_all_tests())
      *out <<endl<<li<< "describe op:\n" << Teuchos::describe(op,verbLevel,li,is);
    else
      *out <<endl<<li<< "describe op: " << op.describe() << endl;
  }

  if(out)
    *out <<endl<<li << "Checking the domain and range spaces ... ";

  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
    domain = op.domain(),
    range  = op.range();
  
  if(1) {

    std::ostringstream oss;
    bool these_results = true;

    oss <<endl<<li<< "op.domain().get() != NULL ? ";
    result = domain.get() != NULL;
    if(!result) these_results = false;
    oss << passfail(result) << endl;
    
    oss <<endl<<li<< "op.range().get() != NULL ? ";
    result = range.get() != NULL;
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }

  if( check_linear_properties() ) {

    if(out)	*out <<endl<<li<< "this->check_linear_properties()==true: Checking the linear properties of the forward linear operator ... ";

    std::ostringstream oss;
    bool these_results = true;

    oss <<endl<<li<< "op.opSupported(NOTRANS) == true ? ";
    result = op.opSupported(NOTRANS);
    if(!result) these_results = false;
    oss << passfail(result) << endl;
    
    oss
      <<endl<<li<< "Checking that the forward operator is truly linear:\n"
      <<endl<<li<< "  0.5*op*(v1 + v2) == 0.5*op*v1 + 0.5*op*v2"
      <<endl<<li<< "          \\_____/         \\___/"
      <<endl<<li<< "             v3            v5"
      <<endl<<li<< "  \\_____________/     \\___________________/"
      <<endl<<li<< "         v4                    v5"
      <<endl<<li<< ""
      <<endl<<li<< "           sum(v4) == sum(v5)"
      << endl;

    for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

      oss <<endl<<li<< "Random vector tests = " << rand_vec_i << endl;
      
      oss <<endl<<li<< "v1 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v1 = createMember(domain);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v1 );
      if(dump_all()) oss <<endl<<li<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
      oss <<endl<<li<< "v2 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v2 = createMember(domain);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v2 );
      if(dump_all()) oss <<endl<<li<< "v2 =\n" << describe(*v2,verbLevel,li,is);
      
      oss <<endl<<li<< "v3 = v1 + v2 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v3 = createMember(domain);
      linear_combination( 2, arrayArg<Scalar>(one,one)(), arrayArg<const VectorBase<Scalar>*>(&*v1,&*v2)(), zero, &*v3 );
      if(dump_all()) oss <<endl<<li<< "v3 =\n" << describe(*v3,verbLevel,li,is);

      oss <<endl<<li<< "v4 = 0.5*op*v3 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v4 = createMember(range);
      apply( op, NOTRANS, *v3, &*v4, half );
      if(dump_all()) oss <<endl<<li<< "v4 =\n" << describe(*v4,verbLevel,li,is);

      oss <<endl<<li<< "v5 = op*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v5 = createMember(range);
      apply( op, NOTRANS, *v1, &*v5 );
      if(dump_all()) oss <<endl<<li<< "v5 =\n" << describe(*v5,verbLevel,li,is);

      oss <<endl<<li<< "v5 = 0.5*op*v2 + 0.5*v5 ...\n" ;
      apply( op, NOTRANS, *v2, &*v5, half, half );
      if(dump_all()) oss <<endl<<li<< "v5 =\n" << describe(*v5,verbLevel,li,is);
      
      const Scalar
        sum_v4 = sum(*v4),
        sum_v5 = sum(*v5);

      result = testRelErr(
        "sum(v4)", sum_v4
        ,"sum(v5)", sum_v5
        ,"linear_properties_error_tol()", linear_properties_error_tol()
        ,"linear_properties_warning_tol()", linear_properties_warning_tol()
        ,&oss,li
        );
      if(!result) these_results = false;

    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }
  else {
    if(out) *out <<endl<<li<< "this->check_linear_properties()==false: Skipping the check of the linear properties of the forward operator!\n";
  }

  if( check_linear_properties() && check_adjoint() ) {

    if(out)	*out <<endl<<li<< "(this->check_linear_properties()&&this->check_adjoint())==true: Checking the linear properties of the adjoint operator ... ";

    std::ostringstream oss;
    bool these_results = true;

    oss <<endl<<li<< "op.opSupported(CONJTRANS) == true ? ";
    result = op.opSupported(CONJTRANS);
    if(!result) these_results = false;
    oss << passfail(result) << endl;
    
    oss
      <<endl<<li<< "Checking that the adjoint operator is truly linear:\n"
      <<endl<<li<< "  0.5*op'*(v1 + v2) == 0.5*op'*v1 + 0.5*op'*v2"
      <<endl<<li<< "           \\_____/         \\____/"
      <<endl<<li<< "              v3             v5"
      <<endl<<li<< "  \\_______________/    \\_____________________/"
      <<endl<<li<< "         v4                      v5"
      <<endl<<li<< ""
      <<endl<<li<< "           sum(v4) == sum(v5)"
      << endl;

    for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

      oss <<endl<<li<< "Random vector tests = " << rand_vec_i << endl;
      
      oss <<endl<<li<< "v1 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v1 = createMember(range);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v1 );
      if(dump_all()) oss <<endl<<li<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
      oss <<endl<<li<< "v2 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v2 = createMember(range);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v2 );
      if(dump_all()) oss <<endl<<li<< "v2 =\n" << describe(*v2,verbLevel,li,is);
      
      oss <<endl<<li<< "v3 = v1 + v2 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v3 = createMember(range);
      linear_combination( 2, arrayArg<Scalar>(one,one)(), arrayArg<const VectorBase<Scalar>*>(&*v1,&*v2)(), zero, &*v3 );
      if(dump_all()) oss <<endl<<li<< "v3 =\n" << describe(*v3,verbLevel,li,is);

      oss <<endl<<li<< "v4 = 0.5*op'*v3 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v4 = createMember(domain);
      apply( op, CONJTRANS, *v3, &*v4, half );
      if(dump_all()) oss <<endl<<li<< "v4 =\n" << describe(*v4,verbLevel,li,is);

      oss <<endl<<li<< "v5 = op'*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v5 = createMember(domain);
      apply( op, CONJTRANS, *v1, &*v5 );
      if(dump_all()) oss <<endl<<li<< "v5 =\n" << describe(*v5,verbLevel,li,is);

      oss <<endl<<li<< "v5 = 0.5*op'*v2 + 0.5*v5 ...\n" ;
      apply( op, CONJTRANS, *v2, &*v5, half, half );
      if(dump_all()) oss <<endl<<li<< "v5 =\n" << describe(*v5,verbLevel,li,is);
      
      const Scalar
        sum_v4 = sum(*v4),
        sum_v5 = sum(*v5);

      result = testRelErr(
        "sum(v4)", sum_v4
        ,"sum(v5)", sum_v5
        ,"linear_properties_error_tol()", linear_properties_error_tol()
        ,"linear_properties_warning_tol()", linear_properties_warning_tol()
        ,&oss,li
        );
      if(!result) these_results = false;

    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }
  else {
    if(out) *out <<endl<<li<< "(this->check_linear_properties()&&this->check_adjoint())==false: Skipping the check of the linear properties of the adjoint operator!\n";
  }
  
  if( check_adjoint() ) {

    if(out)	*out <<endl<<li<< "this->check_adjoint()==true: Checking the agreement of the adjoint and forward operators ... ";

    std::ostringstream oss;
    bool these_results = true;
    
    oss <<endl<<li<< "op.opSupported(CONJTRANS) == true ? ";
    result = op.opSupported(CONJTRANS);
    if(!result) these_results = false;
    oss << passfail(result) << endl;
    
    oss
      <<endl<<li<< "Checking that the adjoint agrees with the non-adjoint operator as:\n"
      <<endl<<li<< "  <0.5*op'*v2,v1> == <v2,0.5*op*v1>"
      <<endl<<li<< "   \\________/            \\_______/"
      <<endl<<li<< "       v4                   v3"
      <<endl<<li<< ""
      <<endl<<li<< "         <v4,v1>  == <v2,v3>"
      << endl;
    
    for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {
      
      oss <<endl<<li<< "Random vector tests = " << rand_vec_i << endl;
      
      oss <<endl<<li<< "v1 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v1 = createMember(domain);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v1 );
      if(dump_all()) oss <<endl<<li<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
      oss <<endl<<li<< "v2 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v2 = createMember(range);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v2 );
      if(dump_all()) oss <<endl<<li<< "v2 =\n" << describe(*v2,verbLevel,li,is);
      
      oss <<endl<<li<< "v3 = 0.5*op*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v3 = createMember(range);
      apply( op, NOTRANS, *v1, &*v3, half );
      if(dump_all()) oss <<endl<<li<< "v3 =\n" << describe(*v3,verbLevel,li,is);
      
      oss <<endl<<li<< "v4 = 0.5*op'*v2 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v4 = createMember(domain);
      apply( op, CONJTRANS, *v2, &*v4, half );
      if(dump_all()) oss <<endl<<li<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
      const Scalar
        prod1 = domain->scalarProd(*v4,*v1),
        prod2 = range->scalarProd(*v2,*v3);

      result = testRelErr(
        "<v4,v1>", prod1
        ,"<v2,v3>", prod2
        ,"adjoint_error_tol()", adjoint_error_tol()
        ,"adjoint_warning_tol()", adjoint_warning_tol()
        ,&oss,li
        );
      if(!result) these_results = false;

    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }
  else {
    if(out)	*out <<endl<<li<< "this->check_adjoint()==false: Skipping check for the agreement of the adjoint and forward operators!\n";
  }

  if( check_for_symmetry() ) {

    if(out) *out <<endl<<li<< "this->check_for_symmetry()==true: Performing check of symmetry ... ";

    std::ostringstream oss;
    bool these_results = true;

    oss <<endl<<li<< "op.domain()->isCompatible(*op.range()) == true : ";
    result = op.domain()->isCompatible(*op.range());
    if(!result) success = false;
    oss << passfail(result) << endl;

    oss
      <<endl<<li<< "Checking that the operator is symmetric as:\n"
      <<endl<<li<< "  <0.5*op*v2,v1> == <v2,0.5*op*v1>"
      <<endl<<li<< "   \\_______/            \\_______/"
      <<endl<<li<< "      v4                    v3"
      <<endl<<li<< ""
      <<endl<<li<< "         <v4,v1> == <v2,v3>"
      << endl;

    for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {
      
      oss <<endl<<li<< "Random vector tests = " << rand_vec_i << endl;
      
      if(dump_all()) oss <<endl<<li<< "v1 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v1 = createMember(domain);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v1 );
      if(dump_all()) oss <<endl<<li<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
      if(dump_all()) oss <<endl<<li<< "v2 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v2 = createMember(range);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v2 );
      if(dump_all()) oss <<endl<<li<< "v2 =\n" << describe(*v2,verbLevel,li,is);
      
      if(dump_all()) oss <<endl<<li<< "v3 = 0.5*op*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v3 = createMember(range);
      apply( op, NOTRANS, *v1, &*v3, half );
      if(dump_all()) oss <<endl<<li<< "v3 =\n" << describe(*v3,verbLevel,li,is);
      
      if(dump_all()) oss <<endl<<li<< "v4 = 0.5*op*v2 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v4 = createMember(domain);
      apply( op, NOTRANS, *v2, &*v4, half );
      if(dump_all()) oss <<endl<<li<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
      const Scalar
        prod1 = domain->scalarProd(*v4,*v1),
        prod2 = range->scalarProd(*v2,*v3);

      result = testRelErr(
        "<v4,v1>", prod1
        ,"<v2,v3>", prod2
        ,"symmetry_error_tol()", symmetry_error_tol()
        ,"symmetry_warning_tol()", symmetry_warning_tol()
        ,&oss,li
        );
      if(!result) success = false;

    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }
  else {
    if(out) *out <<endl<<li<< "this->check_for_symmetry()==false: Skipping check of symmetry ... ";
  }
  
  if(out)
    *out <<endl<<li<< "*** Leaving LinearOpTester<"<<ST::name()<<">::check(...)\n";

  return success;
}

template<class Scalar>
bool LinearOpTester<Scalar>::compare(
  const LinearOpBase<Scalar>  &op1
  ,const LinearOpBase<Scalar> &op2
  ,std::ostream               *out
  ,const std::string          &leadingIndent
  ,const std::string          &indentSpacer
  ) const
{

  using std::endl;
  using Teuchos::arrayArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  bool success = true, result;
  const Scalar one = ST::one();
  const Scalar half = Scalar(0.5)*one;
  const std::string &li = leadingIndent, &is = indentSpacer;
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  if(out) {
    *out
      <<endl<<li<< "*** Entering LinearOpTester<"<<ST::name()<<">::compare(op1,op2,...) ...\n";
    if(show_all_tests())
      *out <<endl<<li<< "describe op1:\n" << Teuchos::describe(op1,verbLevel,li,is);
    else
      *out <<endl<<li<< "describe op1: " << op1.describe() << endl;
    if(show_all_tests())
      *out <<endl<<li<< "describe op2:\n" << Teuchos::describe(op2,verbLevel,li,is);
    else
      *out <<endl<<li<< "describe op2: " << op2.describe() << endl;
  }

  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
    domain = op1.domain(),
    range  = op1.range();

  if(out) *out <<endl<<li<< "Checking that range and domain spaces are compatible ... ";

  if(1) {

    std::ostringstream oss;
    bool these_results = true;

    oss <<endl<<li<< "op1.domain()->isCompatible(*op2.domain()) ? ";
    result = op1.domain()->isCompatible(*op2.domain());
    if(!result) these_results = false;
    oss << passfail(result) << endl;
    
    oss <<endl<<li<< "op1.range()->isCompatible(*op2.range()) ? ";
    result = op1.range()->isCompatible(*op2.range());
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }

  if(out) *out <<endl<<li<< "Checking that op1 == op2 ... ";

  if(1) {

    std::ostringstream oss;
    bool these_results = true;

    oss
      <<endl<<li<< "Checking that op1 and op2 produce the same results:\n"
      <<endl<<li<< "  0.5*op1*v1 == 0.5*op2*v1"
      <<endl<<li<< "  \\________/    \\________/"
      <<endl<<li<< "      v2            v3"
      <<endl<<li<< ""
      <<endl<<li<< "   |sum(v2)| == |sum(v3)|"
      << endl;

    for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {
      
      oss <<endl<<li<< "Random vector tests = " << rand_vec_i << endl;
      
      if(dump_all()) oss <<endl<<li<< "v1 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v1 = createMember(domain);
      Thyra::randomize( Scalar(-one), Scalar(+one), &*v1 );
      if(dump_all()) oss <<endl<<li<< "v1 =\n" << *v1;
      
      if(dump_all()) oss <<endl<<li<< "v2 = 0.5*op1*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v2 = createMember(range);
      apply( op1, NOTRANS, *v1, &*v2, half );
      if(dump_all()) oss <<endl<<li<< "v2 =\n" << *v2;
      
      if(dump_all()) oss <<endl<<li<< "v3 = 0.5*op2*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<Scalar> >	v3 = createMember(range);
      apply( op2, NOTRANS, *v1, &*v3, half );
      if(dump_all()) oss <<endl<<li<< "v3 =\n" << *v3;
      
      const Scalar
        sum_v2 = sum(*v2),
        sum_v3 = sum(*v3);
      
      result = testRelErr(
        "sum(v2)", sum_v2
        ,"sum(v3)", sum_v3
        ,"linear_properties_error_tol()", linear_properties_error_tol()
        ,"linear_properties_warning_tol()", linear_properties_warning_tol()
        ,out,li
        );
      if(!result) these_results = false;
      
    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }

  if(out)
    *out <<endl<<li<< "*** Leaving LinearOpTester<"<<ST::name()<<">::compare(...)\n";

  return success;

}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_TESTER_HPP
