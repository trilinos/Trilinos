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
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_ScaledAdjointLinearOp.hpp"
#include "Thyra_describeLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_UniversalMultiVectorRandomizer.hpp"

namespace Thyra {

// SymmetricLinearOpTester (using partial specialization only test symmetry on operators where RangeScalar and DomainScalar are the same)

template<class RangeScalar, class DomainScalar>
class SymmetricLinearOpTester {
public:
  typedef typename Teuchos::PromotionTraits<RangeScalar,DomainScalar>::promote Scalar;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  static void checkSymmetry(
    const LinearOpBase<RangeScalar,DomainScalar>  &op
    ,MultiVectorRandomizerBase<DomainScalar>      *dRand
    ,std::ostream                                 &oss
    ,const std::string                            &li
    ,const std::string                            &is
    ,const int                                    num_random_vectors
    ,const Teuchos::EVerbosityLevel               verbLevel
    ,const bool                                   dump_all
    ,const ScalarMag                              &symmetry_error_tol
    ,const ScalarMag                              &symmetry_warning_tol
    ,bool                                         *these_results
    )
    {
      using std::endl;
      typedef Teuchos::ScalarTraits<RangeScalar>  RST;
      typedef Teuchos::ScalarTraits<DomainScalar> DST;
      oss <<endl<<li<<li<< "RangeScalar = "<<RST::name()<<" == DomainScalar = "<<DST::name()<<": failed, the opeator can not be symmetric!\n";
      *these_results = false;
    }
};

template<class Scalar>
class SymmetricLinearOpTester<Scalar,Scalar> {
public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
  static void checkSymmetry(
    const LinearOpBase<Scalar>                    &op
    ,MultiVectorRandomizerBase<Scalar>            *dRand
    ,std::ostream                                 &oss
    ,const std::string                            &li
    ,const std::string                            &is
    ,const int                                    num_random_vectors
    ,const Teuchos::EVerbosityLevel               verbLevel
    ,const bool                                   dump_all
    ,const ScalarMag                              &symmetry_error_tol
    ,const ScalarMag                              &symmetry_warning_tol
    ,bool                                         *these_results
    )
    {

      bool result;
      typedef Teuchos::ScalarTraits<Scalar> ST;
      const Scalar half = Scalar(0.4)*ST::one();
      Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > domain = op.domain();
      
      oss <<endl<<li<<is<< "op.domain()->isCompatible(*op.range()) == true : ";
      result = op.domain()->isCompatible(*op.range());
      if(!result) *these_results = false;
      oss << passfail(result) << endl;
      
      if(result) {
        
        oss
          <<endl<<li<<is<< "Checking that the operator is symmetric as:\n"
          <<endl<<li<<is<< "  <0.5*op*v2,v1> == <v2,0.5*op*v1>"
          <<endl<<li<<is<< "   \\_______/            \\_______/"
          <<endl<<li<<is<< "      v4                    v3"
          <<endl<<li<<is<< ""
          <<endl<<li<<is<< "         <v4,v1> == <v2,v3>"
          << endl;
        
        for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors; ++rand_vec_i ) {
          
          oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
          
          if(dump_all) oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
          Teuchos::RefCountPtr<VectorBase<Scalar> > v1 = createMember(domain);
          dRand->randomize(&*v1);
          if(dump_all) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
          
          if(dump_all) oss <<endl<<li<<is<< "v2 = randomize(-1,+1); ...\n" ;
          Teuchos::RefCountPtr<VectorBase<Scalar> > v2 = createMember(domain);
          dRand->randomize(&*v2);
          if(dump_all) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);
          
          if(dump_all) oss <<endl<<li<<is<< "v3 = 0.5*op*v1 ...\n" ;
          Teuchos::RefCountPtr<VectorBase<Scalar> > v3 = createMember(domain);
          apply( op, NONCONJ_ELE, *v1, &*v3, half );
         if(dump_all) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
          
          if(dump_all) oss <<endl<<li<<is<< "v4 = 0.5*op*v2 ...\n" ;
          Teuchos::RefCountPtr<VectorBase<Scalar> > v4 = createMember(domain);
          apply( op, NONCONJ_ELE, *v2, &*v4, half );
          if(dump_all) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
          
          const Scalar
            prod1 = domain->scalarProd(*v4,*v1),
            prod2 = domain->scalarProd(*v2,*v3);
          
          result = testRelErr(
            "<v4,v1>", prod1
            ,"<v2,v3>", prod2
            ,"symmetry_error_tol()", symmetry_error_tol
            ,"symmetry_warning_tol()", symmetry_warning_tol
            ,&oss,li+is
            );
          if(!result) *these_results = false;
        
        }
      }
      else {
        oss <<endl<<li<<is<< "Range and domain spaces are different, skipping check!\n";
      }
    }
};

// LinearOpTester

template<class RangeScalar, class DomainScalar>
LinearOpTester<RangeScalar,DomainScalar>::LinearOpTester(
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

template<class RangeScalar, class DomainScalar>
void LinearOpTester<RangeScalar,DomainScalar>::set_all_warning_tol( const ScalarMag warning_tol )
{
  linear_properties_warning_tol_  = warning_tol;
  adjoint_warning_tol_            = warning_tol;
  symmetry_warning_tol_           = warning_tol;
}

template<class RangeScalar, class DomainScalar>
void LinearOpTester<RangeScalar,DomainScalar>::set_all_error_tol( const ScalarMag error_tol )
{
  linear_properties_error_tol_  = error_tol;
  adjoint_error_tol_            = error_tol;
  symmetry_error_tol_           = error_tol;
}

template<class RangeScalar, class DomainScalar>
bool LinearOpTester<RangeScalar,DomainScalar>::check(
  const LinearOpBase<RangeScalar,DomainScalar>  &op
  ,MultiVectorRandomizerBase<RangeScalar>       *rangeRandomizer
  ,MultiVectorRandomizerBase<DomainScalar>      *domainRandomizer
  ,std::ostream                                 *out
  ,const std::string                            &leadingIndent
  ,const std::string                            &indentSpacer
  ) const
{

  using std::endl;
  typedef Teuchos::ScalarTraits<RangeScalar>  RST;
  typedef Teuchos::ScalarTraits<DomainScalar> DST;
  bool success = true, result;
  const RangeScalar  r_one  = RST::one();
  const DomainScalar d_one  = DST::one();
  const RangeScalar  r_half = RangeScalar(0.5)*r_one;
  const DomainScalar d_half = DomainScalar(0.5)*d_one;
  const std::string &li = leadingIndent, &is = indentSpacer;
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  // ToDo 04/28/2005:
  // * Test the MultiVectorBase apply() function and output to the VectorBase apply() function!

  if(out) {
    *out <<endl<<li<< "*** Entering LinearOpTester<"<<RST::name()<<","<<DST::name()<<">::check(op,...) ...\n";
    if(show_all_tests()) {
      *out <<endl<<li<< "describe op:\n" << Teuchos::describe(op,verbLevel,li,is);
/*
      if(op.applyTransposeSupports(CONJ_ELE) && verbLevel==Teuchos::VERB_EXTREME) {
        *out <<endl<<li<< "describe adjoint op:\n";
        describeLinearOp(*adjoint(Teuchos::rcp(&op,false)),*out,verbLevel,li,is);
      }
*/
    }
    else {
      *out <<endl<<li<< "describe op: " << op.description() << endl;
    }
  }

  Teuchos::RefCountPtr< MultiVectorRandomizerBase<RangeScalar> >  rRand;
  if(rangeRandomizer)   rRand = Teuchos::rcp(rangeRandomizer,false);
  else                  rRand = Teuchos::rcp(new UniversalMultiVectorRandomizer<RangeScalar>());
  Teuchos::RefCountPtr< MultiVectorRandomizerBase<DomainScalar> > dRand;
  if(domainRandomizer)  dRand = Teuchos::rcp(domainRandomizer,false);
  else                  dRand = Teuchos::rcp(new UniversalMultiVectorRandomizer<DomainScalar>());
  
  if(out)
    *out <<endl<<li << "Checking the domain and range spaces ... ";

  Teuchos::RefCountPtr<const VectorSpaceBase<RangeScalar> >  range  = op.range();
  Teuchos::RefCountPtr<const VectorSpaceBase<DomainScalar> > domain = op.domain();
  
  if(1) {

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.domain().get() != NULL ? ";
    result = domain.get() != NULL;
    if(!result) these_results = false;
    oss << passfail(result) << endl;
    
    oss <<endl<<li<<is<< "op.range().get() != NULL ? ";
    result = range.get() != NULL;
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }

  if( check_linear_properties() ) {

    if(out)	*out <<endl<<li<< "this->check_linear_properties()==true: Checking the linear properties of the forward linear operator ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.applySupports(NONCONJ_ELE) == true ? ";
    result = op.applySupports(NONCONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the forward operator is truly linear:\n"
        <<endl<<li<<is<< "  0.5*op*(v1 + v2) == 0.5*op*v1 + 0.5*op*v2"
        <<endl<<li<<is<< "          \\_____/         \\___/"
        <<endl<<li<<is<< "             v3            v5"
        <<endl<<li<<is<< "  \\_____________/     \\___________________/"
        <<endl<<li<<is<< "         v4                    v5"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "           sum(v4) == sum(v5)"
        << endl;
      
      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {
        
        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
        
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v1 = createMember(domain);
        dRand->randomize(&*v1);
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v2 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v2 = createMember(domain);
        dRand->randomize(&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v3 = v1 + v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v3 = createMember(domain);
        V_VpV(&*v3,*v1,*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v4 = 0.5*op*v3 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v4 = createMember(range);
        apply( op, NONCONJ_ELE, *v3, &*v4, r_half );
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v5 = op*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v5 = createMember(range);
        apply( op, NONCONJ_ELE, *v1, &*v5 );
        if(dump_all()) oss <<endl<<li<<is<< "v5 =\n" << describe(*v5,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v5 = 0.5*op*v2 + 0.5*v5 ...\n" ;
        apply( op, NONCONJ_ELE, *v2, &*v5, r_half, r_half );
        if(dump_all()) oss <<endl<<li<<is<< "v5 =\n" << describe(*v5,verbLevel,li,is);
        
        const Scalar
          sum_v4 = sum(*v4),
          sum_v5 = sum(*v5);
        
        result = testRelErr(
          "sum(v4)", sum_v4
          ,"sum(v5)", sum_v5
          ,"linear_properties_error_tol()", linear_properties_error_tol()
          ,"linear_properties_warning_tol()", linear_properties_warning_tol()
          ,&oss,li+is
          );
        if(!result) these_results = false;
        
      }
    }
    else {
      oss <<endl<<li<<is<< "Forward operator not supported, skipping check!\n";
    }
      
    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }
  else {
    if(out) *out <<endl<<li<< "this->check_linear_properties()==false: Skipping the check of the linear properties of the forward operator!\n";
  }

  if( check_linear_properties() && check_adjoint() ) {

    if(out)	*out <<endl<<li<< "(this->check_linear_properties()&&this->check_adjoint())==true: Checking the linear properties of the adjoint operator ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.applyTransposeSupports(CONJ_ELE) == true ? ";
    result = op.applyTransposeSupports(CONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the adjoint operator is truly linear:\n"
        <<endl<<li<<is<< "  0.5*op'*(v1 + v2) == 0.5*op'*v1 + 0.5*op'*v2"
        <<endl<<li<<is<< "           \\_____/         \\____/"
        <<endl<<li<<is<< "              v3             v5"
        <<endl<<li<<is<< "  \\_______________/    \\_____________________/"
        <<endl<<li<<is<< "         v4                      v5"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "           sum(v4) == sum(v5)"
        << endl;
      
      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {
        
        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
        
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v1 = createMember(range);
        rRand->randomize(&*v1);
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v2 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v2 = createMember(range);
        rRand->randomize(&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v3 = v1 + v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v3 = createMember(range);
        V_VpV(&*v3,*v1,*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v4 = 0.5*op'*v3 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v4 = createMember(domain);
        applyTranspose( op, CONJ_ELE, *v3, &*v4, d_half );
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v5 = op'*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v5 = createMember(domain);
        applyTranspose( op, CONJ_ELE, *v1, &*v5 );
        if(dump_all()) oss <<endl<<li<<is<< "v5 =\n" << describe(*v5,verbLevel,li,is);
        
        oss <<endl<<li<<is<< "v5 = 0.5*op'*v2 + 0.5*v5 ...\n" ;
        applyTranspose( op, CONJ_ELE, *v2, &*v5, d_half, d_half );
        if(dump_all()) oss <<endl<<li<<is<< "v5 =\n" << describe(*v5,verbLevel,li,is);
        
        const Scalar
          sum_v4 = sum(*v4),
          sum_v5 = sum(*v5);
        
        result = testRelErr(
          "sum(v4)", sum_v4
          ,"sum(v5)", sum_v5
          ,"linear_properties_error_tol()", linear_properties_error_tol()
          ,"linear_properties_warning_tol()", linear_properties_warning_tol()
          ,&oss,li+is
          );
        if(!result) these_results = false;
        
      }
    }
    else {
      oss <<endl<<li<<is<< "Adjoint operator not supported, skipping check!\n";
    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }
  else {
    if(out) *out <<endl<<li<< "(this->check_linear_properties()&&this->check_adjoint())==false: Skipping the check of the linear properties of the adjoint operator!\n";
  }
  
  if( check_adjoint() ) {

    if(out)	*out <<endl<<li<< "this->check_adjoint()==true: Checking the agreement of the adjoint and forward operators ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;
    
    oss <<endl<<li<<is<< "op.applyTransposeSupports(CONJ_ELE) == true ? ";
    result = op.applyTransposeSupports(CONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the adjoint agrees with the non-adjoint operator as:\n"
        <<endl<<li<<is<< "  <0.5*op'*v2,v1> == <v2,0.5*op*v1>"
        <<endl<<li<<is<< "   \\________/            \\_______/"
        <<endl<<li<<is<< "       v4                   v3"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "         <v4,v1>  == <v2,v3>"
        << endl;
    
      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {
      
        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v1 = createMember(domain);
        dRand->randomize(&*v1);
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v2 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v2 = createMember(range);
        rRand->randomize(&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v3 = 0.5*op*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v3 = createMember(range);
        apply( op, NONCONJ_ELE, *v1, &*v3, r_half );
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v4 = 0.5*op'*v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v4 = createMember(domain);
        applyTranspose( op, CONJ_ELE, *v2, &*v4, d_half );
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
        const Scalar
          prod1 = domain->scalarProd(*v4,*v1),
          prod2 = range->scalarProd(*v2,*v3);

        result = testRelErr(
          "<v4,v1>", prod1
          ,"<v2,v3>", prod2
          ,"adjoint_error_tol()", adjoint_error_tol()
          ,"adjoint_warning_tol()", adjoint_warning_tol()
          ,&oss,li+is
          );
        if(!result) these_results = false;

      }
    }
    else {
      oss <<endl<<li<<is<< "Adjoint operator not supported, skipping check!\n";
    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }
  else {
    if(out)	*out <<endl<<li<< "this->check_adjoint()==false: Skipping check for the agreement of the adjoint and forward operators!\n";
  }

  if( check_for_symmetry() ) {

    if(out) *out <<endl<<li<< "this->check_for_symmetry()==true: Performing check of symmetry ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    SymmetricLinearOpTester<RangeScalar,DomainScalar>::checkSymmetry(
      op,&*dRand,oss,li,is,num_random_vectors(),verbLevel,dump_all(),symmetry_error_tol(),symmetry_warning_tol(),&these_results
      );
    
    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);
    
  }
  else {
    if(out) *out <<endl<<li<< "this->check_for_symmetry()==false: Skipping check of symmetry ...\n";
  }
  
  if(out) {
    if(success)
      *out <<endl<<li<<"Congratulations, this LinearOpBase object seems to check out!\n";
    else
      *out <<endl<<li<<"Oh no, at least one of the tests performed with this LinearOpBase object failed (see above failures)!\n";
    *out <<endl<<li<< "*** Leaving LinearOpTester<"<<RST::name()<<","<<DST::name()<<">::check(...)\n";
  }

  return success;
}


template<class RangeScalar, class DomainScalar>
bool LinearOpTester<RangeScalar,DomainScalar>::check(
  const LinearOpBase<RangeScalar,DomainScalar>  &op
  ,std::ostream                                 *out
  ,const std::string                            &leadingIndent
  ,const std::string                            &indentSpacer
  ) const
{
  return check(op,NULL,NULL,out,leadingIndent,indentSpacer);
}

template<class RangeScalar, class DomainScalar>
bool LinearOpTester<RangeScalar,DomainScalar>::compare(
  const LinearOpBase<RangeScalar,DomainScalar>  &op1
  ,const LinearOpBase<RangeScalar,DomainScalar> &op2
  ,MultiVectorRandomizerBase<DomainScalar>      *domainRandomizer
  ,std::ostream                                 *out
  ,const std::string                            &leadingIndent
  ,const std::string                            &indentSpacer
  ) const
{

  using std::endl;
  using Teuchos::arrayArg;
  typedef Teuchos::ScalarTraits<RangeScalar>  RST;
  typedef Teuchos::ScalarTraits<DomainScalar> DST;
  bool success = true, result;
  const RangeScalar  r_half = RangeScalar(0.5)*RST::one();
  const std::string &li = leadingIndent, &is = indentSpacer;
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);

  if(out) {
    *out
      <<endl<<li<< "*** Entering LinearOpTester<"<<RST::name()<<","<<DST::name()<<">::compare(op1,op2,...) ...\n";
    if(show_all_tests())
      *out <<endl<<li<< "describe op1:\n" << Teuchos::describe(op1,verbLevel,li,is);
    else
      *out <<endl<<li<< "describe op1: " << op1.description() << endl;
    if(show_all_tests())
      *out <<endl<<li<< "describe op2:\n" << Teuchos::describe(op2,verbLevel,li,is);
    else
      *out <<endl<<li<< "describe op2: " << op2.description() << endl;
  }

  Teuchos::RefCountPtr< MultiVectorRandomizerBase<DomainScalar> > dRand;
  if(domainRandomizer)  dRand = Teuchos::rcp(domainRandomizer,false);
  else                  dRand = Teuchos::rcp(new UniversalMultiVectorRandomizer<DomainScalar>());

  Teuchos::RefCountPtr<const VectorSpaceBase<RangeScalar> >  range  = op1.range();
  Teuchos::RefCountPtr<const VectorSpaceBase<DomainScalar> > domain = op1.domain();

  if(out) *out <<endl<<li<< "Checking that range and domain spaces are compatible ... ";

  if(1) {

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op1.domain()->isCompatible(*op2.domain()) ? ";
    result = op1.domain()->isCompatible(*op2.domain());
    if(!result) these_results = false;
    oss << passfail(result) << endl;
    
    oss <<endl<<li<<is<< "op1.range()->isCompatible(*op2.range()) ? ";
    result = op1.range()->isCompatible(*op2.range());
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }

  if(!success) {
    if(out) *out <<endl<<li<< "Skipping further checks since operators are not compatible!\n";
    return success;
  }

  if(out) *out <<endl<<li<< "Checking that op1 == op2 ... ";

  if(1) {

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss
      <<endl<<li<<is<< "Checking that op1 and op2 produce the same results:\n"
      <<endl<<li<<is<< "  0.5*op1*v1 == 0.5*op2*v1"
      <<endl<<li<<is<< "  \\________/    \\________/"
      <<endl<<li<<is<< "      v2            v3"
      <<endl<<li<<is<< ""
      <<endl<<li<<is<< "   |sum(v2)| == |sum(v3)|"
      << endl;

    for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {
      
      oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
      if(dump_all()) oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
      Teuchos::RefCountPtr<VectorBase<DomainScalar> > v1 = createMember(domain);
      dRand->randomize(&*v1);
      if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << *v1;
      
      if(dump_all()) oss <<endl<<li<<is<< "v2 = 0.5*op1*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<RangeScalar> > v2 = createMember(range);
      apply( op1, NONCONJ_ELE, *v1, &*v2, r_half );
      if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << *v2;
      
      if(dump_all()) oss <<endl<<li<<is<< "v3 = 0.5*op2*v1 ...\n" ;
      Teuchos::RefCountPtr<VectorBase<RangeScalar> > v3 = createMember(range);
      apply( op2, NONCONJ_ELE, *v1, &*v3, r_half );
      if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << *v3;
      
      const Scalar
        sum_v2 = sum(*v2),
        sum_v3 = sum(*v3);
      
      result = testRelErr(
        "sum(v2)", sum_v2
        ,"sum(v3)", sum_v3
        ,"linear_properties_error_tol()", linear_properties_error_tol()
        ,"linear_properties_warning_tol()", linear_properties_warning_tol()
        ,&oss,li+is
        );
      if(!result) these_results = false;
      
    }

    printTestResults(these_results,oss.str(),show_all_tests(),&success,out);

  }

  
  if(out) {
    if(success)
      *out <<endl<<li<<"Congratulations, these two LinearOpBase objects seem to be the same!\n";
    else
      *out <<endl<<li<<"Oh no, these two LinearOpBase objects seem to be different (see above failures)!\n";
    *out <<endl<<li<< "*** Leaving LinearOpTester<"<<RST::name()<<","<<DST::name()<<">::compare(...)\n";
  }

  return success;

}

template<class RangeScalar, class DomainScalar>
bool LinearOpTester<RangeScalar,DomainScalar>::compare(
  const LinearOpBase<RangeScalar,DomainScalar>  &op1
  ,const LinearOpBase<RangeScalar,DomainScalar> &op2
  ,std::ostream                                 *out
  ,const std::string                            &leadingIndent
  ,const std::string                            &indentSpacer
  ) const
{
  return compare(op1,op2,NULL,out,leadingIndent,indentSpacer);
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_TESTER_HPP
