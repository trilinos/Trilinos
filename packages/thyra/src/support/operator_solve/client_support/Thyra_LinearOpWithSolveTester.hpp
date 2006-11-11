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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_TESTER_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_TESTER_HPP

#include "Thyra_LinearOpWithSolveTesterDecl.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_describeLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

template<class RangeScalar, class DomainScalar>
LinearOpWithSolveTester<RangeScalar,DomainScalar>::LinearOpWithSolveTester(
  const bool                 check_forward_default
  ,const RangeScalarMag      forward_default_residual_warning_tol
  ,const RangeScalarMag      forward_default_residual_error_tol
  ,const DomainScalarMag     forward_default_solution_error_warning_tol
  ,const DomainScalarMag     forward_default_solution_error_error_tol
  ,const bool                check_forward_residual
  ,const RangeScalarMag      forward_residual_solve_tol
  ,const RangeScalarMag      forward_residual_slack_warning_tol
  ,const RangeScalarMag      forward_residual_slack_error_tol
  ,const bool                check_forward_solution_error
  ,const RangeScalarMag      forward_solution_error_solve_tol
  ,const RangeScalarMag      forward_solution_error_slack_warning_tol
  ,const RangeScalarMag      forward_solution_error_slack_error_tol
  ,const bool                check_adjoint_default
  ,const DomainScalarMag     adjoint_default_residual_warning_tol
  ,const DomainScalarMag     adjoint_default_residual_error_tol
  ,const RangeScalarMag      adjoint_default_solution_error_warning_tol
  ,const RangeScalarMag      adjoint_default_solution_error_error_tol
  ,const bool                check_adjoint_residual
  ,const DomainScalarMag     adjoint_residual_solve_tol
  ,const DomainScalarMag     adjoint_residual_slack_warning_tol
  ,const DomainScalarMag     adjoint_residual_slack_error_tol
  ,const bool                check_adjoint_solution_error
  ,const DomainScalarMag     adjoint_solution_error_solve_tol
  ,const DomainScalarMag     adjoint_solution_error_slack_warning_tol
  ,const DomainScalarMag     adjoint_solution_error_slack_error_tol
  ,const int                 num_random_vectors
  ,const bool                show_all_tests
  ,const bool                dump_all
  ,const int                 num_rhs
  )
  :check_forward_default_(check_forward_default)
  ,forward_default_residual_warning_tol_(forward_default_residual_warning_tol)
  ,forward_default_residual_error_tol_(forward_default_residual_error_tol)
  ,forward_default_solution_error_warning_tol_(forward_default_solution_error_warning_tol)
  ,forward_default_solution_error_error_tol_(forward_default_solution_error_error_tol)
  ,check_forward_residual_(check_forward_residual)
  ,forward_residual_solve_tol_(forward_residual_solve_tol)
  ,forward_residual_slack_warning_tol_(forward_residual_slack_warning_tol)
  ,forward_residual_slack_error_tol_(forward_residual_slack_error_tol)
  ,check_forward_solution_error_(check_forward_solution_error)
  ,forward_solution_error_solve_tol_(forward_solution_error_solve_tol)
  ,forward_solution_error_slack_warning_tol_(forward_solution_error_slack_warning_tol)
  ,forward_solution_error_slack_error_tol_(forward_solution_error_slack_error_tol)
  ,check_adjoint_default_(check_adjoint_default)
  ,adjoint_default_residual_warning_tol_(adjoint_default_residual_warning_tol)
  ,adjoint_default_residual_error_tol_(adjoint_default_residual_error_tol)
  ,adjoint_default_solution_error_warning_tol_(adjoint_default_solution_error_warning_tol)
  ,adjoint_default_solution_error_error_tol_(adjoint_default_solution_error_error_tol)
  ,check_adjoint_residual_(check_adjoint_residual)
  ,adjoint_residual_solve_tol_(adjoint_residual_solve_tol)
  ,adjoint_residual_slack_warning_tol_(adjoint_residual_slack_warning_tol)
  ,adjoint_residual_slack_error_tol_(adjoint_residual_slack_error_tol)
  ,check_adjoint_solution_error_(check_adjoint_solution_error)
  ,adjoint_solution_error_solve_tol_(adjoint_solution_error_solve_tol)
  ,adjoint_solution_error_slack_warning_tol_(adjoint_solution_error_slack_warning_tol)
  ,adjoint_solution_error_slack_error_tol_(adjoint_solution_error_slack_error_tol)
  ,num_random_vectors_(num_random_vectors)
  ,show_all_tests_(show_all_tests)
  ,dump_all_(dump_all)
  ,num_rhs_(num_rhs)
{}

template<class RangeScalar, class DomainScalar>
void LinearOpWithSolveTester<RangeScalar,DomainScalar>::turn_off_all_tests()
{
  check_forward_default_         = false;
  check_forward_residual_        = false;
  check_forward_solution_error_  = false;
  check_adjoint_default_         = false;
  check_adjoint_residual_        = false;
  check_adjoint_solution_error_  = false;
}

template<class RangeScalar, class DomainScalar>
void LinearOpWithSolveTester<RangeScalar,DomainScalar>::set_all_solve_tol( const ScalarMag solve_tol )
{
  forward_residual_solve_tol_ = solve_tol;
  forward_residual_solve_tol_ = solve_tol;
  forward_solution_error_solve_tol_ = solve_tol;
  adjoint_residual_solve_tol_ = solve_tol;
  adjoint_solution_error_solve_tol_ = solve_tol;
}

template<class RangeScalar, class DomainScalar>
void LinearOpWithSolveTester<RangeScalar,DomainScalar>::set_all_slack_warning_tol( const ScalarMag slack_warning_tol )
{
  forward_default_residual_warning_tol_ = slack_warning_tol;
  forward_default_solution_error_warning_tol_ = slack_warning_tol;
  forward_residual_slack_warning_tol_ = slack_warning_tol;
  forward_solution_error_slack_warning_tol_ = slack_warning_tol;
  adjoint_default_residual_warning_tol_ = slack_warning_tol;
  adjoint_default_solution_error_warning_tol_ = slack_warning_tol;
  adjoint_residual_slack_warning_tol_ = slack_warning_tol;
  adjoint_solution_error_slack_warning_tol_ = slack_warning_tol;
}

template<class RangeScalar, class DomainScalar>
void LinearOpWithSolveTester<RangeScalar,DomainScalar>::set_all_slack_error_tol( const ScalarMag slack_error_tol )
{
  forward_default_residual_error_tol_ = slack_error_tol;
  forward_default_solution_error_error_tol_ = slack_error_tol;
  forward_residual_slack_error_tol_ = slack_error_tol;
  forward_solution_error_slack_error_tol_ = slack_error_tol;
  adjoint_default_residual_error_tol_ = slack_error_tol;
  adjoint_default_solution_error_error_tol_ = slack_error_tol;
  adjoint_residual_slack_error_tol_ = slack_error_tol;
  adjoint_solution_error_slack_error_tol_ = slack_error_tol;
}
  
template<class RangeScalar, class DomainScalar>
bool LinearOpWithSolveTester<RangeScalar,DomainScalar>::check(
  const LinearOpWithSolveBase<RangeScalar,DomainScalar>   &op
  ,Teuchos::FancyOStream                                  *out_arg
  ) const
{

  using std::endl;
  using Teuchos::arrayArg;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  typedef Teuchos::VerboseObjectTempState<LinearOpWithSolveBase<RangeScalar,DomainScalar> > VOTS;
  typedef Teuchos::ScalarTraits<Scalar>        ST;
  typedef Teuchos::ScalarTraits<RangeScalar>   DST;
  typedef Teuchos::ScalarTraits<DomainScalar>  RST;
  bool success = true, result;
  const int num_rhs = this->num_rhs();
  Teuchos::RefCountPtr<FancyOStream> out = Teuchos::rcp(out_arg,false);
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);
  Teuchos::Time timer("");

  OSTab tab(out,1,"THYRA");

  if(out.get()) {
    *out <<endl<< "*** Entering LinearOpWithSolveTester<"<<ST::name()<<">::check(op,...) ...\n";
    if(show_all_tests()) {
      *out <<endl<< "describe forward op:\n" << Teuchos::describe(op,verbLevel);
      if(op.applyTransposeSupports(CONJ_ELE) && verbLevel==Teuchos::VERB_EXTREME) {
        *out <<endl<< "describe adjoint op:\n";
        describeLinearOp<RangeScalar,DomainScalar>(
          *adjoint(Teuchos::rcp_implicit_cast< const LinearOpBase<RangeScalar,DomainScalar> >(Teuchos::rcp(&op,false)))
          ,*out,verbLevel
          );
      }
    }
    else {
      *out <<endl<< "describe op: " << op.description() << endl;
    }
  }
  
  Teuchos::RefCountPtr<const VectorSpaceBase<RangeScalar> >   range  = op.range();
  Teuchos::RefCountPtr<const VectorSpaceBase<DomainScalar> >  domain = op.domain();
  
  if( check_forward_default() ) {

    if(out.get())	*out <<endl<< "this->check_forward_default()==true: Checking the default forward solve ... ";

    std::ostringstream ossStore;
    Teuchos::RefCountPtr<FancyOStream> oss = Teuchos::rcp(new FancyOStream(Teuchos::rcp(&ossStore,false)));
    if(out.get()) ossStore.copyfmt(*out);
    bool these_results = true;

    *oss <<endl<< "op.solveSupportsConj(NONCONJ_ELE) == true ? ";
    result = op.solveSupportsConj(NONCONJ_ELE);
    if(!result) these_results = false;
    *oss << passfail(result) << endl;

    if(result) {
    
      *oss
        <<endl<< "Checking that the forward default solve matches the forward operator:\n"
        <<endl<< "  inv(Op)*Op*v1 == v1"
        <<endl<< "          \\___/"
        <<endl<< "           v2"
        <<endl<< "  \\___________/"
        <<endl<< "         v3"
        <<endl<< ""
        <<endl<< "  v4 = v3-v1"
        <<endl<< "  v5 = Op*v3-v2"
        <<endl<< ""
        <<endl<< "  norm(v4)/norm(v1) <= forward_default_solution_error_error_tol()"
        <<endl<< "  norm(v5)/norm(v2) <= forward_default_residual_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        OSTab tab(oss);

        *oss <<endl<< "Random vector tests = " << rand_vec_i << endl;

        tab.incrTab();
      
        *oss <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v1 = createMembers(domain,num_rhs);
        Thyra::randomize( DomainScalar(-1.0), DomainScalar(+1.0), &*v1 );
        if(dump_all()) *oss <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *oss <<endl<< "v2 = Op*v1 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v2 = createMembers(range,num_rhs);
        timer.start(true);
        op.apply(NONCONJ_ELE,*v1,&*v2);
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *oss <<endl<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v3 = createMembers(domain,num_rhs);
        assign(&*v3,DST::zero());
        SolveStatus<Scalar> solveStatus;
        if(1){
          VOTS lowsTempState(Teuchos::rcp(&op,false),oss,verbLevel);
          timer.start(true);
          solveStatus = solve(op,NONCONJ_ELE,*v2,&*v3,static_cast<const SolveCriteria<Scalar>*>(0));
          timer.stop();
          OSTab(oss).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *oss <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *oss
          <<endl<< "solve status:\n";
        OSTab(oss).o() << solveStatus;

        *oss <<endl<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v4 = createMembers(domain,num_rhs);
        V_VmV( &*v4, *v3, *v1 );
        if(dump_all()) *oss <<endl<< "v4 =\n" << describe(*v4,verbLevel);
      
        *oss <<endl<< "v5 = Op*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v5 = createMembers(range,num_rhs);
        assign( &*v5, *v2 );
        timer.start(true);
        op.apply(NONCONJ_ELE,*v3,&*v5,Scalar(1.0),Scalar(-1.0));
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v5 =\n" << describe(*v5,verbLevel);

        std::vector<DomainScalarMag> norms_v1(num_rhs), norms_v4(num_rhs), norms_v4_norms_v1(num_rhs);
        norms(*v1,&norms_v1[0]);
        norms(*v4,&norms_v4[0]);
        std::transform(
          norms_v4.begin(),norms_v4.end(),norms_v1.begin(),norms_v4_norms_v1.begin()
          ,std::divides<DomainScalarMag>()
          );

        result = testMaxErrors(
          num_rhs, "norm(v4)/norm(v1)", &norms_v4_norms_v1[0]
          ,"forward_default_solution_error_error_tol()", forward_default_solution_error_error_tol()
          ,"forward_default_solution_error_warning_tol()", forward_default_solution_error_warning_tol()
          ,&*oss
          );
        if(!result) these_results = false;

        std::vector<RangeScalarMag> norms_v2(num_rhs), norms_v5(num_rhs), norms_v5_norms_v2(num_rhs);
        norms(*v2,&norms_v2[0]);
        norms(*v5,&norms_v5[0]);
        std::transform(
          norms_v5.begin(),norms_v5.end(),norms_v2.begin(),norms_v5_norms_v2.begin()
          ,std::divides<RangeScalarMag>()
          );

        result = testMaxErrors(
          num_rhs, "norm(v5)/norm(v2)", &norms_v5_norms_v2[0]
          ,"forward_default_residual_error_tol()", forward_default_residual_error_tol()
          ,"forward_default_residual_warning_tol()", forward_default_residual_warning_tol()
          ,&*oss
          );
        if(!result) these_results = false;
        
      }
    }
    else {
      *oss <<endl<< "Forward operator not supported, skipping check!\n";
    }

    printTestResults(these_results,ossStore.str(),show_all_tests(),&success,OSTab(out).get());

  }
  else {
    if(out.get()) *out <<endl<< "this->check_forward_default()==false: Skipping the check of the default forward solve!\n";
  }
  
  if( check_forward_residual() ) {

    if(out.get())	*out <<endl<< "this->check_forward_residual()==true: Checking the forward solve with a tolerance on the residual ... ";

    std::ostringstream ossStore;
    Teuchos::RefCountPtr<FancyOStream> oss = Teuchos::rcp(new FancyOStream(Teuchos::rcp(&ossStore,false)));
    if(out.get()) ossStore.copyfmt(*out);
    bool these_results = true;

    *oss <<endl<< "op.solveSupportsConj(NONCONJ_ELE) == true ? ";
    result = op.solveSupportsConj(NONCONJ_ELE);
    if(!result) these_results = false;
    *oss << passfail(result) << endl;

    if(result) {
    
      *oss
        <<endl<< "Checking that the forward solve matches the forward operator to a residual tolerance:\n"
        <<endl<< "  v3 = inv(Op)*Op*v1"
        <<endl<< "               \\___/"
        <<endl<< "                 v2"
        <<endl<< ""
        <<endl<< "  v4 = Op*v3-v2"
        <<endl<< ""
        <<endl<< "  norm(v4)/norm(v2) <= forward_residual_solve_tol() + forward_residual_slack_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        OSTab tab(oss);

        *oss <<endl<< "Random vector tests = " << rand_vec_i << endl;

        tab.incrTab();
      
        *oss <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v1 = createMembers(domain,num_rhs);
        Thyra::randomize( DomainScalar(-1.0), DomainScalar(+1.0), &*v1 );
        if(dump_all()) *oss <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *oss <<endl<< "v2 = Op*v1 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v2 = createMembers(range,num_rhs);
        timer.start(true);
        op.apply(NONCONJ_ELE,*v1,&*v2);
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *oss <<endl<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v3 = createMembers(domain,num_rhs);
        SolveCriteria<Scalar> solveCriteria(
          SolveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
          ,forward_residual_solve_tol()
          );
        assign(&*v3,DST::zero());
        SolveStatus<Scalar> solveStatus;
        if(1){
          VOTS lowsTempState(Teuchos::rcp(&op,false),oss,verbLevel);
          timer.start(true);
          solveStatus = solve<RangeScalar,DomainScalar>(op,NONCONJ_ELE,*v2,&*v3,&solveCriteria);
          timer.stop();
          OSTab(oss).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *oss <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *oss
          <<endl<< "solve status:\n";
        OSTab(oss).o() << solveStatus;
        result = solveStatus.solveStatus==SOLVE_STATUS_CONVERGED;
        if(!result) these_results = false;
        *oss
          <<endl<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(result)<<endl;
      
        *oss <<endl<< "v4 = Op*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v4 = createMembers(range,num_rhs);
        assign( &*v4, *v2 );
        timer.start(true);
        op.apply(NONCONJ_ELE,*v3,&*v4,Scalar(1.0),Scalar(-1.0));
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v4 =\n" << describe(*v4,verbLevel);

        std::vector<DomainScalarMag> norms_v2(num_rhs), norms_v4(num_rhs), norms_v4_norms_v2(num_rhs);
        norms(*v2,&norms_v2[0]);
        norms(*v4,&norms_v4[0]);
        std::transform(
          norms_v4.begin(),norms_v4.end(),norms_v2.begin(),norms_v4_norms_v2.begin()
          ,std::divides<DomainScalarMag>()
          );

        result = testMaxErrors(
          num_rhs, "norm(v4)/norm(v2)", &norms_v4_norms_v2[0]
          ,"forward_residual_solve_tol()+forward_residual_slack_error_tol()", DomainScalarMag(forward_residual_solve_tol()+forward_residual_slack_error_tol())
          ,"forward_residual_solve_tol()_slack_warning_tol()", DomainScalarMag(forward_residual_solve_tol()+forward_residual_slack_warning_tol())
          ,&*oss
          );
        if(!result) these_results = false;

      }
    }
    else {
      *oss <<endl<< "Forward operator not supported, skipping check!\n";
    }

    printTestResults(these_results,ossStore.str(),show_all_tests(),&success,OSTab(out).get());

  }
  else {
    if(out.get()) *out <<endl<< "this->check_forward_residual()==false: Skipping the check of the forward solve with a tolerance on the residual!\n";
  }
  
  if( check_adjoint_default() ) {

    if(out.get())	*out <<endl<< "this->check_adjoint_default()==true: Checking the default adjoint solve ... ";

    std::ostringstream ossStore;
    Teuchos::RefCountPtr<FancyOStream> oss = Teuchos::rcp(new FancyOStream(Teuchos::rcp(&ossStore,false)));
    if(out.get()) ossStore.copyfmt(*out);
    bool these_results = true;

    *oss <<endl<< "op.solveTransposeSupportsConj(CONJ_ELE) == true ? ";
    result = op.solveTransposeSupportsConj(CONJ_ELE);
    if(!result) these_results = false;
    *oss << passfail(result) << endl;

    if(result) {
    
      *oss
        <<endl<< "Checking that the adjoint default solve matches the adjoint operator:\n"
        <<endl<< "  inv(Op')*Op'*v1 == v1"
        <<endl<< "           \\____/"
        <<endl<< "             v2"
        <<endl<< "  \\_____________/"
        <<endl<< "         v3"
        <<endl<< ""
        <<endl<< "  v4 = v3-v1"
        <<endl<< "  v5 = Op'*v3-v2"
        <<endl<< ""
        <<endl<< "  norm(v4)/norm(v1) <= adjoint_default_solution_error_error_tol()"
        <<endl<< "  norm(v5)/norm(v2) <= adjoint_default_residual_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        OSTab tab(oss);

        *oss <<endl<< "Random vector tests = " << rand_vec_i << endl;
        
        tab.incrTab();
      
        *oss <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v1 = createMembers(range,num_rhs);
        Thyra::randomize( RangeScalar(-1.0), RangeScalar(+1.0), &*v1 );
        if(dump_all()) *oss <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *oss <<endl<< "v2 = Op'*v1 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v2 = createMembers(domain,num_rhs);
        timer.start(true);
        op.applyTranspose(CONJ_ELE,*v1,&*v2);
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *oss <<endl<< "v3 = inv(Op')*v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v3 = createMembers(range,num_rhs);
        assign(&*v3,DST::zero());
        SolveStatus<Scalar> solveStatus;
        if(1){
          VOTS lowsTempState(Teuchos::rcp(&op,false),oss,verbLevel);
          timer.start(true);
          solveStatus = solveTranspose(op,CONJ_ELE,*v2,&*v3,static_cast<const SolveCriteria<Scalar>*>(0));
          timer.stop();
          OSTab(oss).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *oss <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *oss
          <<endl<< "solve status:\n";
        OSTab(oss).o() << solveStatus;

        *oss <<endl<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v4 = createMembers(range,num_rhs);
        V_VmV( &*v4, *v3, *v1 );
        if(dump_all()) *oss <<endl<< "v4 =\n" << describe(*v4,verbLevel);
      
        *oss <<endl<< "v5 = Op'*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v5 = createMembers(domain,num_rhs);
        assign( &*v5, *v2 );
        timer.start(true);
        op.applyTranspose(CONJ_ELE,*v3,&*v5,Scalar(1.0),Scalar(-1.0));
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v5 =\n" << describe(*v5,verbLevel);

        std::vector<RangeScalarMag> norms_v1(num_rhs), norms_v4(num_rhs), norms_v4_norms_v1(num_rhs);
        norms(*v1,&norms_v1[0]);
        norms(*v4,&norms_v4[0]);
        std::transform(
          norms_v4.begin(),norms_v4.end(),norms_v1.begin(),norms_v4_norms_v1.begin()
          ,std::divides<RangeScalarMag>()
          );

        result = testMaxErrors(
          num_rhs, "norm(v4)/norm(v1)", &norms_v4_norms_v1[0]
          ,"adjoint_default_solution_error_error_tol()", adjoint_default_solution_error_error_tol()
          ,"adjoint_default_solution_error_warning_tol()", adjoint_default_solution_error_warning_tol()
          ,&*oss
          );
        if(!result) these_results = false;

        std::vector<DomainScalarMag> norms_v2(num_rhs), norms_v5(num_rhs), norms_v5_norms_v2(num_rhs);
        norms(*v2,&norms_v2[0]);
        norms(*v5,&norms_v5[0]);
        std::transform(
          norms_v5.begin(),norms_v5.end(),norms_v2.begin(),norms_v5_norms_v2.begin()
          ,std::divides<DomainScalarMag>()
          );

        result = testMaxErrors(
          num_rhs, "norm(v5)/norm(v2)", &norms_v5_norms_v2[0]
          ,"adjoint_default_residual_error_tol()", adjoint_default_residual_error_tol()
          ,"adjoint_default_residual_warning_tol()", adjoint_default_residual_warning_tol()
          ,&*oss
          );
        if(!result) these_results = false;

      }
    }
    else {
      *oss <<endl<< "Adjoint operator not supported, skipping check!\n";
    }

    printTestResults(these_results,ossStore.str(),show_all_tests(),&success,OSTab(out).get());

  }
  else {
    if(out.get()) *out <<endl<< "this->check_adjoint_default()==false: Skipping the check of the adjoint solve with a default tolerance!\n";
  }
  
  if( check_adjoint_residual() ) {

    if(out.get())	*out <<endl<< "this->check_adjoint_residual()==true: Checking the adjoint solve with a tolerance on the residual ... ";
    
    std::ostringstream ossStore;
    Teuchos::RefCountPtr<FancyOStream> oss = Teuchos::rcp(new FancyOStream(Teuchos::rcp(&ossStore,false)));
    if(out.get()) ossStore.copyfmt(*out);
    bool these_results = true;

    *oss <<endl<< "op.solveTransposeSupportsConj(CONJ_ELE) == true ? ";
    result = op.solveTransposeSupportsConj(CONJ_ELE);
    if(!result) these_results = false;
    *oss << passfail(result) << endl;

    if(result) {
    
      *oss
        <<endl<< "Checking that the adjoint solve matches the adjoint operator to a residual tolerance:\n"
        <<endl<< "  v3 = inv(Op')*Op'*v1"
        <<endl<< "                \\____/"
        <<endl<< "                  v2"
        <<endl<< ""
        <<endl<< "  v4 = Op'*v3-v2"
        <<endl<< ""
        <<endl<< "  norm(v4)/norm(v2) <= adjoint_residual_solve_tol() + adjoint_residual_slack_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        OSTab tab(oss);

        *oss <<endl<< "Random vector tests = " << rand_vec_i << endl;

        tab.incrTab();
      
        *oss <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v1 = createMembers(range,num_rhs);
        Thyra::randomize( RangeScalar(-1.0), RangeScalar(+1.0), &*v1 );
        if(dump_all()) *oss <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *oss <<endl<< "v2 = Op'*v1 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v2 = createMembers(domain,num_rhs);
        timer.start(true);
        op.applyTranspose(CONJ_ELE,*v1,&*v2);
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *oss <<endl<< "v3 = inv(Op')*v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<RangeScalar> > v3 = createMembers(range,num_rhs);
        SolveCriteria<Scalar> solveCriteria(
          SolveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
          ,adjoint_residual_solve_tol()
          );
        assign(&*v3,RST::zero());
        SolveStatus<Scalar> solveStatus;
        if(1){
          VOTS lowsTempState(Teuchos::rcp(&op,false),oss,verbLevel);
          timer.start(true);
          solveStatus = solveTranspose(op,CONJ_ELE,*v2,&*v3,&solveCriteria);
          timer.stop();
          OSTab(oss).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *oss <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *oss
          <<endl<< "solve status:\n";
        OSTab(oss).o() << solveStatus;
        result = solveStatus.solveStatus==SOLVE_STATUS_CONVERGED;
        if(!result) these_results = false;
        *oss
          <<endl<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(result)<<endl;
        if(solveStatus.achievedTol==SolveStatus<RangeScalar>::unknownTolerance())
          *oss <<endl<<"achievedTol==unknownTolerance(): Setting achievedTol = adjoint_residual_solve_tol() = "<<adjoint_residual_solve_tol()<<endl;
      
        *oss <<endl<< "v4 = Op'*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<MultiVectorBase<DomainScalar> > v4 = createMembers(domain,num_rhs);
        assign( &*v4, *v2 );
        timer.start(true);
        op.applyTranspose(CONJ_ELE,*v3,&*v4,Scalar(1.0),Scalar(-1.0));
        timer.stop();
        OSTab(oss).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *oss <<endl<< "v4 =\n" << describe(*v4,verbLevel);

        std::vector<RangeScalarMag> norms_v2(num_rhs), norms_v4(num_rhs), norms_v4_norms_v2(num_rhs);
        norms(*v2,&norms_v2[0]);
        norms(*v4,&norms_v4[0]);
        std::transform(
          norms_v4.begin(),norms_v4.end(),norms_v2.begin(),norms_v4_norms_v2.begin()
          ,std::divides<RangeScalarMag>()
          );

        result = testMaxErrors(
          num_rhs, "norm(v4)/norm(v2)", &norms_v4_norms_v2[0]
          ,"adjoint_residual_solve_tol()+adjoint_residual_slack_error_tol()", RangeScalarMag(adjoint_residual_solve_tol()+adjoint_residual_slack_error_tol())
          ,"adjoint_residual_solve_tol()_slack_warning_tol()", RangeScalarMag(adjoint_residual_solve_tol()+adjoint_residual_slack_warning_tol())
          ,&*oss
          );
        if(!result) these_results = false;

      }
    }
    else {
      *oss <<endl<< "Adjoint operator not supported, skipping check!\n";
    }

    printTestResults(these_results,ossStore.str(),show_all_tests(),&success,OSTab(out).get());

  }
  else {
    if(out.get()) *out <<endl<< "this->check_adjoint_residual()==false: Skipping the check of the adjoint solve with a tolerance on the residual!\n";
  }
  
  if(out.get()) {
    if(success)
      *out <<endl<<"Congratulations, this LinearOpWithSolveBase object seems to check out!\n";
    else
      *out <<endl<<"Oh no, at least one of the tests performed with this LinearOpWithSolveBase object failed (see above failures)!\n";
    *out <<endl<< "*** Leaving LinearOpWithSolveTester<"<<ST::name()<<">::check(...)\n";
  }
  
  return success;
  
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_TESTER_HPP
