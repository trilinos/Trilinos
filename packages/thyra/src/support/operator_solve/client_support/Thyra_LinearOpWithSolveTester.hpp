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
  ,std::ostream                                           *out
  ,const std::string                                      &leadingIndent
  ,const std::string                                      &indentSpacer
  ) const
{

  using std::endl;
  using Teuchos::arrayArg;
  typedef Teuchos::ScalarTraits<Scalar>        ST;
  typedef Teuchos::ScalarTraits<RangeScalar>   DST;
  typedef Teuchos::ScalarTraits<DomainScalar>  RST;
  bool success = true, result;
  const std::string &li = leadingIndent, &is = indentSpacer;
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);
  
  if(out) {
    *out <<endl<<li<< "*** Entering LinearOpWithSolveTester<"<<ST::name()<<">::check(op,...) ...\n";
    if(show_all_tests()) {
      *out <<endl<<li<< "describe forward op:\n" << Teuchos::describe(op,verbLevel,li,is);
      if(op.applyTransposeSupports(CONJ_ELE) && verbLevel==Teuchos::VERB_EXTREME) {
        *out <<endl<<li<< "describe adjoint op:\n";
        describeLinearOp(
          *adjoint(Teuchos::rcp_implicit_cast< const LinearOpBase<RangeScalar,DomainScalar> >(Teuchos::rcp(&op,false)))
          ,*out,verbLevel,li,is
          );
      }
    }
    else {
      *out <<endl<<li<< "describe op: " << op.description() << endl;
    }
  }
  
  Teuchos::RefCountPtr<const VectorSpaceBase<RangeScalar> >   range  = op.range();
  Teuchos::RefCountPtr<const VectorSpaceBase<DomainScalar> >  domain = op.domain();
  
  if( check_forward_default() ) {

    if(out)	*out <<endl<<li<< "this->check_forward_default()==true: Checking the default forward solve ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.solveSupportsConj(NONCONJ_ELE) == true ? ";
    result = op.solveSupportsConj(NONCONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the forward default solve matches the forward operator:\n"
        <<endl<<li<<is<< "  inv(Op)*Op*v1 == v1"
        <<endl<<li<<is<< "          \\___/"
        <<endl<<li<<is<< "           v2"
        <<endl<<li<<is<< "  \\___________/"
        <<endl<<li<<is<< "         v3"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  v4 = v3-v1"
        <<endl<<li<<is<< "  v5 = Op*v3-v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  norm(v4)/norm(v1) <= forward_default_solution_error_error_tol()"
        <<endl<<li<<is<< "  norm(v5)/norm(v2) <= forward_default_residual_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v1 = createMember(domain);
        Thyra::randomize( DomainScalar(-1.0), DomainScalar(+1.0), &*v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v2 = Op*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v2 = createMember(range);
        op.apply(NONCONJ_ELE,*v1,&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);

        oss <<endl<<li<<is<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v3 = createMember(domain);
        assign(&*v3,DST::zero());
        SolveStatus<Scalar> solveStatus = solve(op,NONCONJ_ELE,*v2,&*v3,static_cast<const SolveCriteria<Scalar>*>(0));
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        oss
          <<endl<<li<<is<< "solve status:"
          <<endl<<li<<is<<is<< "solveStatus = " << toString(solveStatus.solveStatus)
          <<endl<<li<<is<<is<< "achievedTol = " << SolveStatus<RangeScalar>::achievedTolToString(solveStatus.achievedTol)
          <<endl<<li<<is<<is<< "iterations  = " << solveStatus.iterations
          <<endl<<li<<is<<is<< "message     = \"" << solveStatus.message << "\""
          <<endl;

        oss <<endl<<li<<is<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v4 = createMember(domain);
        V_VmV( &*v4, *v3, *v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v5 = Op*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v5 = createMember(range);
        assign( &*v5, *v2 );
        op.apply(NONCONJ_ELE,*v3,&*v5,Scalar(1.0),Scalar(-1.0));
        if(dump_all()) oss <<endl<<li<<is<< "v5 =\n" << describe(*v5,verbLevel,li,is);
      
        const DomainScalarMag
          norm_v1 = norm(*v1),
          norm_v4 = norm(*v4),
          norm_v4_norm_v1 = norm_v4/norm_v1;

        result = testMaxErr(
          "norm(v4)/norm(v1)", norm_v4_norm_v1
          ,"forward_default_solution_error_error_tol()", forward_default_solution_error_error_tol()
          ,"forward_default_solution_error_warning_tol()", forward_default_solution_error_warning_tol()
          ,&oss,li+is
          );
        if(!result) these_results = false;
      
        const RangeScalarMag
          norm_v2 = norm(*v2),
          norm_v5 = norm(*v5),
          norm_v5_norm_v2 = norm_v5/norm_v2;

        result = testMaxErr(
          "norm(v5)/norm(v2)", norm_v5_norm_v2
          ,"forward_default_residual_error_tol()", forward_default_residual_error_tol()
          ,"forward_default_residual_warning_tol()", forward_default_residual_warning_tol()
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
    if(out) *out <<endl<<li<< "this->check_forward_default()==false: Skipping the check of the default forward solve!\n";
  }
  
  if( check_forward_residual() ) {

    if(out)	*out <<endl<<li<< "this->check_forward_residual()==true: Checking the forward solve with a tolerance on the residual ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.solveSupportsConj(NONCONJ_ELE) == true ? ";
    result = op.solveSupportsConj(NONCONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the forward solve matches the forward operator to a residual tolerance:\n"
        <<endl<<li<<is<< "  v3 = inv(Op)*Op*v1"
        <<endl<<li<<is<< "               \\___/"
        <<endl<<li<<is<< "                 v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  v4 = Op*v3-v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  norm(v4)/norm(v2) <= forward_residual_solve_tol() + forward_residual_slack_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v1 = createMember(domain);
        Thyra::randomize( DomainScalar(-1.0), DomainScalar(+1.0), &*v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v2 = Op*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v2 = createMember(range);
        op.apply(NONCONJ_ELE,*v1,&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);

        oss <<endl<<li<<is<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v3 = createMember(domain);
        SolveCriteria<Scalar> solveCriteria(SOLVE_TOL_REL_RESIDUAL_NORM,forward_residual_solve_tol());
        assign(&*v3,DST::zero());
        SolveStatus<Scalar> solveStatus = solve<RangeScalar,DomainScalar>(op,NONCONJ_ELE,*v2,&*v3,&solveCriteria);
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        oss
          <<endl<<li<<is<< "solve status:"
          <<endl<<li<<is<<is<< "solveStatus = " << toString(solveStatus.solveStatus)
          <<endl<<li<<is<<is<< "achievedTol = " << SolveStatus<RangeScalar>::achievedTolToString(solveStatus.achievedTol)
          <<endl<<li<<is<<is<< "iterations  = " << solveStatus.iterations
          <<endl<<li<<is<<is<< "message     = \"" << solveStatus.message << "\""
          <<endl;
        oss
          <<endl<<li<<is<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(solveStatus.solveStatus==SOLVE_STATUS_CONVERGED)<<endl;
      
        oss <<endl<<li<<is<< "v4 = Op*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v4 = createMember(range);
        assign( &*v4, *v2 );
        op.apply(NONCONJ_ELE,*v3,&*v4,Scalar(1.0),Scalar(-1.0));
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
        const RangeScalarMag
          norm_v2 = norm(*v2),
          norm_v4 = norm(*v4),
          norm_v4_norm_v2 = norm_v4/norm_v2;

        result = testMaxErr(
          "norm(v4)/norm(v2)", norm_v4_norm_v2
          ,"forward_residual_solve_tol()+forward_residual_slack_error_tol()", RangeScalarMag(forward_residual_solve_tol()+forward_residual_slack_error_tol())
          ,"forward_residual_solve_tol()_slack_warning_tol()", RangeScalarMag(forward_residual_solve_tol()+forward_residual_slack_warning_tol())
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
    if(out) *out <<endl<<li<< "this->check_forward_residual()==false: Skipping the check of the forward solve with a tolerance on the residual!\n";
  }
  
  if( check_forward_solution_error() ) {

    if(out)	*out <<endl<<li<< "this->check_forward_solution_error()==true: Checking the forward solve with a tolerance on the solution error ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.solveSupportsConj(NONCONJ_ELE) == true ? ";
    result = op.solveSupportsConj(NONCONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the forward solve matches the forward operator to a solution error tolerance:\n"
        <<endl<<li<<is<< "  v3 = inv(Op)*Op*v1"
        <<endl<<li<<is<< "               \\___/"
        <<endl<<li<<is<< "                 v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  v4 = v3-v1"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  norm(v4)/norm(v1) <= forward_solution_error_solve_tol() + forward_solution_error_slack_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v1 = createMember(domain);
        Thyra::randomize( DomainScalar(-1.0), DomainScalar(+1.0), &*v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v2 = Op*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v2 = createMember(range);
        op.apply(NONCONJ_ELE,*v1,&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);

        oss <<endl<<li<<is<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v3 = createMember(domain);
        SolveCriteria<Scalar> solveCriteria(SOLVE_TOL_REL_SOLUTION_ERR_NORM,forward_solution_error_solve_tol());
        assign(&*v3,DST::zero());
        SolveStatus<Scalar> solveStatus = solve(op,NONCONJ_ELE,*v2,&*v3,&solveCriteria);
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        oss
          <<endl<<li<<is<< "solve status:"
          <<endl<<li<<is<<is<< "solveStatus = " << toString(solveStatus.solveStatus)
          <<endl<<li<<is<<is<< "achievedTol = " << SolveStatus<RangeScalar>::achievedTolToString(solveStatus.achievedTol)
          <<endl<<li<<is<<is<< "iterations  = " << solveStatus.iterations
          <<endl<<li<<is<<is<< "message     = \"" << solveStatus.message << "\""
          <<endl;
        oss
          <<endl<<li<<is<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(solveStatus.solveStatus==SOLVE_STATUS_CONVERGED)<<endl;
      
        oss <<endl<<li<<is<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v4 = createMember(domain);
        V_VmV( &*v4, *v3, *v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
        const RangeScalarMag
          norm_v1 = norm(*v1),
          norm_v4 = norm(*v4),
          norm_v4_norm_v1 = norm_v4/norm_v1;

        result = testMaxErr(
          "norm(v4)/norm(v1)", norm_v4_norm_v1
          ,"forward_solution_error_solve_tol()+forward_solution_error_slack_error_tol()", RangeScalarMag(forward_solution_error_solve_tol()+forward_solution_error_slack_error_tol())
          ,"forward_solution_error_solve_tol()+forward_solution_error_slack_warning_tol()", RangeScalarMag(forward_solution_error_solve_tol()+forward_solution_error_slack_warning_tol())
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
    if(out) *out <<endl<<li<< "this->check_forward_solution_error()==false: Skipping the check of the forward solve with a tolerance on the solution_error!\n";
  }
  
  if( check_adjoint_default() ) {

    if(out)	*out <<endl<<li<< "this->check_adjoint_default()==true: Checking the default adjoint solve ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.solveTransposeSupportsConj(CONJ_ELE) == true ? ";
    result = op.solveTransposeSupportsConj(CONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the adjoint default solve matches the adjoint operator:\n"
        <<endl<<li<<is<< "  inv(Op')*Op'*v1 == v1"
        <<endl<<li<<is<< "           \\____/"
        <<endl<<li<<is<< "             v2"
        <<endl<<li<<is<< "  \\_____________/"
        <<endl<<li<<is<< "         v3"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  v4 = v3-v1"
        <<endl<<li<<is<< "  v5 = Op'*v3-v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  norm(v4)/norm(v1) <= adjoint_default_solution_error_error_tol()"
        <<endl<<li<<is<< "  norm(v5)/norm(v2) <= adjoint_default_residual_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v1 = createMember(range);
        Thyra::randomize( RangeScalar(-1.0), RangeScalar(+1.0), &*v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v2 = Op*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v2 = createMember(domain);
        op.apply(NONCONJ_ELE,*v1,&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);

        oss <<endl<<li<<is<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v3 = createMember(range);
        assign(&*v3,DST::zero());
        SolveStatus<Scalar> solveStatus = solve(op,NONCONJ_ELE,*v2,&*v3,static_cast<const SolveCriteria<Scalar>*>(0));
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        oss
          <<endl<<li<<is<< "solve status:"
          <<endl<<li<<is<<is<< "solveStatus = " << toString(solveStatus.solveStatus)
          <<endl<<li<<is<<is<< "achievedTol = " << SolveStatus<DomainScalar>::achievedTolToString(solveStatus.achievedTol)
          <<endl<<li<<is<<is<< "iterations  = " << solveStatus.iterations
          <<endl<<li<<is<<is<< "message     = \"" << solveStatus.message << "\""
          <<endl;

        oss <<endl<<li<<is<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v4 = createMember(range);
        V_VmV( &*v4, *v3, *v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v5 = Op*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v5 = createMember(domain);
        assign( &*v5, *v2 );
        op.apply(NONCONJ_ELE,*v3,&*v5,Scalar(1.0),Scalar(-1.0));
        if(dump_all()) oss <<endl<<li<<is<< "v5 =\n" << describe(*v5,verbLevel,li,is);
      
        const RangeScalarMag
          norm_v1 = norm(*v1),
          norm_v4 = norm(*v4),
          norm_v4_norm_v1 = norm_v4/norm_v1;

        result = testMaxErr(
          "norm(v4)/norm(v1)", norm_v4_norm_v1
          ,"adjoint_default_solution_error_error_tol()", adjoint_default_solution_error_error_tol()
          ,"adjoint_default_solution_error_warning_tol()", adjoint_default_solution_error_warning_tol()
          ,&oss,li+is
          );
        if(!result) these_results = false;
      
        const DomainScalarMag
          norm_v2 = norm(*v2),
          norm_v5 = norm(*v5),
          norm_v5_norm_v2 = norm_v5/norm_v2;

        result = testMaxErr(
          "norm(v5)/norm(v2)", norm_v5_norm_v2
          ,"adjoint_default_residual_error_tol()", adjoint_default_residual_error_tol()
          ,"adjoint_default_residual_warning_tol()", adjoint_default_residual_warning_tol()
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
    if(out) *out <<endl<<li<< "this->check_adjoint_default()==false: Skipping the check of the adjoint solve with a default tolerance!\n";
  }
  
  if( check_adjoint_residual() ) {

    if(out)	*out <<endl<<li<< "this->check_adjoint_residual()==true: Checking the adjoint solve with a tolerance on the residual ... ";
    
    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.solveTransposeSupportsConj(CONJ_ELE) == true ? ";
    result = op.solveTransposeSupportsConj(CONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the adjoint solve matches the adjoint operator to a residual tolerance:\n"
        <<endl<<li<<is<< "  v3 = inv(Op')*Op'*v1"
        <<endl<<li<<is<< "                \\____/"
        <<endl<<li<<is<< "                  v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  v4 = Op'*v3-v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  norm(v4)/norm(v2) <= adjoint_residual_solve_tol() + adjoint_residual_slack_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v1 = createMember(range);
        Thyra::randomize( RangeScalar(-1.0), RangeScalar(+1.0), &*v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v2 = Op'*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v2 = createMember(domain);
        op.applyTranspose(CONJ_ELE,*v1,&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);

        oss <<endl<<li<<is<< "v3 = inv(Op')*v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v3 = createMember(range);
        SolveCriteria<Scalar> solveCriteria(SOLVE_TOL_REL_RESIDUAL_NORM,adjoint_residual_solve_tol());
        assign(&*v3,RST::zero());
        SolveStatus<Scalar> solveStatus = solveTranspose(op,CONJ_ELE,*v2,&*v3,&solveCriteria);
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        oss
          <<endl<<li<<is<< "solve status:"
          <<endl<<li<<is<<is<< "solveStatus = " << toString(solveStatus.solveStatus)
          <<endl<<li<<is<<is<< "achievedTol = " << SolveStatus<RangeScalar>::achievedTolToString(solveStatus.achievedTol)
          <<endl<<li<<is<<is<< "iterations  = " << solveStatus.iterations
          <<endl<<li<<is<<is<< "message     = \"" << solveStatus.message << "\""
          <<endl;
        oss
          <<endl<<li<<is<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(solveStatus.solveStatus==SOLVE_STATUS_CONVERGED)<<endl;
        if(solveStatus.achievedTol==SolveStatus<RangeScalar>::unknownTolerance())
          oss <<endl<<li<<is<<"achievedTol==unknownTolerance(): Setting achievedTol = adjoint_residual_solve_tol() = "<<adjoint_residual_solve_tol()<<endl;
      
        oss <<endl<<li<<is<< "v4 = Op'*v3 - v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v4 = createMember(domain);
        assign( &*v4, *v2 );
        op.applyTranspose(CONJ_ELE,*v3,&*v4,Scalar(1.0),Scalar(-1.0));
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
        const DomainScalarMag
          norm_v2 = norm(*v2),
          norm_v4 = norm(*v4),
          norm_v4_norm_v2 = norm_v4/norm_v2;

        result = testMaxErr(
          "norm(v4)/norm(v2)", norm_v4_norm_v2
          ,"adjoint_residual_solve_tol()+adjoint_residual_slack_error_tol()", DomainScalarMag(adjoint_residual_solve_tol()+adjoint_residual_slack_error_tol())
          ,"adjoint_residual_solve_tol()+adjoint_residual_slack_warning_tol()", DomainScalarMag(adjoint_residual_solve_tol()+adjoint_residual_slack_warning_tol())
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
    if(out) *out <<endl<<li<< "this->check_adjoint_residual()==false: Skipping the check of the adjoint solve with a tolerance on the residual!\n";
  }
  
  if( check_adjoint_solution_error() ) {

    if(out)	*out <<endl<<li<< "this->check_adjoint_solution_error()==true: Checking the adjoint solve with a tolerance on the solution error ... ";

    std::ostringstream oss;
    if(out) oss.copyfmt(*out);
    bool these_results = true;

    oss <<endl<<li<<is<< "op.solveTransposeSupportsConj(CONJ_ELE) == true ? ";
    result = op.solveTransposeSupportsConj(CONJ_ELE);
    if(!result) these_results = false;
    oss << passfail(result) << endl;

    if(result) {
    
      oss
        <<endl<<li<<is<< "Checking that the adjoint solve matches the adjoint operator to a solution error tolerance:\n"
        <<endl<<li<<is<< "  v3 = inv(Op')*Op'*v1"
        <<endl<<li<<is<< "                \\____/"
        <<endl<<li<<is<< "                  v2"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  v4 = v3-v1"
        <<endl<<li<<is<< ""
        <<endl<<li<<is<< "  norm(v4)/norm(v1) <= adjoint_solution_error_solve_tol() + adjoint_solution_error_slack_error_tol()"
        <<endl;

      for( int rand_vec_i = 1; rand_vec_i <= num_random_vectors(); ++rand_vec_i ) {

        oss <<endl<<li<<is<< "Random vector tests = " << rand_vec_i << endl;
      
        oss <<endl<<li<<is<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v1 = createMember(range);
        Thyra::randomize( RangeScalar(-1.0), RangeScalar(+1.0), &*v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v1 =\n" << describe(*v1,verbLevel,li,is);
      
        oss <<endl<<li<<is<< "v2 = Op*v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v2 = createMember(domain);
        op.applyTranspose(CONJ_ELE,*v1,&*v2);
        if(dump_all()) oss <<endl<<li<<is<< "v2 =\n" << describe(*v2,verbLevel,li,is);

        oss <<endl<<li<<is<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<RangeScalar> > v3 = createMember(range);
        SolveCriteria<Scalar> solveCriteria(SOLVE_TOL_REL_SOLUTION_ERR_NORM,adjoint_solution_error_solve_tol());
        assign(&*v3,RST::zero());
        SolveStatus<Scalar> solveStatus = solveTranspose(op,CONJ_ELE,*v2,&*v3,&solveCriteria);
        if(dump_all()) oss <<endl<<li<<is<< "v3 =\n" << describe(*v3,verbLevel,li,is);
        oss
          <<endl<<li<<is<< "solve status:"
          <<endl<<li<<is<<is<< "solveStatus = " << toString(solveStatus.solveStatus)
          <<endl<<li<<is<<is<< "achievedTol = " << SolveStatus<DomainScalar>::achievedTolToString(solveStatus.achievedTol)
          <<endl<<li<<is<<is<< "iterations  = " << solveStatus.iterations
          <<endl<<li<<is<<is<< "message     = \"" << solveStatus.message << "\""
          <<endl;
        oss
          <<endl<<li<<is<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(solveStatus.solveStatus==SOLVE_STATUS_CONVERGED)<<endl;
        if(solveStatus.achievedTol==SolveStatus<DomainScalar>::unknownTolerance())
          oss <<endl<<li<<is<<"achievedTol==unknownTolerance(): Setting achievedTol = adjoint_solution_error_solve_tol() = "<<adjoint_solution_error_solve_tol()<<endl;
      
        oss <<endl<<li<<is<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RefCountPtr<VectorBase<DomainScalar> > v4 = createMember(range);
        V_VmV( &*v4, *v3, *v1 );
        if(dump_all()) oss <<endl<<li<<is<< "v4 =\n" << describe(*v4,verbLevel,li,is);
      
        const DomainScalarMag
          norm_v1 = norm(*v1),
          norm_v4 = norm(*v4),
          norm_v4_norm_v1 = norm_v4/norm_v1;

        result = testMaxErr(
          "norm(v4)/norm(v1)", norm_v4_norm_v1
          ,"adjoint_solution_error_solve_tol()+adjoint_solution_error_slack_error_tol()", DomainScalarMag(adjoint_solution_error_solve_tol()+adjoint_solution_error_slack_error_tol())
          ,"adjoint_solution_error_solve_tol()+adjoint_solution_error_slack_warning_tol()", DomainScalarMag(adjoint_solution_error_solve_tol()+adjoint_solution_error_slack_warning_tol())
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
    if(out) *out <<endl<<li<< "this->check_adjoint_solution_error()==false: Skipping the check of the adjoint solve with a tolerance on the solution_error!\n";
  }
  
  if(out) {
    if(success)
      *out <<endl<<li<<"Congratulations, this LinearOpWithSolveBase object seems to check out!\n";
    else
      *out <<endl<<li<<"Oh no, at least one of the tests performed with this LinearOpWithSolveBase object failed (see above failures)!\n";
    *out <<endl<<li<< "*** Leaving LinearOpWithSolveTester<"<<ST::name()<<">::check(...)\n";
  }
  
  return success;
  
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_TESTER_HPP
