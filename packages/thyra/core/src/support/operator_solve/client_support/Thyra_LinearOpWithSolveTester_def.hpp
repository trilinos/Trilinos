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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER


#ifndef THYRA_LINEAR_OP_WITH_SOLVE_TESTER_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_TESTER_HPP


#include "Thyra_LinearOpWithSolveTester_decl.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_describeLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TestingHelpers.hpp"

#include <functional>

namespace Thyra {


// Constructors/initializers


template<class Scalar>
LinearOpWithSolveTester<Scalar>::LinearOpWithSolveTester()
  :check_forward_default_(check_forward_default_default_),
   forward_default_residual_warning_tol_(warning_tol_default_),
   forward_default_residual_error_tol_(error_tol_default_),
   forward_default_solution_error_warning_tol_(warning_tol_default_),
   forward_default_solution_error_error_tol_(error_tol_default_),
   check_forward_residual_(check_forward_residual_default_),
   forward_residual_solve_tol_(solve_tol_default_),
   forward_residual_slack_warning_tol_(slack_warning_tol_default_),
   forward_residual_slack_error_tol_(slack_error_tol_default_),
   check_adjoint_default_(check_adjoint_default_default_),
   adjoint_default_residual_warning_tol_(warning_tol_default_),
   adjoint_default_residual_error_tol_(error_tol_default_),
   adjoint_default_solution_error_warning_tol_(warning_tol_default_),
   adjoint_default_solution_error_error_tol_(error_tol_default_),
   check_adjoint_residual_(check_adjoint_residual_default_),
   adjoint_residual_solve_tol_(solve_tol_default_),
   adjoint_residual_slack_warning_tol_(slack_warning_tol_default_),
   adjoint_residual_slack_error_tol_(slack_error_tol_default_),
   num_random_vectors_(num_random_vectors_default_),
   show_all_tests_(show_all_tests_default_),
   dump_all_(dump_all_default_),
   num_rhs_(num_rhs_default_)
{}


template<class Scalar>
void LinearOpWithSolveTester<Scalar>::turn_off_all_tests()
{
  check_forward_default_ = false;
  check_forward_residual_ = false;
  check_adjoint_default_ = false;
  check_adjoint_residual_ = false;
}


template<class Scalar>
void
LinearOpWithSolveTester<Scalar>::set_all_solve_tol(
  const ScalarMag solve_tol )
{
  forward_residual_solve_tol_ = solve_tol;
  forward_residual_solve_tol_ = solve_tol;
  adjoint_residual_solve_tol_ = solve_tol;
}


template<class Scalar>
void
LinearOpWithSolveTester<Scalar>::set_all_slack_warning_tol(
  const ScalarMag slack_warning_tol )
{
  forward_default_residual_warning_tol_ = slack_warning_tol;
  forward_default_solution_error_warning_tol_ = slack_warning_tol;
  forward_residual_slack_warning_tol_ = slack_warning_tol;
  adjoint_default_residual_warning_tol_ = slack_warning_tol;
  adjoint_default_solution_error_warning_tol_ = slack_warning_tol;
  adjoint_residual_slack_warning_tol_ = slack_warning_tol;
}


template<class Scalar>
void
LinearOpWithSolveTester<Scalar>::set_all_slack_error_tol(
  const ScalarMag slack_error_tol )
{
  forward_default_residual_error_tol_ = slack_error_tol;
  forward_default_solution_error_error_tol_ = slack_error_tol;
  forward_residual_slack_error_tol_ = slack_error_tol;
  adjoint_default_residual_error_tol_ = slack_error_tol;
  adjoint_default_solution_error_error_tol_ = slack_error_tol;
  adjoint_residual_slack_error_tol_ = slack_error_tol;
}


// Overridden from ParameterListAcceptor


template<class Scalar>
void LinearOpWithSolveTester<Scalar>::setParameterList(
  const RCP<ParameterList>& paramList )
{
  using Teuchos::getParameter;
  ParameterList &pl = *paramList;
  this->setMyParamList(paramList);
  paramList->validateParametersAndSetDefaults(*getValidParameters());
  set_all_solve_tol(getParameter<ScalarMag>(pl, AllSolveTol_name_));
  set_all_slack_warning_tol(getParameter<ScalarMag>(pl, AllSlackWarningTol_name_));
  set_all_slack_error_tol(getParameter<ScalarMag>(pl, AllSlackErrorTol_name_));
  show_all_tests(getParameter<bool>(pl, ShowAllTests_name_));
  dump_all(getParameter<bool>(pl, DumpAll_name_));
  // ToDo: Add more parameters as you need them
}


template<class Scalar>
RCP<const ParameterList>
LinearOpWithSolveTester<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(AllSolveTol_name_, solve_tol_default_,
      "Sets all of the solve tolerances to the same value.  Note: This option\n"
      "is applied before any specific tolerance which may override it.");
    pl->set(AllSlackWarningTol_name_, slack_warning_tol_default_,
      "Sets all of the slack warning values to the same value.  Note: This option\n"
      "is applied before any specific tolerance which may override it.");
    pl->set(AllSlackErrorTol_name_, slack_error_tol_default_,
      "Sets all of the slack error values to the same value.  Note: This option\n"
      "is applied before any specific tolerance which may override it.");
    pl->set(ShowAllTests_name_, false,
      "If true, then all tests be traced to the output stream.");
    pl->set(DumpAll_name_, false,
      "If true, then all qualtities will be dumped to the output stream.  Warning!\n"
      "only do this to debug smaller problems as this can create a lot of output");
    // ToDo: Add more parameters as you need them (don't forget to test them)
    validPL = pl;
  }
  return validPL;
}


// LOWS testing


template<class Scalar>
bool LinearOpWithSolveTester<Scalar>::check(
  const LinearOpWithSolveBase<Scalar> &op,
  Teuchos::FancyOStream *out_arg ) const
{

  using std::endl;
  using Teuchos::as;
  using Teuchos::optInArg;
  using Teuchos::inoutArg;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  typedef Teuchos::VerboseObjectTempState<LinearOpWithSolveBase<Scalar> > VOTS;
  typedef Teuchos::ScalarTraits<Scalar> ST;

  bool success = true, result;
  const int l_num_rhs = this->num_rhs();
  Teuchos::RCP<FancyOStream> out = Teuchos::rcp(out_arg,false);
  const Teuchos::EVerbosityLevel verbLevel = (dump_all()?Teuchos::VERB_EXTREME:Teuchos::VERB_MEDIUM);
  Teuchos::Time timer("");

  OSTab tab(out,1,"THYRA");

  if(out.get()) {
    *out <<endl<< "*** Entering LinearOpWithSolveTester<"<<ST::name()<<">::check(op,...) ...\n";
    if(show_all_tests()) {
      *out <<endl<< "describe forward op:\n" << Teuchos::describe(op,verbLevel);
      if(opSupported(op, CONJTRANS) && verbLevel==Teuchos::VERB_EXTREME) {
        *out <<endl<< "describe adjoint op:\n";
        describeLinearOp<Scalar>(
          *adjoint(RCP<const LinearOpBase<Scalar> >(Teuchos::rcp(&op,false))),
          *out, verbLevel
          );
      }
    }
    else {
      *out <<endl<< "describe op: " << op.description() << endl;
    }
  }
  
  Teuchos::RCP<const VectorSpaceBase<Scalar> >   range  = op.range();
  Teuchos::RCP<const VectorSpaceBase<Scalar> >  domain = op.domain();
  
  if( check_forward_default() ) {

    if(out.get())	*out <<endl<< "this->check_forward_default()==true: Checking the default forward solve ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    result = true;
    TEUCHOS_TEST_EQUALITY_CONST(op.solveSupports(NOTRANS), true, *testOut, result);
    if (!result) these_results = false;

    if(result) {
    
      *testOut
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

        OSTab tab2(testOut);

        *testOut <<endl<< "Random vector tests = " << rand_vec_i << endl;

        tab.incrTab();
      
        *testOut <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v1 = createMembers(domain,l_num_rhs);
        Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), v1.ptr() );
        if(dump_all()) *testOut <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *testOut <<endl<< "v2 = Op*v1 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v2 = createMembers(range,l_num_rhs);
        timer.start(true);
        Thyra::apply(op, NOTRANS, *v1, v2.ptr());
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *testOut <<endl<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v3 = createMembers(domain,l_num_rhs);
        assign(v3.ptr(), ST::zero());
        SolveStatus<Scalar> solveStatus;
        {
          VOTS lowsTempState(Teuchos::rcp(&op,false),testOut,verbLevel);
          timer.start(true);
          solveStatus = solve<Scalar>(op, NOTRANS, *v2, v3.ptr());
          timer.stop();
          OSTab(testOut).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *testOut <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *testOut
          <<endl<< "solve status:\n";
        OSTab(testOut).o() << solveStatus;

        *testOut <<endl<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v4 = createMembers(domain,l_num_rhs);
        V_VmV( v4.ptr(), *v3, *v1 );
        if(dump_all()) *testOut <<endl<< "v4 =\n" << describe(*v4,verbLevel);
      
        *testOut <<endl<< "v5 = Op*v3 - v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v5 = createMembers(range, l_num_rhs);
        assign( v5.ptr(), *v2 );
        timer.start(true);
        Thyra::apply(op, NOTRANS, *v3, v5.ptr(), as<Scalar>(1.0) ,as<Scalar>(-1.0));
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v5 =\n" << describe(*v5,verbLevel);

        Array<ScalarMag> norms_v1(l_num_rhs), norms_v4(l_num_rhs), norms_v4_norms_v1(l_num_rhs);
        norms(*v1, norms_v1());
        norms(*v4, norms_v4());
        std::transform(
          norms_v4.begin(), norms_v4.end(), norms_v1.begin(), norms_v4_norms_v1.begin()
          ,std::divides<ScalarMag>()
          );

        result = testMaxErrors<Scalar>(
          "norm(v4)/norm(v1)", norms_v4_norms_v1(),
          "forward_default_solution_error_error_tol()",
          forward_default_solution_error_error_tol(),
          "forward_default_solution_error_warning_tol()",
          forward_default_solution_error_warning_tol(),
          testOut.ptr()
          );
        if(!result) these_results = false;

        Array<ScalarMag> norms_v2(l_num_rhs), norms_v5(l_num_rhs), norms_v5_norms_v2(l_num_rhs);
        norms(*v2, norms_v2());
        norms(*v5, norms_v5());
        std::transform(
          norms_v5.begin(), norms_v5.end(), norms_v2.begin(), norms_v5_norms_v2.begin()
          ,std::divides<ScalarMag>()
          );

        result = testMaxErrors<Scalar>(
          "norm(v5)/norm(v2)", norms_v5_norms_v2(),
          "forward_default_residual_error_tol()",
          forward_default_residual_error_tol(),
          "forward_default_residual_warning_tol()",
          forward_default_residual_warning_tol(),
          testOut.ptr()
          );
        if(!result) these_results = false;
        
      }
    }
    else {
      *testOut <<endl<< "Forward operator not supported, skipping check!\n";
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  }
  else {
    if(out.get()) *out <<endl<< "this->check_forward_default()==false: Skipping the check of the default forward solve!\n";
  }
  
  if( check_forward_residual() ) {

    if(out.get())	*out <<endl<< "this->check_forward_residual()==true: Checking the forward solve with a tolerance on the residual ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    result = true;
    TEUCHOS_TEST_EQUALITY_CONST(op.solveSupports(NOTRANS), true, *testOut, result);
    if (!result) these_results = false;

    if(result) {
    
      *testOut
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

        OSTab tab2(testOut);

        *testOut <<endl<< "Random vector tests = " << rand_vec_i << endl;

        tab.incrTab();
      
        *testOut <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v1 = createMembers(domain,l_num_rhs);
        Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), v1.ptr() );
        if(dump_all()) *testOut <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *testOut <<endl<< "v2 = Op*v1 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v2 = createMembers(range,l_num_rhs);
        timer.start(true);
        Thyra::apply(op, NOTRANS, *v1, v2.ptr());
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *testOut <<endl<< "v3 = inv(Op)*v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v3 = createMembers(domain,l_num_rhs);
        SolveCriteria<Scalar> solveCriteria(
          SolveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
          ,forward_residual_solve_tol()
          );
        assign(v3.ptr(), ST::zero());
        SolveStatus<Scalar> solveStatus;
        {
          VOTS lowsTempState(Teuchos::rcp(&op,false),testOut,verbLevel);
          timer.start(true);
          solveStatus = solve<Scalar>(op, NOTRANS, *v2, v3.ptr(),
            optInArg(solveCriteria));
          timer.stop();
          OSTab(testOut).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *testOut <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *testOut
          <<endl<< "solve status:\n";
        OSTab(testOut).o() << solveStatus;
        result = solveStatus.solveStatus==SOLVE_STATUS_CONVERGED;
        if(!result) these_results = false;
        *testOut
          <<endl<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(result)<<endl;
      
        *testOut <<endl<< "v4 = Op*v3 - v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v4 = createMembers(range,l_num_rhs);
        assign( v4.ptr(), *v2 );
        timer.start(true);
        Thyra::apply(op, NOTRANS, *v3, v4.ptr(), as<Scalar>(1.0), as<Scalar>(-1.0));
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v4 =\n" << describe(*v4,verbLevel);

        Array<ScalarMag> norms_v2(l_num_rhs), norms_v4(l_num_rhs), norms_v4_norms_v2(l_num_rhs);
        norms(*v2, norms_v2());
        norms(*v4, norms_v4());
        std::transform(
          norms_v4.begin(),norms_v4.end(),norms_v2.begin(),norms_v4_norms_v2.begin()
          ,std::divides<ScalarMag>()
          );

        result = testMaxErrors<Scalar>(
          "norm(v4)/norm(v2)", norms_v4_norms_v2(),
          "forward_residual_solve_tol()+forward_residual_slack_error_tol()",
          as<ScalarMag>(forward_residual_solve_tol()+forward_residual_slack_error_tol()),
          "forward_residual_solve_tol()_slack_warning_tol()",
          as<ScalarMag>(forward_residual_solve_tol()+forward_residual_slack_warning_tol()),
          testOut.ptr()
          );
        if(!result) these_results = false;

      }
    }
    else {
      *testOut <<endl<< "Forward operator not supported, skipping check!\n";
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  }
  else {
    if(out.get()) *out <<endl<< "this->check_forward_residual()==false: Skipping the check of the forward solve with a tolerance on the residual!\n";
  }
  
  if( check_adjoint_default() ) {

    if(out.get())	*out <<endl<< "this->check_adjoint_default()==true: Checking the default adjoint solve ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    result = true;
    TEUCHOS_TEST_EQUALITY_CONST(op.solveSupports(CONJTRANS), true, *testOut, result);
    if (!result) these_results = false;

    if(result) {
    
      *testOut
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

        OSTab tab2(testOut);

        *testOut <<endl<< "Random vector tests = " << rand_vec_i << endl;
        
        tab.incrTab();
      
        *testOut <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v1 = createMembers(range,l_num_rhs);
        Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), v1.ptr() );
        if(dump_all()) *testOut <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *testOut <<endl<< "v2 = Op'*v1 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v2 = createMembers(domain,l_num_rhs);
        timer.start(true);
        Thyra::apply(op, CONJTRANS, *v1, v2.ptr());
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *testOut <<endl<< "v3 = inv(Op')*v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v3 = createMembers(range,l_num_rhs);
        assign(v3.ptr(), ST::zero());
        SolveStatus<Scalar> solveStatus;
        {
          VOTS lowsTempState(Teuchos::rcp(&op,false),testOut,verbLevel);
          timer.start(true);
          solveStatus = solve<Scalar>(op, CONJTRANS, *v2, v3.ptr());
          timer.stop();
          OSTab(testOut).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *testOut <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *testOut
          <<endl<< "solve status:\n";
        OSTab(testOut).o() << solveStatus;

        *testOut <<endl<< "v4 = v3 - v1 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v4 = createMembers(range,l_num_rhs);
        V_VmV( v4.ptr(), *v3, *v1 );
        if(dump_all()) *testOut <<endl<< "v4 =\n" << describe(*v4,verbLevel);
      
        *testOut <<endl<< "v5 = Op'*v3 - v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v5 = createMembers(domain,l_num_rhs);
        assign( v5.ptr(), *v2 );
        timer.start(true);
        Thyra::apply(op, CONJTRANS, *v3, v5.ptr(), as<Scalar>(1.0), as<Scalar>(-1.0));
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v5 =\n" << describe(*v5,verbLevel);

        Array<ScalarMag> norms_v1(l_num_rhs), norms_v4(l_num_rhs), norms_v4_norms_v1(l_num_rhs);
        norms(*v1, norms_v1());
        norms(*v4, norms_v4());
        std::transform(
          norms_v4.begin(),norms_v4.end(),norms_v1.begin(),norms_v4_norms_v1.begin()
          ,std::divides<ScalarMag>()
          );

        result = testMaxErrors<Scalar>(
          "norm(v4)/norm(v1)", norms_v4_norms_v1(),
          "adjoint_default_solution_error_error_tol()",
          adjoint_default_solution_error_error_tol(),
          "adjoint_default_solution_error_warning_tol()",
          adjoint_default_solution_error_warning_tol(),
          testOut.ptr()
          );
        if(!result) these_results = false;

        Array<ScalarMag> norms_v2(l_num_rhs), norms_v5(l_num_rhs), norms_v5_norms_v2(l_num_rhs);
        norms(*v2, norms_v2());
        norms(*v5, norms_v5());
        std::transform(
          norms_v5.begin(), norms_v5.end(), norms_v2.begin(), norms_v5_norms_v2.begin()
          ,std::divides<ScalarMag>()
          );

        result = testMaxErrors<Scalar>(
          "norm(v5)/norm(v2)", norms_v5_norms_v2(),
          "adjoint_default_residual_error_tol()",
          adjoint_default_residual_error_tol(),
          "adjoint_default_residual_warning_tol()",
          adjoint_default_residual_warning_tol(),
          testOut.ptr()
          );
        if(!result) these_results = false;

      }
    }
    else {
      *testOut <<endl<< "Adjoint operator not supported, skipping check!\n";
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  }
  else {
    if(out.get()) *out <<endl<< "this->check_adjoint_default()==false: Skipping the check of the adjoint solve with a default tolerance!\n";
  }
  
  if( check_adjoint_residual() ) {

    if(out.get())	*out <<endl<< "this->check_adjoint_residual()==true: Checking the adjoint solve with a tolerance on the residual ... ";

    TestResultsPrinter testResultsPrinter(out, show_all_tests());
    const RCP<FancyOStream> testOut = testResultsPrinter.getTestOStream();

    bool these_results = true;

    result = true;
    TEUCHOS_TEST_EQUALITY_CONST(op.solveSupports(CONJTRANS), true, *testOut, result);
    if (!result) these_results = false;

    if(result) {
    
      *testOut
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

        OSTab tab2(testOut);

        *testOut <<endl<< "Random vector tests = " << rand_vec_i << endl;

        tab.incrTab();
      
        *testOut <<endl<< "v1 = randomize(-1,+1); ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v1 = createMembers(range,l_num_rhs);
        Thyra::randomize( as<Scalar>(-1.0), as<Scalar>(+1.0), v1.ptr() );
        if(dump_all()) *testOut <<endl<< "v1 =\n" << describe(*v1,verbLevel);
      
        *testOut <<endl<< "v2 = Op'*v1 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v2 = createMembers(domain,l_num_rhs);
        timer.start(true);
        Thyra::apply(op, CONJTRANS, *v1, v2.ptr());
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v2 =\n" << describe(*v2,verbLevel);

        *testOut <<endl<< "v3 = inv(Op')*v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v3 = createMembers(range,l_num_rhs);
        SolveCriteria<Scalar> solveCriteria(
          SolveMeasureType(SOLVE_MEASURE_NORM_RESIDUAL,SOLVE_MEASURE_NORM_RHS)
          ,adjoint_residual_solve_tol()
          );
        assign(v3.ptr(), ST::zero());
        SolveStatus<Scalar> solveStatus;
        {
          VOTS lowsTempState(Teuchos::rcp(&op,false),testOut,verbLevel);
          timer.start(true);
          solveStatus = solve<Scalar>(op, CONJTRANS, *v2, v3.ptr(), optInArg(solveCriteria));
          timer.stop();
          OSTab(testOut).o() <<"\n=> Solve time = "<<timer.totalElapsedTime()<<" sec\n";
        }
        if(dump_all()) *testOut <<endl<< "v3 =\n" << describe(*v3,verbLevel);
        *testOut
          <<endl<< "solve status:\n";
        OSTab(testOut).o() << solveStatus;
        result = solveStatus.solveStatus==SOLVE_STATUS_CONVERGED;
        if(!result) these_results = false;
        *testOut
          <<endl<< "check: solveStatus = " << toString(solveStatus.solveStatus) << " == SOLVE_STATUS_CONVERGED : "
          << passfail(result)<<endl;
        if(solveStatus.achievedTol==SolveStatus<Scalar>::unknownTolerance())
          *testOut <<endl<<"achievedTol==unknownTolerance(): Setting achievedTol = adjoint_residual_solve_tol() = "<<adjoint_residual_solve_tol()<<endl;
      
        *testOut <<endl<< "v4 = Op'*v3 - v2 ...\n" ;
        Teuchos::RCP<MultiVectorBase<Scalar> > v4 = createMembers(domain,l_num_rhs);
        assign( v4.ptr(), *v2 );
        timer.start(true);
        Thyra::apply(op, CONJTRANS, *v3, v4.ptr(), as<Scalar>(1.0), as<Scalar>(-1.0));
        timer.stop();
        OSTab(testOut).o() <<"\n=> Apply time = "<<timer.totalElapsedTime()<<" sec\n";
        if(dump_all()) *testOut <<endl<< "v4 =\n" << describe(*v4,verbLevel);

        Array<ScalarMag> norms_v2(l_num_rhs), norms_v4(l_num_rhs), norms_v4_norms_v2(l_num_rhs);
        norms(*v2, norms_v2());
        norms(*v4, norms_v4());
        std::transform(
          norms_v4.begin(),norms_v4.end(),norms_v2.begin(),norms_v4_norms_v2.begin()
          ,std::divides<ScalarMag>()
          );

        result = testMaxErrors<Scalar>(
          "norm(v4)/norm(v2)", norms_v4_norms_v2(),
          "adjoint_residual_solve_tol()+adjoint_residual_slack_error_tol()",
          as<ScalarMag>(adjoint_residual_solve_tol()+adjoint_residual_slack_error_tol()),
          "adjoint_residual_solve_tol()_slack_warning_tol()",
          as<ScalarMag>(adjoint_residual_solve_tol()+adjoint_residual_slack_warning_tol()),
          testOut.ptr()
          );
        if(!result) these_results = false;

      }
    }
    else {
      *testOut <<endl<< "Adjoint operator not supported, skipping check!\n";
    }

    testResultsPrinter.printTestResults(these_results, inoutArg(success));

  }
  else {
    if(out.get())
      *out <<endl<< "this->check_adjoint_residual()==false: Skipping the check of the adjoint solve with a tolerance on the residual!\n";
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


// private static members


// Define local macros used in the next few lines and then undefined

#define LOWST_DEFINE_RAW_STATIC_MEMBER( DATA_TYPE, MEMBER_NAME, DEFAULT_VALUE ) \
template<class Scalar> \
const DATA_TYPE \
LinearOpWithSolveTester<Scalar>::MEMBER_NAME = DEFAULT_VALUE

#define LOWST_DEFINE_MTD_STATIC_MEMBER( DATA_TYPE, MEMBER_NAME, DEFAULT_VALUE ) \
template<class Scalar> \
const typename LinearOpWithSolveTester<Scalar>::DATA_TYPE \
LinearOpWithSolveTester<Scalar>::MEMBER_NAME = DEFAULT_VALUE

LOWST_DEFINE_RAW_STATIC_MEMBER(bool, check_forward_default_default_, true);
LOWST_DEFINE_RAW_STATIC_MEMBER(bool, check_forward_residual_default_, true);
LOWST_DEFINE_RAW_STATIC_MEMBER(bool, check_adjoint_default_default_, true);
LOWST_DEFINE_RAW_STATIC_MEMBER(bool, check_adjoint_residual_default_, true);

LOWST_DEFINE_MTD_STATIC_MEMBER(ScalarMag, warning_tol_default_, 1e-6);
LOWST_DEFINE_MTD_STATIC_MEMBER(ScalarMag, error_tol_default_, 1e-5);
LOWST_DEFINE_MTD_STATIC_MEMBER(ScalarMag, solve_tol_default_, 1e-5);
LOWST_DEFINE_MTD_STATIC_MEMBER(ScalarMag, slack_warning_tol_default_, 1e-6);
LOWST_DEFINE_MTD_STATIC_MEMBER(ScalarMag, slack_error_tol_default_, 1e-5);

LOWST_DEFINE_RAW_STATIC_MEMBER(int, num_random_vectors_default_, 1);
LOWST_DEFINE_RAW_STATIC_MEMBER(bool, show_all_tests_default_, false);
LOWST_DEFINE_RAW_STATIC_MEMBER(bool, dump_all_default_, false);
LOWST_DEFINE_RAW_STATIC_MEMBER(int, num_rhs_default_, 1);

LOWST_DEFINE_RAW_STATIC_MEMBER(std::string, AllSolveTol_name_, "All Solve Tol");
LOWST_DEFINE_RAW_STATIC_MEMBER(std::string, AllSlackWarningTol_name_, "All Slack Warning Tol");
LOWST_DEFINE_RAW_STATIC_MEMBER(std::string, AllSlackErrorTol_name_, "All Slack Error Tol");
LOWST_DEFINE_RAW_STATIC_MEMBER(std::string, ShowAllTests_name_, "Show All Tests");
LOWST_DEFINE_RAW_STATIC_MEMBER(std::string, DumpAll_name_, "Dump All");

#undef LOWST_DEFINE_MTD_STATIC_MEMBER

#undef LOWST_DEFINE_RAW_STATIC_MEMBER


} // namespace Thyra


#endif // THYRA_LINEAR_OP_WITH_SOLVE_TESTER_HPP
