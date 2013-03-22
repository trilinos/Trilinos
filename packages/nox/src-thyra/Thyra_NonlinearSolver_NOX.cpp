//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "Thyra_NonlinearSolver_NOX.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "NOX.H"
#include "NOX_Thyra.H"
#include "NOX_PrePostOperator_Vector.H"
#include "NOX_PrePostOperator_RowSumScaling.H"
#include "NOX_MeritFunction_Weighted.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// ****************************************************************
// ****************************************************************
Thyra::NOXNonlinearSolver::NOXNonlinearSolver()
{
  param_list_ = Teuchos::rcp(new Teuchos::ParameterList);
  valid_param_list_ = Teuchos::rcp(new Teuchos::ParameterList);
}

// ****************************************************************
// ****************************************************************
Thyra::NOXNonlinearSolver::~NOXNonlinearSolver()
{}

// ****************************************************************
// ****************************************************************
void Thyra::NOXNonlinearSolver::
setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& p)
{
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(p));
  param_list_ = p;
  this->resetSolver();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<Teuchos::ParameterList>
Thyra::NOXNonlinearSolver::getNonconstParameterList()
{
  return param_list_;
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<Teuchos::ParameterList>
Thyra::NOXNonlinearSolver::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> _paramList = param_list_;
  param_list_ = Teuchos::null;
  return _paramList;
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const Teuchos::ParameterList>
Thyra::NOXNonlinearSolver::getParameterList() const
{
  return param_list_;
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const Teuchos::ParameterList>
Thyra::NOXNonlinearSolver::getValidParameters() const
{
  return valid_param_list_;
}

// ****************************************************************
// ****************************************************************
void Thyra::NOXNonlinearSolver::
setModel(const Teuchos::RCP<const ModelEvaluator<double> >& model)
{
  TEUCHOS_TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  basePoint_ = model_->createInArgs();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP< const Thyra::ModelEvaluator<double> >
Thyra::NOXNonlinearSolver::getModel() const
{
  return model_;
}

// ****************************************************************
// ****************************************************************
void
Thyra::NOXNonlinearSolver::setBasePoint(
    const Thyra::ModelEvaluatorBase::InArgs<double> &modelInArgs)
{
  basePoint_ = modelInArgs;
}

// ****************************************************************
// ****************************************************************
Thyra::SolveStatus<double> Thyra::NOXNonlinearSolver::
solve(VectorBase<double> *x,
      const SolveCriteria<double> *solveCriteria,
      VectorBase<double> *delta)
{

#ifdef ENABLE_NOX_THYRA_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Thyra::NOXNonlinearSolver::solve");
#endif

  TEUCHOS_ASSERT(nonnull(model_));
  TEUCHOS_ASSERT(nonnull(param_list_));

  NOX::Thyra::Vector initial_guess(Teuchos::rcp(x, false));  // View of x

  if (Teuchos::is_null(solver_)) {


    this->validateAndParseThyraGroupOptions(param_list_->sublist("Thyra Group Options"));

    if (function_scaling_ != "None") {

      if (function_scaling_ == "Row Sum")
	this->setupRowSumScalingObjects();

      TEUCHOS_ASSERT(nonnull(scaling_vector_));
    }
    else {
      TEUCHOS_ASSERT(is_null(scaling_vector_));
    }

    nox_group_ = Teuchos::rcp(new NOX::Thyra::Group(initial_guess, model_, scaling_vector_));
    nox_group_->getNonconstInArgs() = this->basePoint_;

    status_test_ = this->buildStatusTests(*param_list_);
    solver_ = NOX::Solver::buildSolver(nox_group_, status_test_, param_list_);
  }
  else
    solver_->reset(initial_guess);

  NOX::StatusTest::StatusType solvStatus = solver_->solve();

  Thyra::SolveStatus<double> t_status;

  if (solvStatus == NOX::StatusTest::Converged)
    t_status.solveStatus = SOLVE_STATUS_CONVERGED;
  else if (solvStatus == NOX::StatusTest::Unconverged)
    t_status.solveStatus = SOLVE_STATUS_UNCONVERGED;
  else if (solvStatus == NOX::StatusTest::Failed)
    t_status.solveStatus = SOLVE_STATUS_UNCONVERGED;
  else
    t_status.solveStatus = SOLVE_STATUS_UNCONVERGED;

  // Get the solution and update
  const NOX::Abstract::Group& final_group = solver_->getSolutionGroup();
  const NOX::Abstract::Vector& final_solution = final_group.getX();

  const NOX::Thyra::Vector& vec =
    dynamic_cast<const NOX::Thyra::Vector&>(final_solution);

  const ::Thyra::VectorBase<double>& new_x =
    vec.getThyraVector();

  if (delta)
    ::Thyra::V_StVpStV<double>(Teuchos::ptr(delta),1.0,new_x,-1.0,*x);

  //*x = new_x;
  ::Thyra::assign(Teuchos::ptr(x), new_x);

  return t_status;

}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const Thyra::VectorBase<double> >
Thyra::NOXNonlinearSolver::get_current_x() const
{
  return nox_group_->get_current_x();
}

// ****************************************************************
// ****************************************************************
bool Thyra::NOXNonlinearSolver::is_W_current() const
{
  return nox_group_->isJacobian();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP< Thyra::LinearOpWithSolveBase<double> >
Thyra::NOXNonlinearSolver::get_nonconst_W(const bool forceUpToDate)
{
  if (forceUpToDate && !nox_group_->isJacobian())
    nox_group_->computeJacobian();
  return nox_group_->getNonconstJacobian();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const Thyra::LinearOpWithSolveBase<double> >
Thyra::NOXNonlinearSolver::get_W() const
{
  return nox_group_->getJacobian();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> Thyra::NOXNonlinearSolver::
buildStatusTests(Teuchos::ParameterList& p)
{
  Teuchos::RCP<NOX::StatusTest::Generic> status_test;

  NOX::Utils utils(p.sublist("Printing"));

  if (p.isSublist("Status Tests")) {
    status_test =
      NOX::StatusTest::buildStatusTests(p.sublist("Status Tests"), utils);
  }
  else { // Default status test
    Teuchos::RCP<NOX::StatusTest::NormF> absresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
      Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    Teuchos::RCP<NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(wrms);
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    status_test = combo;
  }

  return status_test;
}

// ****************************************************************
// ****************************************************************
void Thyra::NOXNonlinearSolver::resetSolver()
{
  nox_group_ = Teuchos::null;
  status_test_ = Teuchos::null;
  solver_ = Teuchos::null;
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const NOX::Solver::Generic>
Thyra::NOXNonlinearSolver::getNOXSolver() const
{
  return solver_;
}

// ****************************************************************
// ****************************************************************
void Thyra::NOXNonlinearSolver::
validateAndParseThyraGroupOptions(Teuchos::ParameterList& thyra_group_options_sublist)
{
  using Teuchos::ParameterList;

  ParameterList validParams;
  {
    Teuchos::setStringToIntegralParameter<int>(
      "Function Scaling",
      "None",
      "Determines if function scaling of residual, Jacobian, etc. should be used.",
      Teuchos::tuple<std::string>("None","Row Sum", "User Defined"),
      &validParams
      );

    Teuchos::setStringToIntegralParameter<int>(
      "Update Row Sum Scaling",
      "Before Each Nonlinear Solve",
      "Determines if function scaling of residual, Jacobian, etc. should be used.",
      Teuchos::tuple<std::string>("Before Each Nonlinear Solve","Before Each Nonlinear Iteration"),
      &validParams
      );

    validParams.set<Teuchos::RCP< ::Thyra::VectorBase<double> > >("User Defined Scaling", Teuchos::null);
  }

  thyra_group_options_sublist.validateParametersAndSetDefaults(validParams);

  function_scaling_ = thyra_group_options_sublist.get<std::string>("Function Scaling");

  if (function_scaling_ =="Row Sum")
    do_row_sum_scaling_ = true;
  else
    do_row_sum_scaling_ = false;

  std::string string_when_to_update = thyra_group_options_sublist.get<std::string>("Update Row Sum Scaling");
  if (string_when_to_update == "Before Each Nonlinear Solve")
    when_to_update_ = NOX::RowSumScaling::UpdateInvRowSumVectorAtBeginningOfSolve;
  else if (string_when_to_update == "Before Each Nonlinear Iteration")
    when_to_update_ = NOX::RowSumScaling::UpdateInvRowSumVectorAtBeginningOfIteration;

  if (function_scaling_ =="User Defined")
    scaling_vector_ = thyra_group_options_sublist.get<Teuchos::RCP< ::Thyra::VectorBase<double> > >("User Defined Scaling");

}

// ****************************************************************
// ****************************************************************
void Thyra::NOXNonlinearSolver::setupRowSumScalingObjects()
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  scaling_vector_ = ::Thyra::createMember(model_->get_f_space());

  ::Thyra::V_S(scaling_vector_.ptr(),1.0);

  RCP<NOX::Abstract::PrePostOperator> row_sum_observer =
    rcp(new NOX::RowSumScaling(scaling_vector_, when_to_update_));

  Teuchos::ParameterList& nox_parameters = *param_list_;

  if (nox_parameters.sublist("Solver Options").
      isType< RCP<NOX::Abstract::PrePostOperator> >("User Defined Pre/Post Operator")) {

    RCP<NOX::Abstract::PrePostOperator> user_observer =
      nox_parameters.sublist("Solver Options").get< RCP<NOX::Abstract::PrePostOperator> >("User Defined Pre/Post Operator");

    // NOTE: the row_sum_observer should be evalauted after any user
    // oberservers to make sure that the jacobian is accurate.  This
    // is needed, for example, if we have a model evaluator decorator
    // that adds extra input parameters to the model such as a
    // predictor or previous time step solution to be used for
    // semi-implicit models.  The row sum would accidentally use the
    // previous predicted value which would be bad.
    RCP<NOX::PrePostOperatorVector> observer_vector = Teuchos::rcp(new NOX::PrePostOperatorVector);
    observer_vector->pushBack(user_observer);
    observer_vector->pushBack(row_sum_observer);

    nox_parameters.sublist("Solver Options").set< RCP<NOX::Abstract::PrePostOperator> >("User Defined Pre/Post Operator", observer_vector);

  }
  else
    nox_parameters.sublist("Solver Options").set< RCP<NOX::Abstract::PrePostOperator> >("User Defined Pre/Post Operator", row_sum_observer);


  // Set the weighted merit function.  Throw error if a user defined
  // merit funciton is present.
  TEUCHOS_ASSERT( !(nox_parameters.sublist("Solver Options").isType<RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function")));

  RCP<NOX::MeritFunction::Generic> mf = rcp(new NOX::Thyra::WeightedMeritFunction(scaling_vector_));

  nox_parameters.sublist("Solver Options").set<RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function",mf);

}

// ****************************************************************
// ****************************************************************
