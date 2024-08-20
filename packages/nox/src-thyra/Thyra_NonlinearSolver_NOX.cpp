// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_NonlinearSolver_NOX.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "NOX.H"
#include "NOX_Thyra.H"
#include "NOX_PrePostOperator_Vector.H"
#include "NOX_PrePostOperator_RowSumScaling.H"
#include "NOX_MeritFunction_Weighted.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include <typeinfo>

// ****************************************************************
// ****************************************************************
Thyra::NOXNonlinearSolver::NOXNonlinearSolver():
  do_row_sum_scaling_(false),
  when_to_update_(NOX::RowSumScaling::UpdateInvRowSumVectorAtBeginningOfSolve),
  rebuild_solver_(true),
  updatePreconditioner_(true),
  use_base_point_(false)
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
  //TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(p));
  rebuild_solver_ = true;
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
  //TEUCHOS_TEST_FOR_EXCEPT(model.get()==NULL);
  rebuild_solver_ = true;
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
  use_base_point_ = true;
  basePoint_ = modelInArgs;
}

// ****************************************************************
// ****************************************************************
void 
Thyra::NOXNonlinearSolver::
setPrecOp(const Teuchos::RCP< ::Thyra::PreconditionerBase<double>>& precOp,
          const Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double>>& precFactory,
          const bool updatePreconditioner)
{
  TEUCHOS_TEST_FOR_EXCEPTION(nonnull(solver_),std::runtime_error,
                             "ERROR: Preconditioners must be set on Thyra::NonlinearSolver::NOX object before solver is constructed!");
  precOp_ = precOp;
  precFactory_ = precFactory;
  updatePreconditioner_ = updatePreconditioner;
}

// ****************************************************************
// ****************************************************************
void 
Thyra::NOXNonlinearSolver::
setGroup(const Teuchos::RCP<NOX::Abstract::Group>& group)
{
  user_defined_nox_group_ = group;
}

// ****************************************************************
// ****************************************************************
Thyra::SolveStatus<double> Thyra::NOXNonlinearSolver::
solve(VectorBase<double> *x,
      const SolveCriteria<double> * /* solveCriteria */,
      VectorBase<double> *delta)
{

#ifdef ENABLE_NOX_THYRA_TIMERS
  TEUCHOS_FUNC_TIME_MONITOR("Thyra::NOXNonlinearSolver::solve");
#endif

  TEUCHOS_ASSERT(nonnull(model_));
  TEUCHOS_ASSERT(nonnull(param_list_));

  NOX::Thyra::Vector initial_guess(Teuchos::rcp(x, false));  // View of x

  if (Teuchos::is_null(solver_) || rebuild_solver_) {

    rebuild_solver_ = false;

    this->validateAndParseThyraGroupOptions(param_list_->sublist("Thyra Group Options"));

    if (function_scaling_ != "None") {

      if (function_scaling_ == "Row Sum")
	this->setupRowSumScalingObjects();
      
      TEUCHOS_ASSERT(nonnull(scaling_vector_));
    }
    else {
      TEUCHOS_ASSERT(is_null(scaling_vector_));
    }

    right_scaling_vector_     = Teuchos::null;
    inv_right_scaling_vector_ = Teuchos::null;
    if(param_list_->isParameter("Right Scaling Vector")){
      Teuchos::RCP<NOX::Abstract::Vector> abstract_vec = param_list_->get<RCP<NOX::Abstract::Vector> >("Right Scaling Vector");
      right_scaling_vector_ = Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(abstract_vec)->getThyraRCPVector();
    }

    if(param_list_->isParameter("Inverse Right Scaling Vector")){
      Teuchos::RCP<NOX::Abstract::Vector> abstract_inv_vec = param_list_->get<RCP<NOX::Abstract::Vector> >("Inverse Right Scaling Vector");
      inv_right_scaling_vector_ = Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(abstract_inv_vec)->getThyraRCPVector();
    }

    if (is_null(user_defined_nox_group_)) {
      if (is_null(precOp_))
        nox_group_ = Teuchos::rcp(new NOX::Thyra::Group(initial_guess, model_, scaling_vector_, right_scaling_vector_,  inv_right_scaling_vector_, rightScalingFirst_));
      else {
        auto lowsFactory = model_->get_W_factory();
        auto linOp = model_->create_W_op();
        nox_group_ = Teuchos::rcp(new NOX::Thyra::Group(initial_guess, model_, linOp, lowsFactory, precOp_, precFactory_, scaling_vector_, right_scaling_vector_, inv_right_scaling_vector_, rightScalingFirst_, updatePreconditioner_));
      }
    }
    else
      nox_group_ = user_defined_nox_group_;

    if (use_base_point_)
      this->getThyraGroupNonConst(nox_group_)->setBasePoint(this->basePoint_);

    status_test_ = this->buildStatusTests(*param_list_);
    solver_ = NOX::Solver::buildSolver(nox_group_, status_test_, param_list_);
  }
  else {
    if (use_base_point_)
      this->getThyraGroupNonConst(nox_group_)->setBasePoint(this->basePoint_);

    const auto thyra_group = this->getThyraGroupConst(solver_->getSolutionGroupPtr());
    auto nonconst_thyra_group = Teuchos::rcp_const_cast<NOX::Thyra::Group>(thyra_group);
    nonconst_thyra_group->setX(initial_guess);

    solver_->reset();
  }

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

  if (is_null(t_status.extraParameters))
    t_status.extraParameters = Teuchos::parameterList("NOX Solver Status");

  t_status.extraParameters->set("Number of Iterations",
                solver_->getNumIterations());

  // Get the solution and update
  const Teuchos::RCP<const NOX::Abstract::Group> final_group = solver_->getSolutionGroupPtr();
  const Teuchos::RCP<const NOX::Thyra::Group> final_thyra_group = this->getThyraGroupConst(final_group);
  const ::Thyra::VectorBase<double>& new_x = *(final_thyra_group->get_current_x());

  if (delta)
    ::Thyra::V_StVpStV<double>(Teuchos::ptr(delta),1.0,new_x,-1.0,*x);

  ::Thyra::assign(Teuchos::ptr(x), new_x);

  return t_status;

}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const Thyra::VectorBase<double> >
Thyra::NOXNonlinearSolver::get_current_x() const
{
  return this->getThyraGroupConst(nox_group_)->get_current_x();
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
  return this->getThyraGroupNonConst(nox_group_)->getNonconstJacobian();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const Thyra::LinearOpWithSolveBase<double> >
Thyra::NOXNonlinearSolver::get_W() const
{
  return this->getThyraGroupConst(nox_group_)->getJacobian();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP< Thyra::LinearOpBase<double> >
Thyra::NOXNonlinearSolver::get_nonconst_W_op(const bool forceUpToDate)
{
  if (forceUpToDate && !nox_group_->isJacobian())
    nox_group_->computeJacobian();
  return this->getThyraGroupNonConst(nox_group_)->getNonconstJacobianOperator();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP< const Thyra::LinearOpBase<double> >
Thyra::NOXNonlinearSolver::get_W_op() const
{
  return this->getThyraGroupConst(nox_group_)->getJacobianOperator();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const Thyra::PreconditionerBase<double> >
Thyra::NOXNonlinearSolver::get_prec_op() const
{
  return this->getThyraGroupConst(nox_group_)->getPreconditioner();
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP< Thyra::PreconditionerBase<double> >
Thyra::NOXNonlinearSolver::get_nonconst_prec_op()
{
  return this->getThyraGroupNonConst(nox_group_)->getNonconstPreconditioner();
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
    validParams.set(
      "Function Scaling", "None",
      "Determines if function scaling of residual, Jacobian, etc. should be used.",
      rcp(new Teuchos::StringValidator(
          Teuchos::tuple<std::string>("None", "Row Sum", "User Defined"))));

    validParams.set(
      "Update Row Sum Scaling", "Before Each Nonlinear Solve",
      "Determines if function scaling of residual, Jacobian, etc. should be used.",
      rcp(new Teuchos::StringValidator(
          Teuchos::tuple<std::string>("Before Each Nonlinear Solve","Before Each Nonlinear Iteration"))));

    validParams.set<Teuchos::RCP< ::Thyra::VectorBase<double> > >("User Defined Scaling", Teuchos::null);
    validParams.set<bool >("Do Right Scaling First", false);
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

  rightScalingFirst_ = thyra_group_options_sublist.get<bool>("Do Right Scaling First");

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
  // ETP 5/23/16 -- Commenting this out because the parameter list may have
  // been reused from previous solves, and so the merit function that has been
  // set below from previous solves may still be there.
  //TEUCHOS_ASSERT( !(nox_parameters.sublist("Solver Options").isType<RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function")));

  RCP<NOX::MeritFunction::Generic> mf = rcp(new NOX::Thyra::WeightedMeritFunction(scaling_vector_));

  nox_parameters.sublist("Solver Options").set<RCP<NOX::MeritFunction::Generic> >("User Defined Merit Function",mf);

}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<NOX::Thyra::Group>
Thyra::NOXNonlinearSolver::
getThyraGroupNonConst(const Teuchos::RCP<NOX::Abstract::Group>& nox_group)
{
  Teuchos::RCP<NOX::Thyra::Group> ntg = Teuchos::rcp_dynamic_cast<NOX::Thyra::Group>(nox_group,false);
  if (nonnull(ntg))
    return ntg;

  auto nested_group  = nox_group->getNestedGroup();
  if (nonnull(nested_group))
    ntg = this->getThyraGroupNonConst(nested_group); // for recursively nested groups

  TEUCHOS_TEST_FOR_EXCEPTION(ntg.is_null(),std::runtime_error,
                             "ERROR: Thyra::NOXNonlinearSolver::getThyraGroupNonConst() failed to cast to a NOX::Thyra::Group. The object is of type \""
                             << Teuchos::demangleName(typeid(*nox_group).name()) << "\".\n");

  return ntg;
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<const NOX::Thyra::Group>
Thyra::NOXNonlinearSolver::
getThyraGroupConst(const Teuchos::RCP<const NOX::Abstract::Group>& nox_group) const
{
  Teuchos::RCP<const NOX::Thyra::Group> ntg = Teuchos::rcp_dynamic_cast<const NOX::Thyra::Group>(nox_group,false);
  if (nonnull(ntg))
    return ntg;

  const auto nested_group  = nox_group->getNestedGroup();
  if (nonnull(nested_group))
    ntg = this->getThyraGroupConst(nested_group); // for recursively nested groups

  TEUCHOS_TEST_FOR_EXCEPTION(ntg.is_null(),std::runtime_error,
                             "ERROR: Thyra::NOXNonlinearSolver::getThyraGroupConst() failed to cast to a NOX::Thyra::Group. The object is of type \"" << Teuchos::demangleName(typeid(*nox_group).name()) << "\".\n");

  return ntg;
}

// ****************************************************************
// ****************************************************************
