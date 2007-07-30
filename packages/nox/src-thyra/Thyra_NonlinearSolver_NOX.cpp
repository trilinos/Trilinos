// @HEADER
// @HEADER

#include "Thyra_NonlinearSolver_NOX.hpp"

#include "NOX.H"
#include "NOX_Thyra.H"

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
{
  
}

// ****************************************************************
// ****************************************************************
void Thyra::NOXNonlinearSolver::
setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& p)
{
  TEST_FOR_EXCEPT(Teuchos::is_null(p));
  param_list_ = p;
}

// ****************************************************************
// ****************************************************************
Teuchos::RCP<Teuchos::ParameterList> 
Thyra::NOXNonlinearSolver::getParameterList()
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
  TEST_FOR_EXCEPT(model.get()==NULL);
  model_ = model;
  J_ = Teuchos::null;
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
Thyra::SolveStatus<double> Thyra::NOXNonlinearSolver::
solve(VectorBase<double> *x,
      const SolveCriteria<double> *solveCriteria,
      VectorBase<double> *delta
      )
{
  TEST_FOR_EXCEPT(model_.get()==NULL);
  TEST_FOR_EXCEPT(param_list_.get()==NULL);
 
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = x->clone_v();
  
  Teuchos::RCP<NOX::Thyra::Group> nox_group = 
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model_));
    
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

  NOX::Solver::Manager solver(nox_group, combo, param_list_);
  NOX::StatusTest::StatusType solvStatus = solver.solve();
  
  Thyra::SolveStatus<double> t_status;

  if (solvStatus == NOX::StatusTest::Converged)
    t_status.solveStatus = SOLVE_STATUS_CONVERGED;
  else if (solvStatus == NOX::StatusTest::Unconverged) 
    t_status.solveStatus = SOLVE_STATUS_UNCONVERGED;
  else if (solvStatus == NOX::StatusTest::Failed) 
    t_status.solveStatus = SOLVE_STATUS_UNCONVERGED;
  else
    t_status.solveStatus = SOLVE_STATUS_UNCONVERGED;


  const NOX::Abstract::Group& final_group = solver.getSolutionGroup();
  const NOX::Abstract::Vector& final_solution = final_group.getX();

  const NOX::Thyra::Vector& vec = 
    dynamic_cast<const NOX::Thyra::Vector&>(final_solution);

  const ::Thyra::VectorBase<double>& new_x = 
    vec.getThyraVector();

  Thyra::V_StVpStV<double>(delta,1.0,new_x,-1.0,*x);

  *x = new_x;

  // Return default status
  return t_status;

}


Teuchos::RCP< Thyra::LinearOpWithSolveBase<double> >
Thyra::NOXNonlinearSolver::get_nonconst_W()
{
  return J_;
}


Teuchos::RCP<const Thyra::LinearOpWithSolveBase<double> >
Thyra::NOXNonlinearSolver::get_W() const
{
  return J_;
}

// ****************************************************************
// ****************************************************************
