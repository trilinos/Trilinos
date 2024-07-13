// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_Direction_ModifiedNewton.H" // class definition
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"

NOX::Direction::ModifiedNewton::
ModifiedNewton(const Teuchos::RCP<NOX::GlobalData>& gd,
           Teuchos::ParameterList& p)
{
  reset(gd, p);
  ageOfJacobian = -1;
  if (p.sublist("Modified-Newton").get("Max Age of Jacobian", 10) < 0)
    p.sublist("Modified-Newton").set("Max Age of Jacobian", 0);
}

NOX::Direction::ModifiedNewton::~ModifiedNewton()
{
}

bool NOX::Direction::ModifiedNewton::
reset(const Teuchos::RCP<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  globalDataPtr = gd;
  utils = gd->getUtils();

  paramsPtr = &params;

  Teuchos::ParameterList& p = params.sublist("Modified-Newton");

  doRescue = p.get("Rescue Bad Newton Solve", true);
  if (!p.sublist("Linear Solver").isParameter("Tolerance"))
    p.sublist("Linear Solver").get("Tolerance", 1.0e-10);
  ageOfJacobian = -1;
  return true;
}

bool NOX::Direction::ModifiedNewton::
compute(NOX::Abstract::Vector& dir,
    NOX::Abstract::Group& soln,
    const NOX::Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType status;

  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok)
    throwError("compute", "Unable to compute F");

  maxAgeOfJacobian = paramsPtr->sublist("Modified-Newton").get("Max Age of Jacobian", 10);

  if (Teuchos::is_null(oldJacobianGrpPtr)) {
    oldJacobianGrpPtr = soln.clone(DeepCopy);
  }
  NOX::Abstract::Group& oldJacobianGrp = *oldJacobianGrpPtr;

  status = NOX::Abstract::Group::Failed;
  while (status != NOX::Abstract::Group::Ok) {
    // Conditionally compute Jacobian at current solution.
    if ( (ageOfJacobian == -1) || (ageOfJacobian == maxAgeOfJacobian) ) {

      if (ageOfJacobian > 0)
        oldJacobianGrp = soln;
      status = oldJacobianGrp.computeJacobian();
      if (status != NOX::Abstract::Group::Ok)
        throwError("compute", "Unable to compute Jacobian");
      ageOfJacobian = 1;
    }
    else
      ageOfJacobian++;

    // Compute the Modified Newton direction
    status = oldJacobianGrp.applyJacobianInverse(paramsPtr->sublist("Modified-Newton").sublist("Linear Solver"), soln.getF(), dir);
    oldJacobianGrp.logLastLinearSolveStats(*globalDataPtr->getNonConstSolverStatistics());
    dir.scale(-1.0);

    // It didn't converge, but maybe we can recover.
    if ((status != NOX::Abstract::Group::Ok) &&
        (doRescue == false)) {
      throwError("compute", "Unable to solve Newton system");
    }
    else if ((status != NOX::Abstract::Group::Ok) &&
             (doRescue == true)) {
      if (utils->isPrintType(NOX::Utils::Warning))
        utils->out() << "WARNING: NOX::Direction::ModifiedNewton::compute() - "
             << "Linear solve failed to achieve convergence - "
             << "using the step anyway since \"Rescue Bad Newton Solve\" "
             << "is true. Also, flagging recompute of Jacobian." << std::endl;
      ageOfJacobian = maxAgeOfJacobian;
      status = NOX::Abstract::Group::Ok;
    }
  }

  return true;
}

bool NOX::Direction::ModifiedNewton::compute(NOX::Abstract::Vector& dir,
                                             NOX::Abstract::Group& soln,
                                             const NOX::Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
}

bool  NOX::Direction::ModifiedNewton::rescueBadNewtonSolve(const NOX::Abstract::Group& grp) const
{
  //! Check if the "rescue" option has been selected
  if (!doRescue)
    return false;

  //! See if the group has compute the accuracy
  double accuracy;
  NOX::Abstract::Group::ReturnType status = oldJacobianGrpPtr->getNormLastLinearSolveResidual(accuracy);

  // If this functionality is not supported in the group, return false
  /* NOTE FROM TAMMY: We could later modify this to acutally caluclate
     the error itself if it's just a matter of the status being
     NotDefined. */
  if (status != NOX::Abstract::Group::Ok)
    return false;

  // Check if there is any improvement in the relative residual
  double normF = grp.getNormF();

  // If we can't reduce the relative norm at all, we're not happy
  if (accuracy >= normF)
    return false;

  // Otherwise, we just print a warning and keep going
  if (utils->isPrintType(NOX::Utils::Warning))
    utils->out() << "WARNING: NOX::Direction::ModifiedNewton::compute - Unable to achieve desired linear solve accuracy." << std::endl;
  return true;

}

void NOX::Direction::ModifiedNewton::throwError(const std::string& functionName, const std::string& errorMsg)
{
  if (utils->isPrintType(NOX::Utils::Error))
    utils->err() << "NOX::Direction::ModifiedNewton::" << functionName << " - " << errorMsg << std::endl;
  throw std::runtime_error("NOX Error");
}

#endif // WITH_PRERELEASE



