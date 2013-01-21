// $Id$ 
// $Source$ 

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
  throw "NOX Error";
}

#endif // WITH_PRERELEASE



