//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorTaylor_impl_hpp
#define Tempus_PhiEvaluatorTaylor_impl_hpp

#include "Tempus_PhiEvaluatorTaylor.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluator_decl.hpp"
#include "Teuchos_RCPDecl.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorBase.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
PhiEvaluatorTaylor<Scalar>::getValidParameters() const
{
  //TODO
  Teuchos::RCP<Teuchos::ParameterList> pl = this->getValidParametersBasic();

  pl->set(
      "PhiEvaluator Type", "Taylor",
      "Method to approximate the phi-function evaluation.");

  pl->set(
      "Taylor method", "CN",
      "'PDF method' determines the partial fraction decomposition used to approximate the exponential.  "
      "'IE' - uses an implicit Euler approximation (order 1).  "
      "'CN' - uses a Crank-Nicolson approximation (order 2).");

  pl->set<std::string>("Stepper Type", "EPI",
      "The type of Tempus stepper that will use this Phi evaluator.");

  pl->set<bool>(
      "Use FSAL", true,
      "The First-Same-As-Last (FSAL) principle is the situation where the\n"
      "last function evaluation, f(x^{n-1},t^{n-1}) [a.k.a. xDot^{n-1}],\n"
      "can be used for the first function evaluation, f(x^n,t^n)\n"
      "[a.k.a. xDot^n].  For RK methods, this applies to the stages.\n"
      "\n"
      "Often the FSAL priniciple can be used to save an evaluation.\n"
      "However there are cases when it cannot be used, e.g., operator\n"
      "splitting where other steppers/operators have modified the solution,\n"
      "x^*, and thus require the function evaluation, f(x^*, t^{n-1}).\n"
      "\n"
      "It should be noted that when the FSAL priniciple can be used\n"
      "(can set useFSAL=true), setting useFSAL=false will give the\n"
      "same solution but at additional expense.  However, the reverse\n"
      "is not true.  When the FSAL priniciple can not be used\n"
      "(need to set useFSAL=false), setting useFSAL=true will produce\n"
      "incorrect solutions.\n"
      "\n"
      "Default in general for explicit and implicit steppers is false,\n"
      "but individual steppers can override this default.");

  pl->set<std::string>(
      "Initial Condition Consistency", "Consistent",
      "This indicates which type of consistency should be applied to\n"
      "the initial conditions (ICs):\n"
      "\n"
      "  'None'   - Do nothing to the ICs provided in the SolutionHistory.\n"
      "  'Zero'   - Set the derivative of the SolutionState to zero in the\n"
      "             SolutionHistory provided, e.g., xDot^0 = 0, or \n"
      "             xDotDot^0 = 0.\n"
      "  'App'    - Use the application's ICs, e.g., getNominalValues().\n"
      "  'Consistent' - Make the initial conditions for x and xDot\n"
      "             consistent with the governing equations, e.g.,\n"
      "             xDot = f(x,t), and f(x, xDot, t) = 0.  For implicit\n"
      "             ODEs, this requires a solve of f(x, xDot, t) = 0 for\n"
      "             xDot, and another Jacobian and residual may be\n"
      "             needed, e.g., boundary conditions on xDot may need\n"
      "             to replace boundary conditions on x.\n"
      "\n"
      "In general for explicit steppers, the default is 'Consistent',\n"
      "because it is fairly cheap with just one residual evaluation.\n"
      "In general for implicit steppers, the default is 'None', because\n"
      "the application often knows its IC and can set it the initial\n"
      "SolutionState.  Also, as noted above, 'Consistent' may require\n"
      "another Jacobian from the application.  Individual steppers may\n"
      "override these defaults.");

  pl->set<bool>(
      "Initial Condition Consistency Check", true,
      "Check if the initial condition, x and xDot, is consistent with the\n"
      "governing equations, xDot = f(x,t), or f(x, xDot, t) = 0.\n"
      "\n"
      "In general for explicit and implicit steppers, the default is true,\n"
      "because it is fairly cheap with just one residual evaluation.\n"
      "Individual steppers may override this default.");

  pl->set<int>(
      "Taylor Expansion Order", 2,
      "The order of the Taylor expansion used in the EPI stepper.\n"
      "\n"
      "The default is 2.");

  pl->set<std::string>(
      "Phi Evaluator Type", "Taylor",
      "The type of Phi evaluator used in the EPI stepper.\n"
      "\n"
      "The default is 'Taylor'.  Other options are 'PFD' and 'Leja'.");

  //pl->set("?", *member_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
void PhiEvaluatorTaylor<Scalar>::setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs)
{
  inArgs_lin_ = Teuchos::rcpFromRef(inArgs);
}

template <class Scalar>
Thyra::SolveStatus<Scalar> PhiEvaluatorTaylor<Scalar>::computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
							       int k, Scalar cdt, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> rhs_b)
{
  this->phiLinSolv_->computeMassMatrix(*inArgs_lin_);
  this->phiLinSolv_->computeJacobian(*inArgs_lin_);
  this->phiLinSolv_->buildK(k);
  this->phiLinSolv_->buildb(k, rhs_b);
  this->phiLinSolv_->buildATilde(cdt);
  this->phiLinSolv_->buildv();
  auto vec = this->phiLinSolv_->matrixExponential(taylorExpOrder_);

  // Get the first block of the multi-vector calculated from 2x2 multi-matrix
  auto pv = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<Scalar>>(vec, true);
  // phiv = pv->getVectorBlock(0);  // V block
  auto v0 = pv->getVectorBlock(0);  // V block (const)
  Thyra::copy(*v0, phiv.ptr());
  // return vecV;
  Thyra::SolveStatus<Scalar> sStatus;
  return sStatus;
}
  
// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar> > createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto phi = Teuchos::rcp(new PhiEvaluatorTaylor<Scalar>());
  phi->setName("From createPhiEvaluatorTaylor");

  if (pl == Teuchos::null || pl->numParams() == 0) return phi;

  pl->validateParametersAndSetDefaults(*phi->getValidParameters());

  phi->setTaylorExpansionOrder(pl->get<int>("Taylor Expansion Order", 2));
  phi->setName(pl->name());
  //phi->setThing(pl->get("Thing", "default"));

  return phi;
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorTaylor_impl_hpp
