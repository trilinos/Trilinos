//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_BASIC_DISCRETE_ADJOINT_STEPPER_TESTER_DEF_H
#define Rythmos_BASIC_DISCRETE_ADJOINT_STEPPER_TESTER_DEF_H


#include "Rythmos_BasicDiscreteAdjointStepperTester_decl.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Rythmos_DefaultIntegrator.hpp" // ToDo: Remove when we can!
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"


template<class Scalar>
Teuchos::RCP<Rythmos::BasicDiscreteAdjointStepperTester<Scalar> >
Rythmos::basicDiscreteAdjointStepperTester()
{
  return Teuchos::rcp(new BasicDiscreteAdjointStepperTester<Scalar>);
}


template<class Scalar>
Teuchos::RCP<Rythmos::BasicDiscreteAdjointStepperTester<Scalar> >
Rythmos::basicDiscreteAdjointStepperTester(const RCP<ParameterList> &paramList)
{
  const RCP<Rythmos::BasicDiscreteAdjointStepperTester<Scalar> > bdast =
    basicDiscreteAdjointStepperTester<Scalar>();
  bdast->setParameterList(paramList);
  return bdast;
}


namespace Rythmos {


// Overridden from ParameterListAcceptor


template<class Scalar>
void BasicDiscreteAdjointStepperTester<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  namespace BDASTU = BasicDiscreteAdjointStepperTesterUtils;
  paramList->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(paramList);
  errorTol_ = Teuchos::getParameter<double>(*paramList, BDASTU::ErrorTol_name);
}


template<class Scalar>
RCP<const ParameterList>
BasicDiscreteAdjointStepperTester<Scalar>::getValidParameters() const
{
  namespace BDASTU = BasicDiscreteAdjointStepperTesterUtils;
  static RCP<const ParameterList> validPL;
  if (is_null(validPL)) {
    const RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(BDASTU::ErrorTol_name, BDASTU::ErrorTol_default);
    validPL = pl;
  }
  return validPL;
}


template<class Scalar>
bool BasicDiscreteAdjointStepperTester<Scalar>::testAdjointStepper(
  Thyra::ModelEvaluator<Scalar>& adjointModel,
  const Ptr<IntegratorBase<Scalar> >& forwardIntegrator
  )
{

  using Teuchos::describe;
  using Teuchos::outArg;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;
  using Teuchos::dyn_cast;
  using Teuchos::OSTab;
  using Teuchos::sublist;
  typedef Thyra::ModelEvaluatorBase MEB;
  using Thyra::createMember;
  namespace BDASTU = BasicDiscreteAdjointStepperTesterUtils;

  const RCP<Teuchos::FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  const RCP<ParameterList> paramList = this->getMyNonconstParamList();

  bool success = true;

  OSTab tab(out);

  //
  *out << "\n*** Entering BasicDiscreteAdjointStepperTester<Scalar>::testAdjointStepper(...) ...\n";
  //

  const TimeRange<Scalar> fwdTimeRange = forwardIntegrator->getFwdTimeRange();
  const RCP<Rythmos::StepperBase<Scalar> > fwdStepper =
    forwardIntegrator->getNonconstStepper();
  const RCP<const Thyra::ModelEvaluator<Scalar> > fwdModel =
    fwdStepper->getModel();
  const RCP<const Thyra::VectorSpaceBase<Scalar> > f_space = fwdModel->get_f_space();

  //
  *out << "\nA) Construct the IC basis B ...\n";
  //

  const RCP<Thyra::MultiVectorBase<Scalar> > B =
    createMembers(fwdModel->get_x_space(), 1, "B");
  Thyra::seed_randomize<Scalar>(0);
  Thyra::randomize<Scalar>( -1.0, +1.0, B.ptr() );
  
  *out << "\nB: " << describe(*B, verbLevel);

  //
  *out << "\nB) Construct the forward sensitivity integrator ...\n";
  //

  const RCP<ForwardSensitivityStepper<Scalar> > fwdSensStepper =
    forwardSensitivityStepper<Scalar>();
  fwdSensStepper->initializeSyncedSteppersInitCondOnly(
    fwdModel,
    B->domain(), // p_space
    fwdStepper->getInitialCondition(),
    fwdStepper,
    dyn_cast<SolverAcceptingStepperBase<Scalar> >(*fwdStepper).getNonconstSolver()
    );
  *out << "\nfwdSensStepper: " << describe(*fwdSensStepper, verbLevel);

  const MEB::InArgs<Scalar> fwdIC = 
    forwardIntegrator->getStepper()->getInitialCondition();

  MEB::InArgs<Scalar> fwdSensIC = 
    createStateAndSensInitialCondition<Scalar>(*fwdSensStepper, fwdIC, B);
  fwdSensStepper->setInitialCondition(fwdSensIC);

  const RCP<IntegratorBase<Scalar> > fwdSensIntegrator =
    forwardIntegrator->cloneIntegrator();
  fwdSensIntegrator->setStepper(fwdSensStepper,
    forwardIntegrator->getFwdTimeRange().upper());
  *out << "\nfwdSensIntegrator: " << describe(*fwdSensIntegrator, verbLevel);

  //
  *out << "\nC) Construct the adjoint sensitivity integrator ...\n";
  //

  RCP<AdjointModelEvaluator<Scalar> > adjModel =
    adjointModelEvaluator<Scalar>(fwdModel, fwdTimeRange);
  adjModel->setFwdStateSolutionBuffer(
    dyn_cast<DefaultIntegrator<Scalar> >(*forwardIntegrator).getTrailingInterpolationBuffer()
    );
  /** \todo 2009/09/16: above: Replace dynamic cast to get trailing interp
   * buffer with the base integator itself.  ToDo this, we will need to get
   * DefaultIntegrator::getPoints(...) to get points off off the
   * trailingInterpBuffer object also.
   */
  *out << "\nadjModel: " << describe(*adjModel, verbLevel);
  
  // lambda(t_final) = 0.0 (for now)
  const RCP<Thyra::VectorBase<Scalar> > lambda_ic = createMember(f_space);
  V_S( lambda_ic.ptr(), 0.0 );
  
  // lambda_dot(t_final,i) = 0.0
  const RCP<Thyra::VectorBase<Scalar> > lambda_dot_ic = createMember(f_space);
  Thyra::V_S( lambda_dot_ic.ptr(), 0.0 );
  
  MEB::InArgs<Scalar> adj_ic = adjModel->getNominalValues();
  adj_ic.set_x(lambda_ic);
  adj_ic.set_x_dot(lambda_dot_ic);
  *out << "\nadj_ic: " << describe(adj_ic, Teuchos::VERB_EXTREME);
  
  const RCP<Thyra::LinearNonlinearSolver<Scalar> > adjTimeStepSolver =
    Thyra::linearNonlinearSolver<Scalar>();
  const RCP<Rythmos::StepperBase<Scalar> > adjStepper =
    forwardIntegrator->getStepper()->cloneStepperAlgorithm();
  *out << "\nadjStepper: " << describe(*adjStepper, verbLevel);
 
  adjStepper->setInitialCondition(adj_ic);
 
  const RCP<IntegratorBase<Scalar> > adjIntegrator = forwardIntegrator->cloneIntegrator();
  adjIntegrator->setStepper(adjStepper, forwardIntegrator->getFwdTimeRange().length());

  //
  *out << "\nD) Solve the forawrd problem storing the state along the way ...\n";
  //

  const double t_final = fwdTimeRange.upper();
  RCP<const Thyra::VectorBase<Scalar> > x_final, x_dot_final;
  get_fwd_x_and_x_dot( *forwardIntegrator, t_final, outArg(x_final), outArg(x_dot_final) );

  *out << "\nt_final = " << t_final << "\n";
  *out << "\nx_final: " << *x_final;
  *out << "\nx_dot_final: " << *x_dot_final;

  //
  *out << "\nE) Solve the forawrd sensitivity equations and compute fwd reduced sens ...\n";
  //

  // Set the initial condition
  TEST_FOR_EXCEPT(true);

  // Solve the fwd sens equations
  TEST_FOR_EXCEPT(true);

  // Compute the reduced gradient at t_f
  TEST_FOR_EXCEPT(true);

  //
  *out << "\nF) Solve the adjoint equations and compute adj reduced sens ...\n";
  //

  // Compute and set the adjoint initial condition
  TEST_FOR_EXCEPT(true);

  // Solve the adjoint equations
  TEST_FOR_EXCEPT(true);

  // Compute the reduced gradient at t_0
  TEST_FOR_EXCEPT(true);

  //
  *out << "\nG) Compare forward and adjoint reduced sens ...\n";
  //

  TEST_FOR_EXCEPT(true);

  //
  *out << "\n*** Leaving BasicDiscreteAdjointStepperTester<Scalar>::testAdjointStepper(...) ...\n";
  //

  return success;

}


// private:


template<class Scalar>
BasicDiscreteAdjointStepperTester<Scalar>::BasicDiscreteAdjointStepperTester()
  : errorTol_(BasicDiscreteAdjointStepperTesterUtils::ErrorTol_default)
{}


} // namespace Rythmos


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define FORWARD_SENSITIVITY_STEPPER_INSTANT(SCALAR) \
  \
  template class BasicDiscreteAdjointStepperTester<SCALAR >; \
  \
  template RCP<BasicDiscreteAdjointStepperTester<SCALAR > > \
  basicDiscreteAdjointStepperTester(); \
  \
  template RCP<BasicDiscreteAdjointStepperTester<SCALAR> > \
  basicDiscreteAdjointStepperTester(const RCP<ParameterList> &paramList);


#endif //Rythmos_BASIC_DISCRETE_ADJOINT_STEPPER_TESTER_DEF_H
