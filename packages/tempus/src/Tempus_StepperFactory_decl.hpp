// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperFactory_decl_hpp
#define Tempus_StepperFactory_decl_hpp

#include "Teuchos_ParameterList.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
#include "Tempus_StepperExplicitRK.hpp"
#include "Tempus_StepperBDF2.hpp"
#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitDForm.hpp"
#include "Tempus_StepperNewmarkExplicitAForm.hpp"
#include "Tempus_StepperHHTAlpha.hpp"
#include "Tempus_StepperRKButcherTableau.hpp"
#include "Tempus_StepperIMEX_RK.hpp"
#include "Tempus_StepperIMEX_RK_Partition.hpp"
#include "Tempus_StepperLeapfrog.hpp"
#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_StepperTrapezoidal.hpp"
#include "Tempus_StepperSubcycling.hpp"

#include "NOX_Thyra.H"

namespace Tempus {

/** \brief Stepper factory.
 *
 */
template<class Scalar>
class StepperFactory
{
public:

  /// Constructor
  StepperFactory(){}

  /// Destructor
  virtual ~StepperFactory() {}


/// \name Stepper creation methods
//@{
  /// Create stepper from stepper type.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    std::string stepperType = "Forward Euler",
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
      model = Teuchos::null);

  /// Create stepper from a ParameterList.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
      model = Teuchos::null);

  /// Create multi-stepper from ParameterList.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models);
//@}


/// \name Create individual Stepppers from ParameterList
//@{
  /// Create the Operator-Split Stepper from a ParameterList.
  Teuchos::RCP<StepperOperatorSplit<Scalar> >
  createStepperOperatorSplit(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSubcycling<Scalar> >
  createStepperSubcycling(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperIMEX_RK_Partition<Scalar> >
  createStepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperIMEX_RK<Scalar> >
  createStepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperHHTAlpha<Scalar> >
  createStepperHHTAlpha(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperNewmarkImplicitDForm<Scalar> >
  createStepperNewmarkImplicitDForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> >
  createStepperNewmarkImplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperBDF2<Scalar> >
  createStepperBDF2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperTrapezoidal<Scalar> >
  createStepperTrapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperBackwardEuler<Scalar> >
  createStepperBackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperNewmarkExplicitAForm<Scalar> >
  createStepperNewmarkExplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperLeapfrog<Scalar> >
  createStepperLeapfrog(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperForwardEuler<Scalar> >
  createStepperForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  // RK Explicit Methods

  Teuchos::RCP<StepperERK_General<Scalar> >
  createStepperERK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_ForwardEuler<Scalar> >
  createStepperERK_ForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_4Stage4thOrder<Scalar> >
  createStepperERK_4Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_3_8Rule<Scalar> >
  createStepperERK_3_8Rule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_4Stage3rdOrderRunge<Scalar> >
  createStepperERK_4Stage3rdOrderRunge(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_5Stage3rdOrderKandG<Scalar> >
  createStepperERK_5Stage3rdOrderKandG(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_3Stage3rdOrder<Scalar> >
  createStepperERK_3Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_3Stage3rdOrderTVD<Scalar> >
  createStepperERK_3Stage3rdOrderTVD(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    const std::string stepperType = "RK Explicit 3 Stage 3rd order TVD");

  Teuchos::RCP<StepperERK_SSPERK54<Scalar> >
  createStepperERK_SSPERK54(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_3Stage3rdOrderHeun<Scalar> >
  createStepperERK_3Stage3rdOrderHeun(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_Midpoint<Scalar> >
  createStepperERK_Midpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_Trapezoidal<Scalar> >
  createStepperERK_Trapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::string stepperType = "RK Explicit Trapezoidal");

  Teuchos::RCP<StepperERK_Ralston<Scalar> >
  createStepperERK_Ralston(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::string stepperType = "RK Explicit Ralston");

  Teuchos::RCP<StepperERK_BogackiShampine32<Scalar> >
  createStepperERK_BogackiShampine32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperERK_Merson45<Scalar> >
  createStepperERK_Merson45(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  // RK Implicit Methods

  Teuchos::RCP<StepperDIRK_General<Scalar> >
  createStepperDIRK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperDIRK_BackwardEuler<Scalar> >
  createStepperDIRK_BackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_2Stage2ndOrder<Scalar> >
  createStepperSDIRK_2Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_3Stage2ndOrder<Scalar> >
  createStepperSDIRK_3Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_SSPDIRK22<Scalar> >
  createStepperSDIRK_SSPDIRK22(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_SSPDIRK32<Scalar> >
  createStepperSDIRK_SSPDIRK32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_SSPDIRK23<Scalar> >
  createStepperSDIRK_SSPDIRK23(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_SSPDIRK33<Scalar> >
  createStepperSDIRK_SSPDIRK33(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_2Stage3rdOrder<Scalar> >
  createStepperSDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperEDIRK_2Stage3rdOrder<Scalar> >
  createStepperEDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperDIRK_1StageTheta<Scalar> >
  createStepperDIRK_1StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperEDIRK_2StageTheta<Scalar> >
  createStepperEDIRK_2StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperEDIRK_TrapezoidalRule<Scalar> >
  createStepperEDIRK_TrapezoidalRule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_ImplicitMidpoint<Scalar> >
  createStepperSDIRK_ImplicitMidpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperDIRK_1Stage1stOrderRadauIA<Scalar> >
  createStepperDIRK_1Stage1stOrderRadauIA(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperDIRK_2Stage2ndOrderLobattoIIIB<Scalar> >
  createStepperDIRK_2Stage2ndOrderLobattoIIIB(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_5Stage4thOrder<Scalar> >
  createStepperSDIRK_5Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_3Stage4thOrder<Scalar> >
  createStepperSDIRK_3Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_5Stage5thOrder<Scalar> >
  createStepperSDIRK_5Stage5thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  Teuchos::RCP<StepperSDIRK_21Pair<Scalar> >
  createStepperSDIRK_21Pair(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  void createSubSteppers(
    Teuchos::RCP<StepperOperatorSplit<Scalar> > stepper,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);
//@}


/// \name Methods to set member data from ParameterList
//@{
  /// Set Stepper member data from the ParameterList.
  void setStepperValues(
    Teuchos::RCP<Stepper<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  /// Create a tableau from the ParameterList.
  Teuchos::RCP<RKButcherTableau<Scalar> > createTableau(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  /// Set StepperExplicit member data from the ParameterList.
  void setStepperExplicitValues(
    Teuchos::RCP<StepperExplicit<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  /// Set StepperImplicit member data from the ParameterList.
  void setStepperImplicitValues(
    Teuchos::RCP<StepperImplicit<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  /// Set StepperRK member data from the ParameterList.
  void setStepperRKValues(
    Teuchos::RCP<StepperExplicitRK<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  /// Set solver from ParameterList.
  void setStepperSolverValues(
    Teuchos::RCP<StepperImplicit<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  /// Set StepperDIRK member data from the ParameterList.
  void setStepperDIRKValues(
    Teuchos::RCP<StepperDIRK<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL);

  void setTableausPartition(
    Teuchos::RCP<StepperIMEX_RK_Partition<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::string stepperType);

  void setTableaus(Teuchos::RCP<StepperIMEX_RK<Scalar> > stepper,
                   Teuchos::RCP<Teuchos::ParameterList> stepperPL,
                   std::string stepperType);
//@}

private:

  /// Stepper Factory.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);

};


} // namespace Tempus
#endif // Tempus_StepperFactory_decl_hpp
