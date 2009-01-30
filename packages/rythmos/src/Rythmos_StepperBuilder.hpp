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

#ifndef Rythmos_STEPPER_BUILDER_H
#define Rythmos_STEPPER_BUILDER_H

#include "Rythmos_Types.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitRKStepper.hpp"
#include "Rythmos_ExplicitTaylorPolynomialStepper.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace Rythmos {

inline const std::string BackwardEuler_name() { return  "Backward Euler"; }
inline const std::string BackwardEulerSettings_name() { return  "Backward Euler Settings"; }

inline const std::string ImplicitBDF_name() { return  "Implictit BDF"; }
inline const std::string ImplicitBDFSettings_name() { return  "Implicit BDF Settings"; }

inline const std::string ForwardEuler_name() { return  "Forward Euler"; }
inline const std::string ForwardEulerSettings_name() { return  "Forward Euler Settings"; }

inline const std::string ExplicitRK_name() { return  "Explicit RK"; }
inline const std::string ExplicitRKSettings_name() { return  "Explicit RK Settings"; }

inline const std::string ImplicitRK_name() { return  "Implicit RK"; }
inline const std::string ImplicitRKSettings_name() { return  "Implicit RK Settings"; }

inline const std::string ExplicitTP_name() { return  "Explicit Taylor Polynomial"; }
inline const std::string ExplicitTPSettings_name() { return  "Explicit Taylor Polynomial Settings"; }

} // namespace Rythmos

namespace {

inline const std::string StepperType_name() { return  "Stepper Type"; }
inline const std::string StepperType_default() { return  Rythmos::BackwardEuler_name(); }

enum E_StepperBuilderSelectionTypes {
  RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_BACKWARDEULER,
  RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_IMPLICITBDF,
  RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_FORWARDEULER,
  RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_EXPLICITRK,
  RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_IMPLICITRK,
  RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_EXPLICITTP
};

Teuchos::Array<std::string>
  S_StepperBuilderSelectionTypes = Teuchos::tuple<std::string>(
      Rythmos::BackwardEuler_name(),
      Rythmos::ImplicitBDF_name(),
      Rythmos::ForwardEuler_name(),
      Rythmos::ExplicitRK_name(),
      Rythmos::ImplicitRK_name(),
      Rythmos::ExplicitTP_name()
      );

  Teuchos::Array<std::string>
  S_StepperBuilderSettingsTypes = Teuchos::tuple<std::string>(
      Rythmos::BackwardEulerSettings_name(),
      Rythmos::ImplicitBDFSettings_name(),
      Rythmos::ForwardEulerSettings_name(),
      Rythmos::ExplicitRKSettings_name(),
      Rythmos::ImplicitRKSettings_name(),
      Rythmos::ExplicitTPSettings_name()
      );

const RCP<Teuchos::StringToIntegralParameterEntryValidator<E_StepperBuilderSelectionTypes> >
  stepperBuilderSelectionTypeValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<E_StepperBuilderSelectionTypes>(
        S_StepperBuilderSelectionTypes,
        Teuchos::tuple<E_StepperBuilderSelectionTypes>(
          RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_BACKWARDEULER,
          RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_IMPLICITBDF,
          RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_FORWARDEULER,
          RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_EXPLICITRK,
          RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_IMPLICITRK,
          RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_EXPLICITTP
          ),
        StepperType_name()
        )
      );

} // namespace

namespace Rythmos {

template<class Scalar>
  class StepperBuilder : virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief . */
  StepperBuilder() {}

  /** \brief . */
  RCP<StepperBase<Scalar> > create();
  
  /** \brief . */
  void setNonlinearSolver(const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver);

  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const RCP<Teuchos::ParameterList> & paramList);
  
  /** \brief . */
  RCP<Teuchos::ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<Teuchos::ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const Teuchos::ParameterList> getValidParameters() const;
 
  //@}
private:
  RCP<ParameterList> paramList_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > nlSolver_;

};

// nonmember constructor
template<class Scalar>
RCP<StepperBuilder<Scalar> > stepperBuilder() {
  return Teuchos::rcp(new StepperBuilder<Scalar>() );
}


template<class Scalar>
void StepperBuilder<Scalar>::setParameterList(const RCP<Teuchos::ParameterList> & paramList)
{
  if (!is_null(paramList)) {
    paramList->validateParameters(*this->getValidParameters());
    paramList_ = paramList;
  }
}

template<class Scalar>
RCP<Teuchos::ParameterList> StepperBuilder<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<const Teuchos::ParameterList> StepperBuilder<Scalar>::getValidParameters() const
{
  static RCP<ParameterList> validPL;
  if (is_null(validPL)) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set( StepperType_name(), StepperType_default(),
        "", stepperBuilderSelectionTypeValidator 
        );
    // Set up Backward Euler settings sublist
    ParameterList& beSettings = pl->sublist(BackwardEulerSettings_name());
    {
      Rythmos::BackwardEulerStepper<double> stepper;
      beSettings.setParameters(
        *stepper.getValidParameters()).disableRecursiveValidation();
    }
    // Set up Implicit BDF settings sublist
    ParameterList& ibdfSettings = pl->sublist(ImplicitBDFSettings_name());
    {
      Rythmos::ImplicitBDFStepper<double> stepper;
      ibdfSettings.setParameters(
        *stepper.getValidParameters()).disableRecursiveValidation();
    }
    // Set up Forward Euler settings sublist
    ParameterList& feSettings = pl->sublist(ForwardEulerSettings_name());
    {
      Rythmos::ForwardEulerStepper<double> stepper;
      feSettings.setParameters(
        *stepper.getValidParameters()).disableRecursiveValidation();
    }
    // Set up Explicit RK settings sublist
    ParameterList& erkSettings = pl->sublist(ExplicitRKSettings_name());
    {
      Rythmos::ExplicitRKStepper<double> stepper;
      erkSettings.setParameters(
        *stepper.getValidParameters()).disableRecursiveValidation();
    }
    // Set up Implicit RK settings sublist
    ParameterList& irkSettings = pl->sublist(ImplicitRKSettings_name());
    {
      Rythmos::ImplicitRKStepper<double> stepper;
      irkSettings.setParameters(
        *stepper.getValidParameters()).disableRecursiveValidation();
    }
    // Set up Explicit TP settings sublist
    ParameterList& etpSettings = pl->sublist(ExplicitTPSettings_name());
    {
      Rythmos::ExplicitTaylorPolynomialStepper<double> stepper;
      etpSettings.setParameters(
        *stepper.getValidParameters()).disableRecursiveValidation();
    }

    validPL = pl;
  }

  return validPL;
}

template<class Scalar>
RCP<Teuchos::ParameterList> StepperBuilder<Scalar>::unsetParameterList()
{
  RCP<ParameterList> oldPL = paramList_;
  paramList_ = Teuchos::null;
  return(oldPL);
}

template<class Scalar>
RCP<StepperBase<Scalar> > StepperBuilder<Scalar>::create()
{
  TEST_FOR_EXCEPTION( is_null(paramList_), std::logic_error, "Error!, you must first set a parameter list on the builder!" );
  RCP<StepperBase<Scalar> > stepper;
  E_StepperBuilderSelectionTypes stepperTypeEnum = stepperBuilderSelectionTypeValidator->getIntegralValue(
      *paramList_, StepperType_name(), StepperType_default()
      );
  if (stepperTypeEnum == RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_BACKWARDEULER) {
    stepper = backwardEulerStepper<Scalar>();
  } else if (stepperTypeEnum == RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_IMPLICITBDF) {
    stepper = implicitBDFStepper<Scalar>();
  } else if (stepperTypeEnum == RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_FORWARDEULER) {
    stepper = forwardEulerStepper<Scalar>();
  } else if (stepperTypeEnum == RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_EXPLICITRK) {
    stepper = explicitRKStepper<Scalar>();
  } else if (stepperTypeEnum == RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_IMPLICITRK) {
    stepper = implicitRKStepper<Scalar>();
  } else if (stepperTypeEnum == RYTHMOS_STEPPERBUILDER_SELECTION_TYPE_EXPLICITTP) {
    stepper = explicitTaylorPolynomialStepper<Scalar>();
  } else {
    TEST_FOR_EXCEPTION( true, std::logic_error,
        "Error!  Code should never get here!"
        );
  }
  // Set the appropriate parameterlist on the stepper
  RCP<ParameterList> settingsPL = Teuchos::sublist(paramList_,S_StepperBuilderSettingsTypes[stepperTypeEnum]);
  stepper->setParameterList(settingsPL);
  // Set Nonlinear Solver on the Stepper if it is of type Rythmos::SolverAcceptingStepperBase
  {
    RCP<SolverAcceptingStepperBase<Scalar> > saSB = Teuchos::rcp_dynamic_cast<SolverAcceptingStepperBase<Scalar> >(stepper,false);
    if (!is_null(saSB)) {
      TEST_FOR_EXCEPTION( is_null(nlSolver_), std::logic_error,
          "Error!  Attempting to build a SolverAcceptingStepperBase object without a valid Thyra::NonlinearSolverBase."
          );
      saSB->setSolver(nlSolver_);
    }
  }
  return(stepper);
}

template<class Scalar>
void StepperBuilder<Scalar>::setNonlinearSolver(const RCP<Thyra::NonlinearSolverBase<Scalar> >& nlSolver)
{
  TEST_FOR_EXCEPTION(is_null(nlSolver), std::logic_error,
      "Error!  Attempting to call setNonlinearSolver with null RCP"
      );
  nlSolver_ = nlSolver;
}

} // namespace Rythmos

#endif //Rythmos_STEPPER_BUILDER_H
