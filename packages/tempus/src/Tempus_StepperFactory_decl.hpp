//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperFactory_decl_hpp
#define Tempus_StepperFactory_decl_hpp

#include "Teuchos_ParameterList.hpp"

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"

namespace Tempus {

/** \brief Stepper factory.
 *
 */
template <class Scalar>
class StepperFactory {
 public:
  /// Constructor
  StepperFactory() {}

  /// Destructor
  virtual ~StepperFactory() {}

  /// \name Stepper constructors
  //@{
  /// Create stepper from stepper type.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
      std::string stepperType = "Forward Euler",
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model =
          Teuchos::null);

  /// Create stepper from a ParameterList.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
      Teuchos::RCP<Teuchos::ParameterList> stepperPL,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model =
          Teuchos::null);

  /// Create multi-stepper from ParameterList.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
      Teuchos::RCP<Teuchos::ParameterList> stepperPL,
      std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models);
  //@}

 private:
  /// Stepper Factory.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
      std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> stepperPL,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);
};

}  // namespace Tempus
#endif  // Tempus_StepperFactory_decl_hpp
