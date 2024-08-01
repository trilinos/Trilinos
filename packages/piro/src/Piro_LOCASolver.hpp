// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_LOCASOLVER_HPP
#define PIRO_LOCASOLVER_HPP

#include "Piro_SteadyStateSolver.hpp"

#include "Piro_ObserverBase.hpp"

#include "LOCA.H"
#include "LOCA_Thyra.H"
#include "LOCA_Thyra_SaveDataStrategy.H"

#include "Teuchos_ParameterList.hpp"

namespace Piro {

/** \brief Thyra-based Model Evaluator for LOCA solves
 *  \ingroup Piro_Thyra_solver_grp
 * */
template <typename Scalar>
class LOCASolver : public SteadyStateSolver<Scalar> {
public:
  /** \name Constructor/Destructor */
  //@{
  /** \brief Constructs a LOCASolver instance given a model and optionally a data saving strategy . */
  LOCASolver(
      const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel = Teuchos::null,
      const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy = Teuchos::null);

  ~LOCASolver();
  //@}

  //! Return the current nonlinear solver pointer.
  Teuchos::RCP<NOX::Solver::Generic>
  getSolver();

  //! Return stepper parameters
  Teuchos::ParameterList &
  getStepperParams();

  //! Return step size parameters
  Teuchos::ParameterList &
  getStepSizeParams();

  //! Return the underlying LOCA stepper
  Teuchos::RCP<LOCA::Stepper>
  getStepper();

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;
  //@}

  Teuchos::RCP<Teuchos::ParameterList> piroParams_;
  Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> saveDataStrategy_;

  Teuchos::RCP<LOCA::GlobalData> globalData_;
  mutable LOCA::ParameterVector paramVector_;
  Teuchos::RCP<LOCA::Thyra::Group> group_;
  Teuchos::RCP<LOCA::StatusTest::Abstract> locaStatusTests_;
  Teuchos::RCP<NOX::StatusTest::Generic> noxStatusTests_;
  NOX::Utils utils_;

  Teuchos::RCP<LOCA::Stepper> stepper_;
  mutable bool first_;
  
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model_; 
};


template <typename Scalar>
Teuchos::RCP<LOCASolver<Scalar> >
observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel = Teuchos::null,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer = Teuchos::null);

} // namespace Piro

#endif /* PIRO_LOCASOLVER_HPP */
