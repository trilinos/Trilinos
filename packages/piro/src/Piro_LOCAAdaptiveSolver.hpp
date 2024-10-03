// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_LOCAADAPTIVESOLVER_HPP
#define PIRO_LOCAADAPTIVESOLVER_HPP

#include "Piro_SteadyStateSolver.hpp"

#include "Piro_ObserverBase.hpp"

#include "LOCA.H"
#include "LOCA_Thyra.H"
#include "LOCA_Thyra_SaveDataStrategy.H"
#include "Thyra_AdaptiveSolutionManager.hpp"

#include "LOCA_AdaptiveStepper.H"

#include "Teuchos_ParameterList.hpp"

namespace Piro {

/** \brief Thyra-based Model Evaluator for LOCAAdaptive solves
 *  \ingroup Piro_Thyra_solver_grp
 * */
template <typename Scalar>
class LOCAAdaptiveSolver : public SteadyStateSolver<Scalar> {
public:
  /** \name Constructor/Destructor */
  //@{
  /** \brief Constructs a LOCAAdaptiveSolver instance given a model and optionally a data saving strategy . */
  LOCAAdaptiveSolver(
      const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &ajdointModel,
      const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
      const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy);

  ~LOCAAdaptiveSolver();
  //@}

  //! Update the final solution to the main solver
  void reportFinalPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint, const bool /* wasSolved */)
       { finalPoint_ = Teuchos::rcpFromRef(finalPoint); }

  //! Returns the underlying stepper
  Teuchos::RCP<LOCA::AdaptiveStepper>
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
  const Teuchos::RCP<Thyra::AdaptiveSolutionManager> solMgr_;
  Teuchos::RCP<LOCA::StatusTest::Abstract> locaStatusTests_;
  Teuchos::RCP<NOX::StatusTest::Generic> noxStatusTests_;
  NOX::Utils utils_;

  Teuchos::RCP<LOCA::AdaptiveStepper> stepper_;

  //! Holds the final solution
  Teuchos::RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> > finalPoint_;
  
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model_; 
};


template <typename Scalar>
Teuchos::RCP<LOCAAdaptiveSolver<Scalar> >
observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel,
    const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer);

} // namespace Piro

#endif /* PIRO_LOCAADAPTIVESOLVER_HPP */
