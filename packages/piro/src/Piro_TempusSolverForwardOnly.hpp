// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_TEMPUSSOLVERFORWARDONLY_H
#define PIRO_TEMPUSSOLVERFORWARDONLY_H

#include "Piro_ConfigDefs.hpp"
#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_IntegratorObserver.hpp"

#include "Piro_ObserverBase.hpp"

#include <map>
#include <string>

namespace Piro {

/** \brief Thyra-based Model Evaluator for Tempus solves.
 */
template <typename Scalar>
class TempusSolverForwardOnly
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{
public:
  /** \name Constructors/initializers */
  //@{

  /** \brief Initializes the internals, though the object is a blank slate. To initialize it call <code>initialize</code> */
  TempusSolverForwardOnly();

  /** \brief Initialize with internally built objects according to the given parameter list. */
  TempusSolverForwardOnly(
      const Teuchos::RCP<Teuchos::ParameterList> &appParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
      const Teuchos::RCP<Tempus::IntegratorObserver<Scalar> > &observer = Teuchos::null);

  /** \brief Initialize using prebuilt objects. */
  TempusSolverForwardOnly(
      const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
      const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
      Scalar finalTime,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &initialConditionModel = Teuchos::null,
      Teuchos::EVerbosityLevel verbosityLevel = Teuchos::VERB_DEFAULT);
  //@}

  /** \brief Initialize using prebuilt objects - supplying initial time value. */
  TempusSolverForwardOnly(
      const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
      const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
      Scalar initialTime,
      Scalar finalTime,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &initialConditionModel = Teuchos::null,
      Teuchos::EVerbosityLevel verbosityLevel = Teuchos::VERB_DEFAULT);
  //@}

  void initialize(
      const Teuchos::RCP<Teuchos::ParameterList> &appParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
      const Teuchos::RCP<Tempus::IntegratorObserver<Scalar> > &observer = Teuchos::null);

  Teuchos::RCP<const Tempus::IntegratorBasic<Scalar> > getTempusIntegrator() const;

  Teuchos::RCP<const Thyra::NonlinearSolverBase<Scalar> > getTimeStepSolver() const;

  /** \name Overridden from Thyra::ModelEvaluatorBase. */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  //@}


private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase. */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;
  //@}

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDp_op_impl(int j, int l) const;

  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidTempusParameters() const;
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > getUnderlyingModel() const;
  Scalar get_t_initial() const;
  Scalar get_t_final() const;
  int get_num_p() const;
  int get_num_g() const;

  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > fwdStateIntegrator;

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > initialConditionModel;

  Teuchos::RCP<Teuchos::FancyOStream> out;
  Teuchos::EVerbosityLevel solnVerbLevel;

  bool isInitialized;
};

/** \brief Non-member constructor function */
template <typename Scalar>
Teuchos::RCP<TempusSolverForwardOnly<Scalar> >
tempusSolverForwardOnly(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<ObserverBase<Scalar> > &piroObserver);

}

/** \class Piro::TempusSolverForwardOnly
 *  \ingroup Piro_Thyra_solver_grp
 * */

#include "Piro_TempusSolverForwardOnly_Def.hpp"

#endif
