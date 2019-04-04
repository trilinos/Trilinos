// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PIRO_TEMPUSSOLVER_H
#define PIRO_TEMPUSSOLVER_H

#include "Piro_ConfigDefs.hpp"
#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorObserverBasic.hpp"

#include "Piro_ObserverBase.hpp"

#include "Piro_TempusStepperFactory.hpp"
#include "Piro_TempusStepControlFactory.hpp"

// This "define" turns on the extended template interface in TempusSolver.
// Is it necessary??
#if defined(HAVE_PIRO_TEMPUS) 
#define ALBANY_BUILD
#endif

#ifdef ALBANY_BUILD
#include "Kokkos_DefaultNode.hpp"
#endif

#include <map>
#include <string>

namespace Piro {

/** \brief Thyra-based Model Evaluator for Tempus solves
 *  \ingroup Piro_Thyra_solver_grp
 * */
#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal = int, typename GlobalOrdinal = LocalOrdinal,
          typename Node = KokkosClassic::DefaultNode::DefaultNodeType>
#else
template <typename Scalar>
#endif
class TempusSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{
public:
  /** \name Constructors/initializers */
  //@{

  /** \brief Initializes the internals, though the object is a blank slate. To initialize it call <code>initialize</code> */
  TempusSolver();

  /** \brief Initialize with internally built objects according to the given parameter list. */
  TempusSolver(
      const Teuchos::RCP<Teuchos::ParameterList> &appParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      bool computeSensitivities,
      const Teuchos::RCP<Piro::ObserverBase<Scalar> > &piroObserver = Teuchos::null);

  /** \brief Initialize using prebuilt objects. */
  TempusSolver(
      const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
      const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      Scalar finalTime,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &initialConditionModel = Teuchos::null,
      Teuchos::EVerbosityLevel verbosityLevel = Teuchos::VERB_DEFAULT);
  //@}

  /** \brief Initialize using prebuilt objects - supplying initial time value. */
  TempusSolver(
      const Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > &stateIntegrator,
      const Teuchos::RCP<Tempus::Stepper<Scalar> > &stateStepper,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      Scalar initialTime,
      Scalar finalTime,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &initialConditionModel = Teuchos::null,
      Teuchos::EVerbosityLevel verbosityLevel = Teuchos::VERB_DEFAULT);
  //@}

  void initialize(
      const Teuchos::RCP<Teuchos::ParameterList> &appParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model);

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

  void addStepperFactory(const std::string & stepperName,
                         const Teuchos::RCP<Piro::TempusStepperFactory<Scalar> > & stepperFactories);

  void addStepControlFactory(const std::string & stepControlName,
                             const Teuchos::RCP<Piro::TempusStepControlFactory<Scalar> > & step_control_strategy);

  //! Set start time for time-integration
  void
  setStartTime(const Scalar start_time);

  //! Get start time for time-integration
  Scalar
  getStartTime() const;

  //! Set final time for time-integration
  void
  setFinalTime(const Scalar final_time);

  //! Get final time for time-integration
  Scalar
  getFinalTime() const;
  
  //! Set initial time step for time-integration
  void
  setInitTimeStep(const Scalar init_time_step);
  
  //! Get initial time step for time-integration
  Scalar
  getInitTimeStep() const;

  //! Set initial time, initial solution, velocity and acceleration
  void setInitialState(Scalar t0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > x0,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdot0 = Teuchos::null,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xdotdot0 = Teuchos::null);

  //! Set initial guess for Newton method
  void setInitialGuess( Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess = Teuchos::null); 

  //! Return RCP to Tempus::SolutionHistory
  Teuchos::RCP<Tempus::SolutionHistory<Scalar> > 
  getSolutionHistory() const; 
 
  //! Return Thyra nonlinear solver underlying Tempus::Stepper object  
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > 
  getSolver() const;

  Tempus::Status 
  getTempusIntegratorStatus() const;

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

  Teuchos::RCP<Tempus::IntegratorBasic<Scalar> > fwdStateIntegrator;
  Teuchos::RCP<Tempus::Stepper<Scalar> > fwdStateStepper;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > fwdTimeStepSolver;

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model;
  Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > thyraModel;
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > initialConditionModel;

  Scalar t_initial;
  Scalar t_final;

  int num_p;
  int num_g;

  bool computeSensitivities_;

  Teuchos::RCP<Teuchos::FancyOStream> out;
  Teuchos::EVerbosityLevel solnVerbLevel;

  // used for adding user defined steppers externally, this gives us "the open-close principal"
  std::map<std::string,Teuchos::RCP<Piro::TempusStepperFactory<Scalar> > > stepperFactories;

  std::map<std::string,Teuchos::RCP<Piro::TempusStepControlFactory<Scalar> > > stepControlFactories;

  bool isInitialized;

  Teuchos::RCP<Piro::ObserverBase<Scalar> > piroObserver_;

  bool supports_x_dotdot_; 
 
  //! Set observer
  void setObserver(); 

  //! Boolean to tell TempusSolver whether or not to abort if a transient solve fails 
  bool abort_on_failure_; 

};

/** \brief Non-member constructor function */
#ifdef ALBANY_BUILD
template <typename Scalar, typename LocalOrdinal = int, typename GlobalOrdinal = LocalOrdinal,
          typename Node = KokkosClassic::DefaultNode::DefaultNodeType>
Teuchos::RCP<TempusSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#else
template <typename Scalar>
Teuchos::RCP<TempusSolver<Scalar> >
#endif
tempusSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<ObserverBase<Scalar> > &piroObserver);

}

/** \class Piro::TempusSolver
 *  \ingroup Piro_Thyra_solver_grp
 * */

#include "Piro_TempusSolver_Def.hpp"

#endif
