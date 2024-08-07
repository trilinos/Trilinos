// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_TRAPEZOIDRULESOLVER_H
#define PIRO_TRAPEZOIDRULESOLVER_H

#include <iostream>

#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Piro_ObserverBase.hpp"

#include "Piro_NOXSolver.hpp"
#include "Thyra_AdaptiveSolutionManager.hpp"

namespace Piro {

template <typename Scalar>
class TrapezoidDecorator
    : public Thyra::ModelEvaluatorDelegatorBase<Scalar> {

  public:

  /** \name Constructors/initializers */
  //@{

  TrapezoidDecorator(
                const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &modelEvaluator
                );

  //@}

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_x() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_x_dot() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > get_x_dotdot() const;

   void reportFinalPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint, const bool wasSolved);



  //! Method to give info to compute xDotDot(x), so that the
  // NOX solver can treat the time dep problem as steady
  void injectData(const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_,
                  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_pred_a_, Scalar fdt2_,
                  const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_pred_v_, Scalar tdt_,
                  Scalar time_ );

  // Resize internal arrays when mesh is adapted
  void resize(const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x_);

  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar> > DMEWSF;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_pred_a;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_pred_v;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_save;
   Scalar fdt2;
   Scalar tdt;
   Scalar time;

   Teuchos::RCP<Teuchos::FancyOStream> out;

};

template <typename Scalar>
class TrapezoidRuleSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar> {


  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  TrapezoidRuleSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams,
                      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
                      const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
                      const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer = Teuchos::null
                      );

  //@}

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
  /** \brief . */
  Teuchos::RCP<Piro::NOXSolver<Scalar> > getNOXSolver() const;
  /** \brief . */
  Teuchos::RCP<Piro::TrapezoidDecorator<Scalar> > getDecorator() const;
  /** \brief . */
  Teuchos::RCP<Thyra::AdaptiveSolutionManager> getSolutionManager() const;
  /** \brief .*/
  void disableCalcInitAccel() { calc_init_accel_ = false; }; 
  /** \brief .*/
  void enableCalcInitAccel() { calc_init_accel_ = true; }; 
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

  /** \brief . */
  void setProblemParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  void setSolverParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidTrapezoidRuleParameters() const;

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP<Piro::TrapezoidDecorator<Scalar> > model;
   Teuchos::RCP<Piro::ObserverBase<Scalar> > observer;
   Teuchos::RCP<Thyra::AdaptiveSolutionManager> solMgr;
   Teuchos::RCP<Teuchos::FancyOStream> out;
   Teuchos::EVerbosityLevel solnVerbLevel;

   Teuchos::RCP<Piro::NOXSolver<Scalar> > noxSolver;

   int num_p;
   int num_g;

   int numTimeSteps;
   Scalar t_init, t_final, delta_t;

   mutable bool calc_init_accel_{true}; 
};

}

#include "Piro_TrapezoidRuleSolver_Def.hpp"

#endif
