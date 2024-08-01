// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_VELOCITYVERLETSOLVER_H
#define PIRO_VELOCITYVERLETSOLVER_H

#include <iostream>

#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Piro_ObserverBase.hpp"

#include "Thyra_AdaptiveSolutionManager.hpp"

namespace Piro {

template <typename Scalar>
class VelocityVerletSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar> {


  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  VelocityVerletSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams,
                       const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
                       const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
                       const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer = Teuchos::null
                      );

  //@}

/** \name Overridden from Thyra::ModelEvaluatorBase. */
  //@{
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_dot() const;
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

  /** \brief . */
  void setProblemParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  void setSolverParamDefaults(Teuchos::ParameterList* appParams_);
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidVelocityVerletParameters() const;

   //These are set in the constructor and used in evalModel
   mutable Teuchos::RCP<Teuchos::ParameterList> appParams;
   Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model;
   Teuchos::RCP<Piro::ObserverBase<Scalar> > observer;
   Teuchos::RCP<Thyra::AdaptiveSolutionManager> solMgr;
   Teuchos::RCP<Teuchos::FancyOStream> out;
   Teuchos::EVerbosityLevel solnVerbLevel;

   int num_p;
   int num_g;

   int numTimeSteps;
   Scalar t_init, t_final, delta_t;

};

}

#include "Piro_VelocityVerletSolver_Def.hpp"

#endif
