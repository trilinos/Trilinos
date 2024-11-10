// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_STEADYSTATESOLVER_HPP
#define PIRO_STEADYSTATESOLVER_HPP

#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Piro_Helpers.hpp"

namespace Piro {

/** \brief Thyra-based abstract Model Evaluator for steady-states solves
 *  \ingroup Piro_Thyra_solver_grp
 * */
template <typename Scalar>
class SteadyStateSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{
  public:

  /** \name Constructors/initializers */
  //@{
  /** \brief . */
  explicit SteadyStateSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &adjointModel);

  /** \brief . */
  SteadyStateSolver(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &adjointModel,
      int numParameters);
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  //@}

  /** \name Overridden from Thyra::ResponseOnlyModelEvaluatorBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;

  //@}

  /** \name Getters for subclasses. */
  //@{

  /** \brief . */
  const Thyra::ModelEvaluator<Scalar> &getModel() const;

  /** \brief . */
  int num_p() const;
  /** \brief . */
  int num_g() const;

  /** \brief . */
  SENS_METHOD getSensitivityMethod();
  //@}

  /** \name Setters for subbclasses */
  /** \brief . */
  void setSensitivityMethod(const std::string& sensitivity_method_string);
  //@}

  //@}

  protected:
  /** \name Service methods for subclasses. */
  //@{

  /** \brief . */
  void evalConvergedModelResponsesAndSensitivities(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs,
      Teuchos::ParameterList& appParams) const;

  /** \brief . */
  void evalReducedHessian(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs,
      Teuchos::ParameterList& appParams) const;
  //@}

  private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDp_op_impl(int j, int l) const;
  //@}

  /** \name Internal implemention methods. */
  //@{
  /** \brief Implementation of createInArgs . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgsImpl() const;
  //@}

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_, adjointModel_;

  int num_p_;
  int num_g_;

  SENS_METHOD sensitivityMethod_;
};

}

#include "Piro_SteadyStateSolver_Def.hpp"
#endif /*PIRO_STEADYSTATESOLVER_HPP*/
