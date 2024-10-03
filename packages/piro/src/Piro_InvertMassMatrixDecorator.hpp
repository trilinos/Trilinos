// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_INVERTMASSMATRIXDECORATOR_H
#define PIRO_INVERTMASSMATRIXDECORATOR_H

#include <iostream>

#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveHelpers.hpp"

#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Piro_config.hpp"

namespace Piro {

template<typename Scalar>
class InvertMassMatrixDecorator
    : public Thyra::ModelEvaluatorDefaultBase<Scalar>
{

  public:

  /** \name Constructors/initializers */
  //@{

  /** \brief Takes the number of elements in the discretization . */
  InvertMassMatrixDecorator(
                Teuchos::RCP<Teuchos::ParameterList> stratParams,
                Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
                bool massMatrixIsConstant=true,
                bool lumpMassMatrix=false,
                bool massMatrixIsCoeffOfSecondDeriv=false
                );

  //@}

  ~InvertMassMatrixDecorator();


  /** \name Overridden from Thyra::ModelEvaluator . */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_multiplier_space(int j) const;

  Teuchos::ArrayView<const std::string> get_g_names(int j) const;

  /** \brief . */
  //Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP< Thyra::LinearOpBase< Scalar > > create_W_op () const;
  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;

  Teuchos::RCP<Thyra::PreconditionerBase<Scalar> > create_W_prec() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgsImpl() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */

  void reportFinalPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& finalPoint, const bool wasSolved);

  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

  private:
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_multiplier_space() const;
  /** \brief . */
  //Teuchos::RCP<const Epetra_Vector> get_x_init() const;
  /** \brief . */
  //Teuchos::RCP<const Epetra_Vector> get_x_dot_init() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;

  Teuchos::RCP<const Teuchos::Array<std::string> >
    get_p_names(int l) const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;

   //! Operator form of df/dp for distributed parameters - for transient adjoint sensitivities
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> create_DfDp_op_impl(int j) const
  {
    return model->create_DfDp_op(j); 
  }


  //@}

  private:

   //These are set in the constructor and used in evalModel
   Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model;
   Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot;

   Teuchos::RCP<Thyra::LinearOpBase<Scalar> > massMatrix;
   Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory;

   bool massMatrixIsConstant; // User Setting
   bool lumpMassMatrix; // User Setting to rowSum Matrix
   bool massMatrixIsCoeffOfSecondDeriv; // Set to true for x_dotdot acceleration problems
   Teuchos::RCP<Thyra::VectorBase<Scalar> > invDiag;

   // The following get modified in evalModel and so are mutable
   mutable Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows;
   mutable bool calcMassMatrix; //Internal flag
};
}

#include "Piro_InvertMassMatrixDecorator_Def.hpp"

#endif
