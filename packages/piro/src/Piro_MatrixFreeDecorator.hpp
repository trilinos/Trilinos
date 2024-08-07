// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_MATRIXFREEDECORATOR_HPP
#define PIRO_MATRIXFREEDECORATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

#include "Teuchos_RCP.hpp"

namespace Piro {

/** \brief Decorator class that creates a Jacobian (W) operator
 * using matrix-free directional derivatives.
 *
 * This class wraps a model evaluator supporting residual calculation
 * and adds support for the W matrix that implements a Thyra-based version of
 * NOX_Epetra_MatrixFree conforming to the Thyra::LinerarOp interface.
 *
 * This class supports time-dependent problems (characterized by x_dot != null)
 * and uses the input values of the alpha and beta coefficients.
 */

template <typename Scalar>
class MatrixFreeDecorator : public Thyra::ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers */
  //@{
  /** \brief . */
  explicit MatrixFreeDecorator(Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model, 
                               double lambda_ = 1.0e-6);
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  //@}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  //@}
  // Constant used in formulas to pick perturbation, typically 1.0e-6
  double lambda; 
};

}

#include "Piro_MatrixFreeDecorator_Def.hpp"

#endif /* PIRO_MATRIXFREEDECORATOR_HPP */
