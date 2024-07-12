// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_MATRIXFREELINEAROP_HPP
#define PIRO_MATRIXFREELINEAROP_HPP

#include "Thyra_LinearOpDefaultBase.hpp"

#include "Thyra_ModelEvaluatorBase.hpp"
#include "Thyra_VectorBase.hpp"

namespace Piro {

/** \brief This class implements a matrix-free Jacobian linear operator
 * based on finite difference.
 *
 * This class wraps a model evaluator supporting residual calculation,
 * and computes a finite difference approximation of the
 * Jacobian directional derivative at the chosen base point.
 * It is the Thyra-based functional equivalent of Piro::Epetra::MatrixFreeOperator,
 * and implements the Thyra::LinerarOp interface instead of Epetra_Operator.
 *
 * Time-dependent problems (characterized by x_dot != null) are supported.
 * The class uses the input values of the alpha and beta coefficients
 * that appear in the definition of the Jacobian operator:
 * W = alpha * (Df/ Dx_dot) + beta * (Df / Dx)
 */

template <typename Scalar>
class MatrixFreeLinearOp : public Thyra::LinearOpDefaultBase<Scalar>
{
public:

  /** \name Constructors/initializers */
  //@{
  /** \brief Construct a partially initialized Jacobian operator
   *         for the specified model evaluator.
   *
   * @pre <tt>model != null</tt>*/
  explicit MatrixFreeLinearOp(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model, 
                              const double lambda);
  //@}

  /** \name Overridden from Thyra::LinearOpBase< Scalar >*/
  //@{
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > range() const;

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > domain() const;
  //@}

  /** \name Underlying model evaluator */
  //@{
  /** \brief Model evaluator whose Jacobian operator is implemented by this object */
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model() const;
  //@}

  /** \name Base point at which the Jacobian operator is computed */
  //@{
  /** \brief Input arguments defining the base point. */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint() const;

  /** \brief Residual evaluated at the base point.
   *
   * The object is fully initialized if and only if <tt>f_base != null</tt>. */
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > f_base() const;

  /** \brief Set the base point.
   *
   * This function must be called after construction to fully initialize the object.
   *
   * The object makes shallow copies of its arguments.
   * The user is responsible for avoiding unwanted aliasing.
   * Also, if the passed pointers do not have (possibly shared) ownership
   * of the references objects, the user must ensure that they are kept alive
   * until this function is called again, or the object is destroyed.
   * */
  void setBase(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> > &f_base);
  //@}

protected:
  /** \name Overridden from Thyra::LinearOpBase< Scalar >*/
  //@{
  /** \brief . */
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const;

  /** \brief . */
  void applyImpl(
      const Thyra::EOpTransp M_trans,
      const Thyra::MultiVectorBase<Scalar> &X,
      const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y,
      const Scalar alpha,
      const Scalar beta) const;
  //@}

private:
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > f_base_;
  
  // Constant used in formulas to pick perturbation, typically 1.0e-6
  double lambda_; 
};

}

#include "Piro_MatrixFreeLinearOp_Def.hpp"

#endif /* PIRO_MATRIXFREELINEAROP_HPP */
