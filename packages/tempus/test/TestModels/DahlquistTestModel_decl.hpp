//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_DAHLQUIST_TEST_MODEL_DECL_HPP
#define TEMPUS_TEST_DAHLQUIST_TEST_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"               // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp"  // Implementation

namespace Tempus_Test {

/** \brief The classic Dahlquist Test Problem.
 *
 * This is the canonical Dahlquist test equation
 *   \f[
 *   \mathbf{\dot{x}} = \lambda \mathbf{x}.
 *   \f]
 * where \f$\lambda \in \mathbf{C}\f$ with the initial condition
 * \f$\mathbf{x}(0) = 1\f$.  The exact solution is
 *   \f[
 *   \mathbf{x}(t) = \exp{\lambda t}.
 *   \f]
 * When \f$\mathbf{R}(\lambda) < 0\f$, \f$ \mathbf{x} \rightarrow 0\f$ as
 * \f$t \rightarrow \infty\f$.  Numerical methods that exhibit this
 * property are said to be A-stable.
 */

template <class Scalar>
class DahlquistTestModel : public Thyra::StateFuncModelEvaluatorBase<Scalar> {
 public:
  // Default Constructor
  DahlquistTestModel();

  // Constructor
  DahlquistTestModel(Scalar lambda, bool includeXDot);

  void constructDahlquistTestModel(Scalar lambda, bool includeXDot);

  /// Default destructor
  ~DahlquistTestModel() = default;

  /// Exact solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

  Scalar getLambda() const { return lambda_; }

  /** \name Public functions overridden from ModelEvaluator. */
  //@{
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  //@}

 private:
  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar) const;
  //@}

 private:
  int dim_;
  int Np_;
  int np_;
  int Ng_;
  int ng_;
  bool haveIC_;
  bool acceptModelParams_;
  Scalar lambda_;
  bool includeXDot_;

  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > DxDp_space_;

  Scalar xIC_;     ///< Initial condition for x.
  Scalar xDotIC_;  ///< Initial condition for xDot.
};

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_DAHLQUIST_TEST_MODEL_DECL_HPP
