//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_STEADY_QUADRATIC_MODEL_DECL_HPP
#define TEMPUS_TEST_STEADY_QUADRATIC_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"               // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp"  // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief Simple quadratic equation with a stable steady-state.
 *
 * This is a simple differential equation
 *   \f[
 *   \mathbf{\dot{x}}=\mathbf{x}^2 - b^2
 *   \f]
 * which has steady state solutions \f$\mathbf{x} = \pm b\f$.  The solution
 * \f$\mathbf{x} = b\f$ is stable if \f$b < 0\f$ and the solution
 * \f$\mathbf{x} = -b\f$ is stable if \f$b > 0\f$.  This model is used to
 * test pseudo-transient sensitivity analysis methods.
 */
template <class Scalar>
class SteadyQuadraticModel : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
                             public Teuchos::ParameterListAcceptorDefaultBase {
 public:
  // Constructor
  SteadyQuadraticModel(
      Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  Scalar getSteadyStateSolution() const;

  // Exact sensitivity solution
  Scalar getSteadyStateSolutionSensitivity() const;

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

  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const &paramList);
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

 private:
  void setupInOutArgs_() const;

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar) const;
  //@}

 private:
  int dim_;                ///< Number of state unknowns (2)
  int Np_;                 ///< Number of parameter vectors (1)
  int np_;                 ///< Number of parameters in this vector (2)
  int Ng_;                 ///< Number of observation functions (1)
  int ng_;                 ///< Number of elements in this observation function (1)
  bool useDfDpAsTangent_;  ///< Treat DfDp OutArg as tangent (df/dx*dx/dp+df/dp)
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > DxDp_space_;

  // Parameters for the model
  Scalar b_;  ///< Model parameter
};

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_STEADY_QUADRATIC_MODEL_DECL_HPP
