//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_NONAUTOSRC_MODEL_DECL_HPP
#define TEMPUS_TEST_NONAUTOSRC_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"               // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp"  // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief Scalar nonautosrc-source model problem.
  *
  * This model solves a single-species implicit ODE:
  *
  * \f[
  *   \dot{x}(t) = S(t)
  * \f]
  *
  * where
  *
  * \f[
  *   S(t) =
  *   \frac{\phi_\mathrm{amp}}{2}
  *   \left[
  *     \tanh(a t - b)
  *     -
  *     \tanh(\kappa t - \beta)
  *   \right]
  *   \frac{\gamma \sigma}{N_A}.
  * \f]
  *
  * The implicit residual is:
  *
  * \f[
  *   f = \dot{x} - S(t).
  * \f]
  */

template <class Scalar>
class NonAutoSrcModel : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
                    public Teuchos::ParameterListAcceptorDefaultBase {
 public:
  // Constructor
  NonAutoSrcModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

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

  //@}

  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const &paramList);
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

 private:
  void setupInOutArgs_() const;
  Scalar sourceFlux_(const Scalar time) const;
  Scalar nonautoSource_(const Scalar time) const;

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar) const;
  //@}


 protected:
  int dim_;  ///< Number of state unknowns, currently 1
  int Np_;   ///< Number of parameter vectors, currently 0
  int np_;   ///< Number of parameters, currently 0
  int Ng_;                  ///< Number of observation functions (1)
  int ng_;                  ///< Number of elements in this observation function (1)
  bool haveIC_;             ///< false => no nominal values are provided (default=true)
  bool acceptModelParams_;  ///< Changes inArgs to require parameters
  bool useDfDpAsTangent_;   ///< Treat DfDp OutArg as tangent (df/dx*dx/dp+df/dp)
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;

  Scalar t0_ic_;  ///< Time value where the initial condition is specified
  Scalar gamma_;
  Scalar xs_;
  Scalar avogadro_;

  Scalar flux_amp_;
  Scalar a_;
  Scalar b_;
  Scalar kappa_;
  Scalar flux_beta_;

  Scalar initial_amp_;
};

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_NONAUTOSRC_MODEL_DECL_HPP
