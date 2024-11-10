//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_WrapperModelEvaluatorBasic_decl_hpp
#define Tempus_WrapperModelEvaluatorBasic_decl_hpp

#include <functional>
#include "Tempus_config.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"

namespace Tempus {

/** \brief A ModelEvaluator for residual evaluations given a state.
 *  This ModelEvaluator takes a state, x, and determines its residual,
 *  \f$ g(x) \f$, which is suitable for a nonlinear solve.  This is
 *  accomplished by computing the time derivative of the state, x_dot,
 *  (through Lambda functions), supplying the current time, and calling
 *  the application application ModelEvaluator, \f$ f(x,\dot{x},t) \f$.
 */
template <typename Scalar>
class WrapperModelEvaluatorBasic
  : public Tempus::WrapperModelEvaluator<Scalar> {
 public:
  /// Constructor
  WrapperModelEvaluatorBasic(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> &appModel)
    : appModel_(appModel), timeDer_(Teuchos::null), evaluationType_(SOLVE_FOR_X)
  {
    using Teuchos::rcp_const_cast;

    p_     = Teuchos::rcp(new ImplicitODEParameters<Scalar>());
    index_ = -1;

    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Scalar> inArgs = appModel_->getNominalValues();
    x_                         = rcp_const_cast<Thyra::VectorBase<Scalar>>(inArgs.get_x());

    if (inArgs.supports(MEB::IN_ARG_x_dot)) {
      xDot_ = rcp_const_cast<Thyra::VectorBase<Scalar>>(inArgs.get_x_dot());
    }
    else {
      xDot_ = Teuchos::null;
    }
  }

  /// Set the underlying application ModelEvaluator
  virtual void setAppModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> &me)
  {
    appModel_ = me;
  }

  /// Get the underlying application model 'f'
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> getAppModel() const
  {
    return appModel_;
  }

  /// Set parameters for application implicit ModelEvaluator solve.
  void setForSolve(
      const Teuchos::RCP<Thyra::VectorBase<Scalar>> &x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar>> &xDot, const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar>> &p,
      const Teuchos::RCP<Thyra::VectorBase<Scalar>> &y = Teuchos::null,
      const int index                                  = -1 /* index and y are for IMEX_RK_Partition */)
  {
    x_     = x;
    xDot_  = xDot;
    time_  = time;
    p_     = p;
    y_     = y;
    index_ = index;

    timeDer_        = p->timeDer_;
    evaluationType_ = p->evaluationType_;
  }

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> create_W_op() const
  {
    return appModel_->create_W_op();
  }

  Teuchos::RCP<Thyra::PreconditionerBase<Scalar>> create_W_prec() const
  {
    return appModel_->create_W_prec();
  }

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar>>
  get_W_factory() const
  {
    return appModel_->get_W_factory();
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_f_space() const
  {
    return appModel_->get_f_space();
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_p_space(int p) const
  {
    return appModel_->get_p_space(p);
  }

  Teuchos::RCP<const Teuchos::Array<std::string>> get_p_names(int p) const
  {
    return appModel_->get_p_names(p);
  }

  Teuchos::ArrayView<const std::string> get_g_names(int g) const
  {
    return appModel_->get_g_names(g);
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_x_space() const
  {
    return appModel_->get_x_space();
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_g_space(int i) const
  {
    return appModel_->get_g_space(i);
  }

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const
  {
    return appModel_->getNominalValues();
  }

  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> create_DfDp_op(int l) const
  {
    return appModel_->create_DfDp_op(l);
  }

  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> create_DgDx_op(int j) const
  {
    return appModel_->create_DgDx_op(j);
  }

  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> create_DgDp_op(int j, int l) const
  {
    return appModel_->create_DgDp_op(j, l);
  }

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;
  //@}

 private:
  /// Default constructor - not allowed
  WrapperModelEvaluatorBasic();

 private:
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> appModel_;

  Teuchos::RCP<Thyra::VectorBase<Scalar>> x_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> xDot_;
  Scalar time_;
  Teuchos::RCP<ImplicitODEParameters<Scalar>> p_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> y_;
  int index_;

  Teuchos::RCP<TimeDerivative<Scalar>> timeDer_;
  EVALUATION_TYPE evaluationType_;
};

}  // namespace Tempus

#endif  // Tempus_WrapperModelEvaluatorBasic_decl_hpp
