//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_WrapperModelEvaluatorSecondOrder_decl_hpp
#define Tempus_WrapperModelEvaluatorSecondOrder_decl_hpp

#include <functional>
#include "Tempus_config.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Thyra_VectorBase.hpp"

//#define VERBOSE_DEBUG_OUTPUT

namespace Tempus {

/** \brief A ModelEvaluator for residual evaluations given a state.
 *  This ModelEvaluator takes a state, x, and determines its residual,
 *  \f$ g(x) \f$, which is suitable for a nonlinear solve.  This is
 *  accomplished by computing the time derivative of the state, x_dot,
 *  (through Lambda functions), supplying the current time, and calling
 *  the application application ModelEvaluator, \f$ f(\dot{x},x,t) \f$.
 *
 *  This class breaks the primary design principle for ModelEvaluators;
 *  it is not stateless!
 */
template <typename Scalar>
class WrapperModelEvaluatorSecondOrder
  : public Tempus::WrapperModelEvaluator<Scalar> {
 public:
  typedef Thyra::VectorBase<Scalar> Vector;

  /// Constructor
  WrapperModelEvaluatorSecondOrder(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &appModel,
      const std::string schemeName)
    : appModel_(appModel), out_(Teuchos::VerboseObjectBase::getDefaultOStream())
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
    *out_ << "    schemeName = " << schemeName << "\n";
#endif
    if (schemeName == "Newmark Implicit a-Form" || schemeName == "HHT-Alpha") {
      schemeType_ = NEWMARK_IMPLICIT_AFORM;
    }
    else if (schemeName == "Newmark Implicit d-Form") {
      schemeType_ = NEWMARK_IMPLICIT_DFORM;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Error: WrapperModelEvaluatorSecondOrder called with unsopported "
              << "schemeName = " << schemeName
              << "!  Supported schemeNames are: 'Newmark Implicit a-Form' and "
              << "'HHT-Alpha'.\n");
    }
  }

  /// Set the underlying application ModelEvaluator
  void setAppModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &me)
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    appModel_ = me;
  }

  /// Get the underlying application model 'f'
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getAppModel() const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_;
  }

  /// Set values needed in evalModelImpl
  void initializeNewmark(Teuchos::RCP<const Vector> v_pred,
                         Teuchos::RCP<const Vector> d_pred, Scalar delta_t,
                         Scalar t, Scalar beta, Scalar gamma)
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    v_pred_  = v_pred;
    d_pred_  = d_pred;
    delta_t_ = delta_t;
    t_       = t;
    beta_    = beta;
    gamma_   = gamma;
  }

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_->create_W_op();
  }

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_->get_W_factory();
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_->get_f_space();
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int p) const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_->get_p_space(p);
  }

  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int p) const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_->get_p_names(p);
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_->get_x_space();
  }

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int i) const
  {
    return appModel_->get_g_space(i);
  }

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const
  {
#ifdef VERBOSE_DEBUG_OUTPUT
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
    return appModel_->getNominalValues();
  }

  /// Get InArgs the wrapper ModelEvalutor.
  virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> getInArgs()
  {
    return wrapperInArgs_;
  }

  /// Get OutArgs the wrapper ModelEvalutor.
  virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar> getOutArgs()
  {
    return wrapperOutArgs_;
  }

  /// Set parameters for application implicit ModelEvaluator solve.
  virtual void setForSolve(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > &x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > &xDot, const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar> > &p,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > &y = Teuchos::null,
      const int index                                   = -1 /* index and y are for IMEX_RK_Partition */)
  {
  }

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;
  //@}

  enum SCHEME_TYPE { NEWMARK_IMPLICIT_AFORM,
                     NEWMARK_IMPLICIT_DFORM };

 private:
  /// Default constructor - not allowed
  WrapperModelEvaluatorSecondOrder();

 private:
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;
  Scalar t_;
  Scalar gamma_;
  Scalar beta_;
  Scalar delta_t_;
  Teuchos::RCP<const Vector> d_pred_;
  Teuchos::RCP<const Vector> v_pred_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  SCHEME_TYPE schemeType_;

  Teuchos::RCP<TimeDerivative<Scalar> > timeDer_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> wrapperInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> wrapperOutArgs_;
};

}  // namespace Tempus

#endif  // Tempus_WrapperModelEvaluatorSecondOrder_hpp
