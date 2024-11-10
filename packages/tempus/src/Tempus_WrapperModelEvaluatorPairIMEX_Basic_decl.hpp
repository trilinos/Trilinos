//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorPairIMEX_Basic_decl_hpp
#define Tempus_ModelEvaluatorPairIMEX_Basic_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"

namespace Tempus {

/** \brief ModelEvaluator pair for implicit and explicit (IMEX) evaulations.
 *
 *  This ModelEvaluator takes a state, x, and determines the explicit and
 *  implicit residuals.  Additionally, it coordinates the explicit and
 *  implicit physics to ensure they are compatible, e.g.,
 *  how to translate between implicit and explicit model in and out
 *  arguments, if needed.
 *
 *  All functions called on WrapperModelEvaluatorPairIMEX_Basic will call
 *  the same function on the implicit Model Evaluator.  This was selected
 *  because the WrapperModelEvaluatorPairIMEX_Basic will be passed to the
 *  solvers which in turn make calls to solve the implicit ODE.
 *
 *  If the explicit version of the Model Evaluator functions are needed,
 *  one should directly call it through the explicit Model Evaluator, e.g.,
 *  getExplicitModel()->get_x_space().
 *
 *  This was taken and modified from Drekar's IMEXModelPair class.
 */
template <typename Scalar>
class WrapperModelEvaluatorPairIMEX_Basic
  : public Tempus::WrapperModelEvaluatorPairIMEX<Scalar> {
 public:
  /// Constructor
  WrapperModelEvaluatorPairIMEX_Basic(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel)
    : timeDer_(Teuchos::null)
  {
    setExplicitModel(explicitModel);
    setImplicitModel(implicitModel);
    initialize();
  }

  /// Destructor
  virtual ~WrapperModelEvaluatorPairIMEX_Basic() {}

  /// Initialize after setting member data.
  virtual void initialize();

  /// \name Overridden from Tempus::WrapperModelEvaluatorPairIMEX
  //@{
  virtual void setAppModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& me);
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getAppModel()
      const;

  /// Set parameters for application implicit ModelEvaluator solve.
  void setForSolve(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar> >& p,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& y = Teuchos::null,
      const int index                                   = -1 /* index and y are for IMEX_RK_Partition */)
  {
    x_       = x;
    xDot_    = xDot;
    time_    = time;
    p_       = p;
    y_       = y;
    index_   = index;
    timeDer_ = p->timeDer_;
  }

  //@}

  /// \name Methods that apply to both explicit and implicit terms.
  //@{
  /// Get the x-solution space
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space()
      const;

  /// Get the g space
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(
      int i) const;

  /// Get the p space
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(
      int i) const;
  //@}

  //@{ \name Accessors
  virtual void setExplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
  {
    explicitModel_ = model;
  }
  virtual void setImplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
  {
    implicitModel_ = model;
  }
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getExplicitModel()
      const
  {
    return explicitModel_;
  }
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getImplicitModel()
      const
  {
    return implicitModel_;
  }
  //@}

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
  virtual Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const
  {
    return implicitModel_->create_W_op();
  }

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const
  {
    return implicitModel_->get_W_factory();
  }

  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space()
      const
  {
    return explicitModel_->get_f_space();
  }

  virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& in,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& out) const;
  //@}

 protected:
  /// Default constructor -- only allowed for derived classes
  WrapperModelEvaluatorPairIMEX_Basic() {}

  /// Setup ME when using default constructor -- for derived classes
  void setup(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel)
  {
    setExplicitModel(explicitModel);
    setImplicitModel(implicitModel);
    initialize();
  }

 protected:
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > explicitModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > implicitModel_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot_;
  Scalar time_;
  Teuchos::RCP<ImplicitODEParameters<Scalar> > p_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > y_;
  int index_;
  Teuchos::RCP<TimeDerivative<Scalar> > timeDer_;
};

}  // namespace Tempus

#endif  // Tempus_ModelEvaluatorPairIMEX_Basic_decl_hpp
