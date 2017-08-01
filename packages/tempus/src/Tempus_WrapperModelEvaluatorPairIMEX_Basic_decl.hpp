// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ModelEvaluatorIMEXPair_Basic_decl_hpp
#define Tempus_ModelEvaluatorIMEXPair_Basic_decl_hpp

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
  : public Tempus::WrapperModelEvaluatorPairIMEX<Scalar>
{
public:

  /// Constructor
  WrapperModelEvaluatorPairIMEX_Basic(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel)
    : explicitModel_(explicitModel), implicitModel_(implicitModel)
  {}


  /// Destructor
  virtual ~WrapperModelEvaluatorPairIMEX_Basic(){}

  /// \name Overridden from Tempus::WrapperModelEvaluator
  //@{
    virtual void setAppModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getAppModel() const;
    virtual void initialize(
      std::function<void (const Thyra::VectorBase<Scalar> &,
                                Thyra::VectorBase<Scalar> &)> computeXDot,
      double t, double alpha, double beta)
    { TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
       "Error - WrapperModelEvaluatorPairIMEX_Basic::initialize is not "
       "implemented yet!\n"); }
  //@}

  /// \name Methods that apply to both explicit and implicit terms.
  //@{
    /// Get the x-solution space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_x_space() const;

    /// Get the g space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_g_space(int i) const;

    /// Get the p space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_p_space(int i) const;

    /// Set values to compute x dot and evaluate application model.
    virtual void initialize2(
      std::function<void (const Thyra::VectorBase<Scalar> &,
                                Thyra::VectorBase<Scalar> &)> computeXDot,
                          Scalar ts, Scalar tHats,
                          Scalar alpha, Scalar beta,
                          Scalar stepSize, Scalar stageNumber)
    {
      computeXDot_ = computeXDot;
      ts_ = ts; tHats_ = tHats;
      alpha_ = alpha; beta_ = beta;
      stepSize_ = stepSize;
      stageNumber_ = stageNumber;
    }
  //@}

  //@{ \name Accessors
    virtual void setExplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model )
    { explicitModel_ = model; }
    virtual void setImplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model )
    { implicitModel_ = model; }
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getExplicitModel() const { return explicitModel_; }
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getImplicitModel() const { return implicitModel_; }
  //@}

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
    virtual Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const
      { return getImplicitModel()->create_W_op(); }

    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      get_W_factory() const { return getImplicitModel()->get_W_factory(); }

    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_f_space() const { return getImplicitModel()->get_f_space(); }

    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
    virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar>createOutArgsImpl() const;

    virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> & in,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> & out) const;
  //@}

  /// Evaluate the explicit model evaluator.
  virtual void evalExplicitModel(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X, Scalar time,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const;

  /// Evaluate the implicit model evaluator.
  virtual void evalImplicitModel(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;


private:
  /// Default constructor - not allowed
  WrapperModelEvaluatorPairIMEX_Basic(){}

protected:

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > explicitModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > implicitModel_;

  std::function<void (const Thyra::VectorBase<Scalar> &,
                            Thyra::VectorBase<Scalar> &)> computeXDot_;
  Scalar ts_;
  Scalar tHats_;
  Scalar alpha_;
  Scalar beta_;

  Scalar stepSize_;
  Scalar stageNumber_;
};

} // namespace Tempus

#endif // Tempus_ModelEvaluatorIMEXPair_Basic_decl_hpp
