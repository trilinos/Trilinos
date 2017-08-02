// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ModelEvaluatorIMEXPair_hpp
#define Tempus_ModelEvaluatorIMEXPair_hpp

//#include "Tempus_ModelEvaluatorPairIMEX.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"


namespace Tempus {

/** \brief ModelEvaluator pair for implicit and explicit (IMEX) evaulations.
 *
 *  This is an interface for a ModelEvaluator that takes a state, x,
 *  and determines the explicit and implicit residuals.
 *
 *  This was taken and modified from Drekar's IMEXModelPair class.
 */
template <typename Scalar>
class WrapperModelEvaluatorPairIMEX
  : public Tempus::WrapperModelEvaluator<Scalar>
{
public:

  /// \name Methods that apply to both explicit and implicit terms.
  //@{
    /// Get the x-solution space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_x_space() const = 0;

    /// Get the g space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_g_space(int i) const = 0;

    /// Get the p space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_p_space(int i) const = 0;

    virtual void initialize2(
      std::function<void (const Thyra::VectorBase<Scalar> &,
                                Thyra::VectorBase<Scalar> &)> computeXDot,
                          Scalar ts, Scalar tHats,
                          Scalar alpha, Scalar beta,
                          Scalar stepSize, Scalar stageNumber) = 0;
  //@}

  //@{ \name Accessors
    virtual void setExplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & ) = 0;
    virtual void setImplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & ) = 0;
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getExplicitModel() const = 0;
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getImplicitModel() const = 0;
  //@}

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
    virtual Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const = 0;

    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      get_W_factory() const = 0;

    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_f_space() const = 0;

    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar>
      getNominalValues() const = 0;

    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const = 0;

    virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar>
      createOutArgsImpl() const = 0;

    virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> & in,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> & out) const = 0;
  //@}

  /// Evaluate the explicit model evaluator.
  virtual void evalExplicitModel(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X, Scalar time,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const = 0;

  /// Evaluate the implicit model evaluator.
  virtual void evalImplicitModel(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const = 0;
};

} // namespace Tempus

#endif // Tempus_ModelEvaluatorIMEXPair_hpp
