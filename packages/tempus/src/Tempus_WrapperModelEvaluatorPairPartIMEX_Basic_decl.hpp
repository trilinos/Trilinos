// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ModelEvaluatorPairPartIMEX_Basic_decl_hpp
#define Tempus_ModelEvaluatorPairPartIMEX_Basic_decl_hpp

#include "Tempus_WrapperModelEvaluatorPairIMEX.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"


namespace Tempus {

/** \brief ModelEvaluator pair for implicit and explicit (IMEX) evaulations.
 *
 *  All functions called on WrapperModelEvaluatorPairPartIMEX_Basic will call
 *  the same function on the implicit Model Evaluator.  This was selected
 *  because the WrapperModelEvaluatorPairPartIMEX_Basic will be passed to the
 *  solvers which in turn make calls to solve the implicit ODE.
 *
 *  If the explicit version of the Model Evaluator functions are needed,
 *  one should directly call it through the explicit Model Evaluator, e.g.,
 *  getExplicitModel()->get_x_space().
 *
 *  The one exception to this rule is for getNominalValues(), which is
 *  controlled by implicitNominalValues.  During the Integrator initialization
 *  this->getNominalValues needs to return
 *  explicitModel_->getNominalValues() [implicitNominalValues=false is the
 *  default], but during the nonlinear solves this->getNominalValues needs
 *  to return implicitModel_->getNominalValues() [implicitNominalValues=true].
 */
template <typename Scalar>
class WrapperModelEvaluatorPairPartIMEX_Basic
  : public Tempus::WrapperModelEvaluatorPairIMEX<Scalar>
{
public:

  /// Default constructor -- Still requires setting the models and running initialize.
  WrapperModelEvaluatorPairPartIMEX_Basic();

  /// Constructor
  WrapperModelEvaluatorPairPartIMEX_Basic(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel,
    int numExplicitOnlyBlocks = 0, int parameterIndex = -1);

  /// Destructor
  virtual ~WrapperModelEvaluatorPairPartIMEX_Basic(){}

  /// Initialize after setting member data.
  virtual void initialize();

  /// \name Overridden from Tempus::WrapperModelEvaluatorPairIMEX
  //@{
    virtual void setAppModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getAppModel() const;

    /// Set InArgs the wrapper ModelEvalutor.
    virtual void setInArgs(Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs)
    { wrapperImplicitInArgs_.setArgs(inArgs); }

    /// Get InArgs the wrapper ModelEvalutor.
    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> getInArgs()
    { return wrapperImplicitInArgs_; }

    /// Set OutArgs the wrapper ModelEvalutor.
    virtual void setOutArgs(Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs)
    { wrapperImplicitOutArgs_.setArgs(outArgs); }

    /// Get OutArgs the wrapper ModelEvalutor.
    virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar> getOutArgs()
    { return wrapperImplicitOutArgs_; }

    /// Set parameters for application implicit ModelEvaluator solve.
    virtual void setForSolve(Teuchos::RCP<TimeDerivative<Scalar> > timeDer,
      Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs,
      Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs,
      EVALUATION_TYPE /* evaluationType */ = SOLVE_FOR_X)
    {
      timeDer_ = timeDer;
      wrapperImplicitInArgs_.setArgs(inArgs);
      wrapperImplicitOutArgs_.setArgs(outArgs);
      useImplicitModel_ = true;
    }

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
  //@}

  //@{ \name Accessors
    virtual void setNumExplicitOnlyBlocks(int numExp)
    {numExplicitOnlyBlocks_ = numExp; }
    virtual void setExplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model )
    { explicitModel_ = model; }
    virtual void setImplicitModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model );
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getExplicitModel() const { return explicitModel_; }
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getImplicitModel() const { return implicitModel_; }
    virtual int getNumExplicitOnlyBlocks() const
    { return numExplicitOnlyBlocks_; }

    /// Extract IMEX vector from a full solution vector
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getIMEXVector(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & full) const;

    /// Extract IMEX vector for reading
    virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> > getIMEXVector(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & full) const;

    /// Extract explicit-only vector from a full solution vector
    virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getExplicitOnlyVector(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> > & full) const;

    /// Extract explicit-only vector for reading
    virtual Teuchos::RCP<const Thyra::VectorBase<Scalar> >getExplicitOnlyVector(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & full) const;

    /// Set the parameter index for explicit-only vector
    virtual void setParameterIndex(int parameterIndex = -1);
    /// Get the parameter index for explicit-only vector
    virtual int getParameterIndex() const { return parameterIndex_; }

    /// Set parameter to switch wrapperME base functions between explicit and implicit functions.
    virtual void setUseImplicitModel(bool tf) { useImplicitModel_ = tf; }
   /// Get parameter to switch wrapperME base functions between explicit and implicit functions.
    virtual bool getUseImplicitModel() const { return useImplicitModel_; }
  //@}

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
    virtual Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const
      { return implicitModel_->create_W_op(); }

    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      get_W_factory() const { return implicitModel_->get_W_factory(); }

    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_f_space() const;

    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
    virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
    virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar>createOutArgsImpl() const;

    virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> & in,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> & out) const;
  //@}

protected:

  /// Setup ME when using default constructor -- for derived classes
  void setup(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel,
    int numExplicitOnlyBlocks = 0, int parameterIndex = -1);

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > explicitModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > implicitModel_;

  Teuchos::RCP<TimeDerivative<Scalar> >              timeDer_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar>          wrapperImplicitInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>         wrapperImplicitOutArgs_;

  int numExplicitOnlyBlocks_;
  int parameterIndex_;    ///< implicit parameter index for explicit-only vector
  bool useImplicitModel_; ///< if true, use implicitModel_ else explicitModel_
};

} // namespace Tempus

#endif // Tempus_ModelEvaluatorPairPartIMEX_Basic_decl_hpp
