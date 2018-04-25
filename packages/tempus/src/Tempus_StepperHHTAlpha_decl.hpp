// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperHHTAlpha_decl_hpp
#define Tempus_StepperHHTAlpha_decl_hpp

#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorSecondOrder.hpp"

namespace Tempus {


/** \brief HHT-Alpha time stepper.
 *
 * Here, we implement the HHT-Alpha scheme; see : Hilber, H.M, Hughes,T.J.R and
 * Taylor, R.L. "Improved Numerical Dissipation for Time Integration Algorithms 
 * in Structural Dynamics" Earthquake Engineering and Structural Dynamics, 5:282-292
 * , 1977. in predictor/corrector form; 
 *
 * There are one parameters in the scheme: \f$\alpha\f$, which must be in the range
 *  \f$[0,1/3]\f$. When \f$\alpha= 0\f$, the scheme reduces to the Newmark Beta
 * scheme with no numerical damping.  \f$\alpha= 1/3\f$ gives max numerical
 * damping.
 */
 
template<class Scalar>
class StepperHHTAlpha : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor
  StepperHHTAlpha(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null){}

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {
      if (gamma_ == 0.5) return 2.0;
      else return 1.0;
    }
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 2.0;}

    virtual bool isExplicit()         const {return false;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
  //@}

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

    void predictVelocity(Thyra::VectorBase<Scalar>& vPred,
                             const Thyra::VectorBase<Scalar>& v,
                             const Thyra::VectorBase<Scalar>& a,
                             const Scalar dt) const;

    void predictDisplacement(Thyra::VectorBase<Scalar>& dPred,
                               const Thyra::VectorBase<Scalar>& d,
                               const Thyra::VectorBase<Scalar>& v,
                               const Thyra::VectorBase<Scalar>& a,
                               const Scalar dt) const;

    void correctVelocity(Thyra::VectorBase<Scalar>& v,
                             const Thyra::VectorBase<Scalar>& vPred,
                             const Thyra::VectorBase<Scalar>& a,
                             const Scalar dt) const;

    void correctDisplacement(Thyra::VectorBase<Scalar>& d,
                               const Thyra::VectorBase<Scalar>& v,
							   const Thyra::VectorBase<Scalar>& a_old,
                               const Thyra::VectorBase<Scalar>& a,
                               const Scalar dt) const;

private:

  /// Default Constructor -- not allowed
  StepperHHTAlpha();

private:

  Scalar beta_;
  Scalar gamma_;
  Scalar alpha_;

  Teuchos::RCP<Teuchos::FancyOStream> out_;

};
} // namespace Tempus

#endif // Tempus_StepperHHTAlpha_decl_hpp
