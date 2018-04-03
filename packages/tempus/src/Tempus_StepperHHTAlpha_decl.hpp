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
 * Here, we implement the HHT-Alpha scheme in predictor/corrector form;
 * see equations (10) and (13)-(19) in: G.M. Hulbert, J. Chung,
 * "Explicit time integration algorithms for structural dynamics with
 * optimal numerical dissipation", Comput. Methods Appl. Mech. Engrg.
 * 137 175-188 (1996).
 *
 * There are four parameters in the scheme: \f$\alpha_m\f$, \f$\alpha_f\f$,
 * \f$\beta\f$ and \f$\gamma\f$, all of which must be in the range \f$[0,1]\f$.
 * When \f$\alpha_m=\alpha_f = 0\f$, the scheme reduces to the Newmark Beta
 * scheme (see Tempus::StepperNewmark for details).  Like the Newmark Beta
 * scheme, the HHT-Alpha scheme can be either first or second order accurate,
 * and either explicit or implicit.
 *
 * Although the general form of the scheme has been implemented in Tempus,
 * it has only been verified for the case when \f$\alpha_m=\alpha_f = 0\f$
 * (corresponding to the Newmark Beta) scheme, so other values for these
 * parameters are not allowed at the present time.  Also, note that, like
 * the Newmark Beta stepper, the linear solve for the explicit version of
 * this scheme has not been optimized (the mass matrix is not lumped).
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

    void predictVelocity_alpha_f(Thyra::VectorBase<Scalar>& vPred,
                                 const Thyra::VectorBase<Scalar>& v) const;

    void predictDisplacement_alpha_f(Thyra::VectorBase<Scalar>& dPred,
                                     const Thyra::VectorBase<Scalar>& d) const;

    void correctAcceleration(Thyra::VectorBase<Scalar>& a_n_plus1,
                              const Thyra::VectorBase<Scalar>& a_n) const;

    void correctVelocity(Thyra::VectorBase<Scalar>& v,
                             const Thyra::VectorBase<Scalar>& vPred,
                             const Thyra::VectorBase<Scalar>& a,
                             const Scalar dt) const;

    void correctDisplacement(Thyra::VectorBase<Scalar>& d,
                               const Thyra::VectorBase<Scalar>& dPred,
                               const Thyra::VectorBase<Scalar>& a,
                               const Scalar dt) const;

private:

  /// Default Constructor -- not allowed
  StepperHHTAlpha();

private:

  Scalar beta_;
  Scalar gamma_;
  Scalar alpha_f_;
  Scalar alpha_m_;

  Teuchos::RCP<Teuchos::FancyOStream> out_;

};
} // namespace Tempus

#endif // Tempus_StepperHHTAlpha_decl_hpp
