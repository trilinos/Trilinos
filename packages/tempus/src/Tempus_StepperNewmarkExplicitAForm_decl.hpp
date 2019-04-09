// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkExplicitAForm_decl_hpp
#define Tempus_StepperNewmarkExplicitAForm_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExplicit.hpp"

namespace Tempus {


/** \brief Newmark Explicit time stepper.
 *
 *  This is the specific case of the more general Newmark time stepper
 *  where this stepper is explicit (\f$\beta = 0\f$) (i.e., no solver used).
 *
 *  The governing equation is solved by this stepper is
 *  \f[
 *    \mathbf{M}\, \ddot{\mathbf{x}} + \mathbf{C}\, \dot{\mathbf{x}}
 *    + \mathbf{K}\, \mathbf{x} = \mathbf{F}(t)
 *  \f]
 *  For the A-form (i.e., solving for the acceleration,
 *  \f$\mathbf{a} = \ddot{\mathbf{x}}\f$), we have the following explicit ODE
 *  \f[
 *     \mathbf{a} = -\mathbf{M}^{-1}\left[ \mathbf{C}\, \mathbf{v}
 *    + \mathbf{K}\, \mathbf{d} - \mathbf{F}(t) \right]
 *    = \bar{\mathbf{f}}(\mathbf{d}, \mathbf{v}, t)
 *  \f]
 *  where \f$\mathbf{v} = \dot{\mathbf{x}}\f$ and \f$\mathbf{d} = \mathbf{x}\f$.
 *
 *  <b> Algorithm </b>
 *  The algorithm for the Newmark explicit A-form is
 *   - if ( !useFSAL )
 *     - \f$\mathbf{a}^{n-1} =
 *          \bar{\mathbf{f}}(\mathbf{d}^{n-1}, \mathbf{v}^{n-1}, t^{n-1})\f$
 *   - \f$\mathbf{d}^{\ast} = \mathbf{d}^{n-1} + \Delta t \mathbf{v}^{n-1}
 *                            + \Delta t^2 \mathbf{a}^{n-1} / 2\f$
 *   - \f$\mathbf{v}^{\ast} =
 *        \mathbf{v}^{n-1} + \Delta t (1-\gamma) \mathbf{a}^{n-1}\f$
 *   - \f$\mathbf{a}^{\ast} =
 *        \bar{\mathbf{f}}(\mathbf{d}^{\ast}, \mathbf{v}^{\ast}, t^{n-1})\f$
 *   - \f$\mathbf{d}^n = \mathbf{d}^{\ast}\f$
 *   - \f$\mathbf{v}^n =
 *        \mathbf{v}^{\ast} + \Delta t \gamma \mathbf{a}^{\ast}\f$
 *   - if ( useFSAL )
 *     - \f$\mathbf{a}^n =
 *          \bar{\mathbf{f}}(\mathbf{d}^n, \mathbf{v}^n, t^n)\f$
 *
 *  The default for Forward Euler is to use FSAL (useFSAL=true).
 */
template<class Scalar>
class StepperNewmarkExplicitAForm
  : virtual public Tempus::StepperExplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Constructs with a default ParameterList.
   *  - Can reset ParameterList with setParameterList().
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
  */
  StepperNewmarkExplicitAForm();

  /// Constructor
  StepperNewmarkExplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > /* obs */ = Teuchos::null){}

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    virtual std::string getStepperType() const
     { return this->stepperPL_->template get<std::string>("Stepper Type"); }

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {
      if (gamma_ == 0.5) return 2.0;
      else return 1.0;
    }
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 2.0;}
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
      {return Scalar(1.0e+99);}

    virtual bool isExplicit()         const {return true;}
    virtual bool isImplicit()         const {return false;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}

    virtual OrderODE getOrderODE()   const {return SECOND_ORDER_ODE;}
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

protected:

  Scalar gamma_;

};
} // namespace Tempus

#endif // Tempus_StepperNewmarkExplicitAForm_decl_hpp
