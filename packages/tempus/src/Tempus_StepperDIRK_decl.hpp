#ifndef Tempus_StepperDIRK_decl_hpp
#define Tempus_StepperDIRK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_ResidualModelEvaluator.hpp"

namespace Tempus {


/** \brief Diagonally Implicit Runge-Kutta (DIRK) time stepper.
 *  The general DIRK method for \f$s\f$-stages, can be written as
 *  \f[
 *    X_{i} - \Delta t\, a\,\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t) = x_{n-1}
 *    + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\bar{f}(X_{j},t_{n-1}+c_{j}\Delta t)
 *  \f]
 *  \f[
 *    x_{n} = x_{n-1}
 *    + \Delta t\,\sum_{i=1}^{s}b_{i}\,\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)
 *  \f]
 *  where \f$\dot{x}=\bar{f}(x,t)\f$ is the explicit form of the
 *  ODE, \f$X_{i}\f$ are intermediate approximations to the solution
 *  at times, \f$t_{n-1}+c_{i}\Delta t\f$, (stage solutions) which may
 *  be correct to a lower order of accuracy than the solution, \f$x_{n}\f$.
 *  We should note that these lower-order approximations are combined
 *  through \f$b_{i}\f$ so that error terms cancel out and produce a
 *  more accurate solution. Note for DIRK methods that \f$a_{ii}=a\f$
 *  for all the diagonal components.  This is also referred to as
 *  Singly Diagonal Implicit Runge-Kutta (SDIRK) methods.
 *
 *  Note that the stage time derivatives,
 *  \f[
 *    \dot{X}_{i} = \bar{f}(X_{i},t_{n-1}+c_{i}\Delta t),
 *  \f]
 *  can be found via
 *  \f[
 *    \dot{X}_{i} = \frac{1}{\Delta t a_{ii}} [ X_{i} - x_{n-1}
 *                  - \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\dot{X}_{j} ]
 *  \f]
 */
template<class Scalar>
class StepperDIRK : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor
  StepperDIRK(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return residualModel_->getTransientModel();}

    void setSolver(
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &
      solver = Teuchos::null);

    void setTableau(std::string stepperType = "");

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> >getDefaultStepperState();
    virtual Scalar getOrder() const{return DIRK_ButcherTableau_->order();}
    virtual Scalar getOrderMin() const{return DIRK_ButcherTableau_->orderMin();}
    virtual Scalar getOrderMax() const{return DIRK_ButcherTableau_->orderMax();}
  //@}

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

private:

  /// Default Constructor -- not allowed
  StepperDIRK();

protected:

  std::string                                       description_;
  Teuchos::RCP<Teuchos::ParameterList>              pList_;
  Teuchos::RCP<ResidualModelEvaluator<Scalar> >     residualModel_;
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>         inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>        outArgs_;

  Teuchos::RCP<const RKButcherTableau<Scalar> >     DIRK_ButcherTableau_;

  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> >     stageXDot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >          stageX_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >          stageXPartial_;

  // Compute the balancing time derivative as a function of x
  std::function<void (const Thyra::VectorBase<Scalar> &,
                            Thyra::VectorBase<Scalar> &)>
  xDotFunction(Scalar s,Teuchos::RCP<const Thyra::VectorBase<Scalar> > stageX);

  Teuchos::RCP<Thyra::VectorBase<Scalar> >          ee_;
};
} // namespace Tempus

#endif // Tempus_StepperDIRK_decl_hpp
