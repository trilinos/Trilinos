#ifndef Tempus_StepperExplicitRK_decl_hpp
#define Tempus_StepperExplicitRK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_RKButcherTableau.hpp"


namespace Tempus {


/** \brief Explicit Runge-Kutta time stepper.
 *  The general explicit Runge-Kutta method for \f$s\f$-stages, can be
 *  written as
 *  \f[
 *    X_{i} = x_{n-1}
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
 *  more accurate solution. Note for explicit RK that \f$a_{ij}=0\f$ for
 *  \f$j \leq i\f$ and does not require any solves.
 */
template<class Scalar>
class StepperExplicitRK : virtual public Tempus::Stepper<Scalar>
{
public:

  /// Constructor
  StepperExplicitRK(
    Teuchos::RCP<Teuchos::ParameterList>                      pList,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return eODEModel_;}

    void setTableau(std::string stepperType = "");

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {return ERK_ButcherTableau_->order();}
    virtual Scalar getOrderMin() const {return ERK_ButcherTableau_->orderMin();}
    virtual Scalar getOrderMax() const {return ERK_ButcherTableau_->orderMax();}
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
  StepperExplicitRK();

  void explicitEvalModel(Teuchos::RCP<SolutionState<Scalar> > currentState);

protected:

  std::string                                        description_;
  Teuchos::RCP<Teuchos::ParameterList>               pList_;
  /// Explicit ODE ModelEvaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > eODEModel_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>          inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>         outArgs_;

  Teuchos::RCP<const RKButcherTableau<Scalar> >      ERK_ButcherTableau_;

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stagef_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >           stageX_;


  Teuchos::RCP<Thyra::VectorBase<Scalar> >           ee_;
};
} // namespace Tempus

#endif // Tempus_StepperExplicitRK_decl_hpp
