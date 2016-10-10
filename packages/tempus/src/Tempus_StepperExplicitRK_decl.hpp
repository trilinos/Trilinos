#ifndef Tempus_StepperExplicitRK_decl_hpp
#define Tempus_StepperExplicitRK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_RKButcherTableau.hpp"


namespace Tempus {


/** \brief Explicit Runge-Kutta time stepper.
 *  Explicit Runge-Kutta time stepper does not require any solver(s).
 */
template<class Scalar>
class StepperExplicitRK : virtual public Tempus::Stepper<Scalar>
{
public:

  /// Constructor
  StepperExplicitRK(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel);

  /// \name Basic stepper methods
  //@{
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
