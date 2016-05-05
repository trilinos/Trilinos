#ifndef TEMPUS_STEPPERFORWARDEULER_HPP
#define TEMPUS_STEPPERFORWARDEULER_HPP

#include "Tempus_Stepper.hpp"

namespace Tempus {


/** \brief Forward Euler time stepper.
 *  Forward Euler is an explicit time stepper (i.e., no solver used).
 */
template<class Scalar>
class StepperForwardEuler : virtual public Tempus::Stepper<Scalar>
{
public:

  /// Constructor
  StepperForwardEuler(
    Teuchos::RCP<Teuchos::ParameterList>                pList_,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model_ );

  /// Destructor
//  virtual ~StepperForwardEuler();

  /// \name Basic stepper methods
  //@{
    /// Take the specified timestep, dt, and return true if successful.
    virtual bool takeStep(const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    virtual void setStepperState(
      const Teuchos::RCP<Tempus::StepperState<Scalar> >& stepperState);

    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getStepperState();
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

protected:

  Teuchos::RCP<Teuchos::ParameterList>               pList;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs;

  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState;
};
} // namespace Tempus

#include "Tempus_StepperForwardEuler_impl.hpp"

#endif // TEMPUS_STEPPERFORWARDEULER_HPP
