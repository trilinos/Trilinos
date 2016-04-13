#ifndef TEMPUS_STEPPERFORWARDEULER_HPP
#define TEMPUS_STEPPERFORWARDEULER_HPP

#include "Tempus_Stepper.hpp"


namespace tempus {
template<class Scalar>
/** \brief Forward Euler time stepper.
 *  Forward Euler is an explicit time stepper (i.e., no solver used).
 */
class StepperForwardEuler : virtual public Stepper
{
  public:

    /// Constructor
    StepperForwardEuler(
      const RCP<Thyra::ModelEvaluator<Scalar> >& model);

    /// Destructor
    virtual ~StepperForwardEuler();
    //@}

    /// \name Basic stepper methods
    //@{
    /// Take the specified timestep, dt, and return true if successful.
    virtual bool takeStep(const Ptr<SolutionState<Scalar> >& workingState);

    //@}

    /// \name ParameterList methods
    //@{
    virtual void setParameterList(RCP<ParameterList> const& pl);
    virtual RCP<ParameterList> getNonconstParameterList();
    virtual RCP<ParameterList> unsetParameterList();
    virtual RCP<const ParameterList> getValidParameters() const;
    //@}

    /// \name Accessor methods
    //@{
    virtual std::string description() const;
    virtual void describe( Teuchos::FancyOStream        & out,
                           const Teuchos::EVerbosityLevel verbLevel) const;
    //@}

    /// \name Error estimation methods
    //@{

    //@}

    /// \name Observer methods
    //@{

    //@}

    /// \name Adjoint methods
    //@{
    //virtual Scalar takeAdjointStep();
    //@}

    /// \name Solution history methods
    //@{

    /// Functionality like InterpolationBuffer for multi-step methods, BDF.

    //@}

  protected:

    RCP<ParameterList>                        pList;
    RCP<const Thyra::ModelEvaluator<Scalar> > model;

    Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs;
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs;

};
} // namespace tempus
#endif // TEMPUS_STEPPERFORWARDEULER_HPP
