// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRKObserverComposite_decl_hpp
#define Tempus_StepperExplicitRKObserverComposite_decl_hpp

#include "Tempus_StepperExplicitRKObserver.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This observer is a composite observer,
 *
 *  which takes other StepperExplicitRKObservers and sequentially calls each
 *  individual observer function.
 */
template<class Scalar>
class StepperExplicitRKObserverComposite
  : virtual public Tempus::StepperExplicitRKObserver<Scalar>
{
public:

  /// Default constructor
  StepperExplicitRKObserverComposite();

  /// Destructor
  virtual ~StepperExplicitRKObserverComposite();

  /// \name Override StepperExplicitRKObserver basic methods
  //@{
  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper) override;

  /// Observe Stepper at beginning of each stage.
  virtual void observeBeginStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK) override;

  /// Observe Stepper before Explicit evaluation of Implicit ODE ME.
  virtual void observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK) override;

  /// Observe Stepper at end of each stage.
  virtual void observeEndStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    StepperExplicitRK<Scalar> & stepperExplicitRK) override;

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper) override;

    // add observer to the composite observer list
    void addObserver(const Teuchos::RCP<StepperExplicitRKObserver<Scalar> > &observer);

    // clear all observer from the composite observer list
    void clearObservers();
  //@}

private:

  std::vector<Teuchos::RCP<StepperExplicitRKObserver<Scalar > > > observers_;

};

} // namespace Tempus
#endif // Tempus_StepperExplicitRKObserverComposite_decl_hpp
