// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKObserverComposite_decl_hpp
#define Tempus_StepperRKObserverComposite_decl_hpp

#include "Tempus_StepperRKObserver.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This observer is a composite observer,
 *
 *  which takes other StepperRKObservers and sequentially calls each
 *  individual observer function.
 *
 *  NOTE: certain RK observers (ERK,DIRK) methods execute 'back-to-back'
 *  without any intermediate code.
 */
template<class Scalar>
class StepperRKObserverComposite
  : virtual public Tempus::StepperRKObserver<Scalar>
{
public:

  /// Default constructor
  StepperRKObserverComposite();

  /// Destructor
  virtual ~StepperRKObserverComposite();

  /// \name Override StepperRKObserver basic methods
  //@{
  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper) override;

  /// Observe Stepper at beginning of each stage.
  virtual void observeBeginStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepperRK) override;

  /// Observe Stepper before Explicit evaluation of Implicit ODE ME.
  virtual void observeBeforeImplicitExplicitly(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// Observe Stepper before nonlinear solve.
  virtual void observeBeforeSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > sh ,
    Stepper<Scalar> & stepperRK) override;

  /// Observe Stepper after nonlinear solve.
  virtual void observeAfterSolve(
    Teuchos::RCP<SolutionHistory<Scalar> >  sh ,
    Stepper<Scalar> & stepperRK) override;

  /// Observe Stepper before  evaluation of Implicit ODE ME.
  virtual void observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepperRK) override;

  /// Observe Stepper at end of each stage.
  virtual void observeEndStage(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepperRK) override;

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper) override;

  // add observer to the composite observer list
  void addObserver(const Teuchos::RCP<StepperRKObserver<Scalar> > &observer);

  // clear all observer from the composite observer list
  void clearObservers();

  // get the number of RK stepper observers present in the composite
  std::size_t getSize() const { return observers_.size(); }
  //@}

private:

  std::vector<Teuchos::RCP<StepperRKObserver<Scalar > > > observers_;

};

} // namespace Tempus
#endif // Tempus_StepperRKObserverComposite_decl_hpp
