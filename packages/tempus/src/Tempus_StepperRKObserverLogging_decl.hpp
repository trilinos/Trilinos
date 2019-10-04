// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKObserverLogging_decl_hpp
#define Tempus_StepperRKObserverLogging_decl_hpp

#include "Tempus_StepperRKObserver.hpp"
#include <list>

namespace Tempus {

/** \brief This observer logs calls to observer functions.
 *  This observer simply logs and counts the calls to each of the
 *  observer functions.  This is useful in monirtoring and debugging
 *  the time integration.
 */
template<class Scalar>
class StepperRKObserverLogging
  : virtual public Tempus::StepperRKObserver<Scalar>
{
public:

  /// Constructor
  StepperRKObserverLogging();

  /// Destructor
  virtual ~StepperRKObserverLogging();

  /// \name Override StepperRKObserver basic methods
  //@{
  /// 1.) Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// 2.) Observe Stepper at beginning of each stage.
  virtual void observeBeginStage(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// 3.) Observe Stepper before Explicit evaluation of Implicit ODE ME (IMEX).
  virtual void observeBeforeImplicitExplicitly(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// 4.) Observe Stepper before nonlinear solve (DIRK/IMEX).
  virtual void observeBeforeSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// 5.) Observe Stepper after nonlinear solve (DIRK/IMEX).
  virtual void observeAfterSolve(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// 6.) Observe Stepper before Explicit evaluation of Implicit ODE ME (IMEX).
  virtual void observeBeforeExplicit(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// 7.) Observe Stepper at end of each stage.
  virtual void observeEndStage(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;

  /// 8.) Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */) override;
  //@}

  void resetLogCounters();

  Teuchos::RCP<const std::map<std::string,int> > getCounters();

  Teuchos::RCP<const std::list<std::string> > getOrder();

  /** \name String names logged in map
      Use these strings to validate a call stack with this observer
  */
  //@{
    const std::string nameObserveBeginTakeStep_; 
    const std::string nameObserveBeginStage_;              
    const std::string nameObserveBeforeImplicitExplicitly_;
    const std::string nameObserveBeforeSolve_;             
    const std::string nameObserveAfterSolve_;              
    const std::string nameObserveBeforeExplicit_;          
    const std::string nameObserveEndStage_;                
    const std::string nameObserveEndTakeStep_;             
  //@}

private:

  /** \brief Asserts next call on the stack is correct and removes from stack

      This is a const method so that it can be called from the
      derived StepperRKObserver methods that are const.
  */
  void logCall(const std::string call) const;

  Teuchos::RCP< std::map<std::string,int> > counters_;
  Teuchos::RCP< std::list<std::string> > order_;

};

} // namespace Tempus
#endif // Tempus_StepperRKObserverLogging_decl_hpp
