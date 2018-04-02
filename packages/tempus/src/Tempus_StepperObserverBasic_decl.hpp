// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperObserverBasic_decl_hpp
#define Tempus_StepperObserverBasic_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperObserver.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
template<class Scalar> class Stepper;

/** \brief StepperObserverBasic class for Stepper class.
 */
template<class Scalar>
class StepperObserverBasic : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Constructor
  StepperObserverBasic();

  /// Destructor
  virtual ~StepperObserverBasic();

  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper);

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Stepper<Scalar> & stepper);
};

} // namespace Tempus

#endif // Tempus_StepperObserverBasic_decl_hpp
