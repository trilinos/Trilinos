// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperSubcyclingObserver_hpp
#define Tempus_StepperSubcyclingObserver_hpp

#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration for recursive includes (this Observer <--> Stepper)
template<class Scalar> class StepperSubcycling;

/** \brief StepperSubcyclingObserver class for StepperSubcycling.
 *
 * This is a means for application developers to perform tasks
 * during the time steps, e.g.,
 *   - Compute specific quantities
 *   - Output information
 *   - "Massage" the working solution state
 *   - ...
 *
 * <b>Design Considerations</b>
 *   - StepperSubcyclingObserver is not stateless!  Developers may touch the
 *     solution state!  Developers need to be careful not to break the
 *     restart (checkpoint) capability.
 */
template<class Scalar>
class StepperSubcyclingObserver
 : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Constructor
  StepperSubcyclingObserver(){}

  /// Destructor
  virtual ~StepperSubcyclingObserver(){}

  /// Observe Stepper at beginning of takeStep.
  virtual void observeBeginTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}

  /// Observe Stepper at end of takeStep.
  virtual void observeEndTakeStep(
    Teuchos::RCP<SolutionHistory<Scalar> > /* sh */,
    Stepper<Scalar> & /* stepper */){}
};


} // namespace Tempus
#endif // Tempus_StepperSubcyclingObserver_hpp
