// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBackwardEulerVertexAppAction_hpp
#define Tempus_StepperBackwardEulerVertexAppAction_hpp

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

// Forward Declaration
template<class Scalar> class StepperBackwardEulerVertex;

/** \brief Application Action for StepperBackwardEulerVertex.
 *
 *  This class provides a means to apply various actions with the BackwardEulerVertex time step.
 *  The data available to this class is solution variables (through
 *  SolutionHistory), and stepper data (through the Stepper).  It allows
 *  the application to just observe this data, i.e., use but not change
 *  any of it (USER BEWARE!).
 *
 *  The locations for these AppAction calls
 *  (StepperBackwardEulerVertexAppAction::ACTION_LOCATION) are shown in the
 *  algorithm documentation of the StepperBackwardEulerVertex.
 */
template<class Scalar>
class StepperBackwardEulerVertexAppAction
{
public:

  /// Indicates the location of application action (see algorithm).
  enum ACTION_LOCATION {
    BEGIN_STEP,     ///< At the beginning of the step.
    BEFORE_SOLVE,   ///< Before the implicit solve.
    AFTER_SOLVE,    ///< After the implicit solve.
    END_STEP        ///< At the end of the step.
  };

  /// Constructor
  StepperBackwardEulerVertexAppAction(){}

  /// Destructor
  virtual ~StepperBackwardEulerVertexAppAction(){}

  /// Execute application action for BackwardEulerVertex Stepper.
  virtual void execute(
    Teuchos::RCP<SolutionHistory<Scalar> > sh,
    Teuchos::RCP<StepperBackwardEulerVertex<Scalar> > stepper,
    const typename StepperBackwardEulerVertexAppAction<Scalar>::ACTION_LOCATION actLoc) = 0;
};

} // namespace Tempus

#endif // Tempus_StepperBackwardEulerVertexAppAction_hpp
