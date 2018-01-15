// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepControlStrategy_hpp
#define Tempus_StepControlStrategy_hpp

#include "Tempus_SolutionHistory.hpp"


namespace Tempus {

template<class Scalar> class TimeStepControl;

/** \brief StepControlStrategy class for TimeStepControl
 *
 */
template<class Scalar>
class StepControlStrategy
{
public:

  /// Constructor
  StepControlStrategy(){}

  /// Destructor
  virtual ~StepControlStrategy(){}

  /// Observe Stepper at beginning of takeStep.
  virtual void getNextTimeStep(TimeStepControl<Scalar> tsc, Teuchos::RCP<SolutionHistory<Scalar> > sh ){}

};
} // namespace Tempus
#endif // Tempus_StepControlStrategy_hpp
