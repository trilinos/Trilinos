// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKBase_hpp
#define Tempus_StepperRKBase_hpp


namespace Tempus {

/** \brief Base class for Runge-Kutta methods, ExplicitRK, DIRK and IMEX.
 *
 */
template<class Scalar>
class StepperRKBase
{

public:

  virtual int getStageNumber() { return stageNumber_; }

  virtual void setStageNumber(int s) { stageNumber_ = s; }

protected:

  int stageNumber_;    //< The Runge-Kutta stage number.
};

} // namespace Tempus

#endif // Tempus_StepperRKBase_hpp
