#ifndef __Piro_RythmosStepControlFactory_hpp__
#define __Piro_RythmosStepControlFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"
//#include "Rythmos_StepControlStrategyAcceptingStepperBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"

namespace Piro {

template <typename Scalar>
class RythmosStepControlFactory : public Rythmos::ImplicitBDFStepperRampingStepControl<Scalar>{
//class RythmosStepControlFactory : public Rythmos::StepControlStrategyState<Scalar>{
public:

  virtual ~RythmosStepControlFactory() {}

  virtual Teuchos::RCP<Rythmos::StepControlStrategyBase<Scalar> > buildStepControl() = 0;
/*
  virtual void setRequestedStepSize(const StepperBase<Scalar>& stepper,
                                    const Scalar& stepSize,
                                    const StepSizeType& stepSizeType) = 0;
*/  
};

}

#endif
