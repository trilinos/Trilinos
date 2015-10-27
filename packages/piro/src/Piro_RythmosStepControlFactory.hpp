#ifndef __Piro_RythmosStepControlFactory_hpp__
#define __Piro_RythmosStepControlFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"

namespace Piro {

template <typename Scalar>
class RythmosStepControlFactory : public Rythmos::ImplicitBDFStepperRampingStepControl<Scalar>{
public:

  virtual ~RythmosStepControlFactory() {}

  virtual Teuchos::RCP<Rythmos::StepControlStrategyBase<Scalar> > buildStepControl() = 0;

};

}

#endif
