#ifndef __Piro_TempusStepControlFactory_hpp__
#define __Piro_TempusStepControlFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Tempus_Stepper.hpp"
#include "Tempus_TimeStepControl.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

namespace Piro {

template <typename Scalar>
class TempusStepControlFactory : public Tempus::TimeStepControl<Scalar>{

public:

  virtual ~TempusStepControlFactory() {}

  virtual Teuchos::RCP<Tempus::TimeStepControl<Scalar> > buildStepControl() = 0;

};

}

#endif
