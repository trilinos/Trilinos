#ifndef __Piro_RythmosStepperFactory_hpp__
#define __Piro_RythmosStepperFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Rythmos_StepperBase.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

namespace Piro {

template <typename Scalar>
class RythmosStepperFactory {
public:

  virtual ~RythmosStepperFactory() {}

  virtual Teuchos::RCP<Rythmos::StepperBase<Scalar> > buildStepper(
                        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
                        const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > & solver,
                        const Teuchos::RCP<Teuchos::ParameterList> & paramList) = 0;
  
};

}

#endif
