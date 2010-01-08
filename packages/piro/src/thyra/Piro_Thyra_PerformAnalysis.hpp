#ifndef PIRO_THYRA_PERFORMANALYSIS
#define PIRO_THYRA_PERFORMANALYSIS

#include "Piro_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"

namespace Piro {
  namespace Thyra {

  void PerformAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& analysisParams
     );

  }
}

#endif //PIRO_THYRA_PERFORMANALYSIS
