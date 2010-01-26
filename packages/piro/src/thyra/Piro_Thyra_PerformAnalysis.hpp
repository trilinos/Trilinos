#ifndef PIRO_THYRA_PERFORMANALYSIS
#define PIRO_THYRA_PERFORMANALYSIS

#include "Piro_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"

namespace Piro {
  namespace Thyra {

  int PerformAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& analysisParams
     );

  int PerformMoochoAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& moochoParams
     );

  int PerformDakotaAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& dakotaParams
     );

  int PerformOptiPackAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& optipackParams,
     Teuchos::ParameterList& globipackParams
     );
  }
}

#endif //PIRO_THYRA_PERFORMANALYSIS
