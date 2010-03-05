#ifndef PIRO_THYRA_PERFORMANALYSIS
#define PIRO_THYRA_PERFORMANALYSIS

#include "Piro_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Piro {
  namespace Thyra {

  int PerformAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& analysisParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

  int PerformMoochoAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& moochoParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

  int PerformDakotaAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& dakotaParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

  int PerformOptiPackAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& optipackParams,
     Teuchos::ParameterList& globipackParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

   Teuchos::RCP<const Teuchos::ParameterList>
     getValidPiroAnalysisParameters();

  }
}

#endif //PIRO_THYRA_PERFORMANALYSIS
