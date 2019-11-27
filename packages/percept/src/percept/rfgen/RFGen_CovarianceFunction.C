#include "RFGen_CovarianceFunction.h"

namespace RFGen
{

Teuchos::RCP<CovarianceFunction> 
buildCovarianceFunction(
  const unsigned covar_type,
  const int sdim,
  std::vector<double> cv_scales)
{
  Teuchos::RCP<CovarianceFunction> covarFunc;

  switch (covar_type)
  {
    case EXP_L2:
      covarFunc = Teuchos::rcp(
        new ExpMagL2CovarianceFunction(sdim, cv_scales));
      break;
    case EXP_L1:
      covarFunc = Teuchos::rcp(
        new ExpMagL1CovarianceFunction(sdim, cv_scales));
      break;
    case EXP_1D_L1:
      covarFunc = Teuchos::rcp(
        new ExpMag1DL1CovarianceFunction(sdim, cv_scales));
      break;
    default:
      //outputP0() << "Unknown covariance function type " << covar_type;
      TEUCHOS_ASSERT(false);
  }  

  return covarFunc;
}

}
