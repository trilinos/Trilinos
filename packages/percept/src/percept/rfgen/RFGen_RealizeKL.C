#include <cmath>

#include "RFGen_RealizeKL.h"

#include "Teuchos_Assert.hpp"

namespace RFGen
{

RealizeKL::RealizeKL(
  const unsigned rv_dim,
  const double rf_mean,
  const double rf_variance)
  : 
  rv_dim_(rv_dim),
  rf_mean_(rf_mean),
  rf_variance_(rf_variance)
{
  rf_coeffs_.resize(rv_dim);
}

void 
RealizeKL::computeRealization(
  const shards::Array<double,shards::NaturalOrder,Eigen> &rv_coeffs,
  const shards::Array<double,shards::NaturalOrder,Eigen> &kl_eigenValues,
  const shards::Array<double,shards::NaturalOrder,Cell,Eigen> &kl_eigenVectors,
  shards::Array<double,shards::NaturalOrder,Cell> &rf_values)
{
  // sanity checks
  TEUCHOS_TEST_FOR_EXCEPT(rv_coeffs.dimension(0)!=(int) rv_dim_);
  TEUCHOS_TEST_FOR_EXCEPT(kl_eigenValues.dimension(0)<(int) rv_dim_);
  TEUCHOS_TEST_FOR_EXCEPT(kl_eigenVectors.dimension(0)!=rf_values.dimension(0));
  TEUCHOS_TEST_FOR_EXCEPT(kl_eigenVectors.dimension(1)<(int) rv_dim_);

  for (unsigned j=0; j<rv_dim_; j++)
  {
    // sanity check
    TEUCHOS_TEST_FOR_EXCEPT(kl_eigenValues(j)<=0);
    rf_coeffs_[j] = rf_variance_ * std::sqrt(kl_eigenValues(j)) * rv_coeffs(j);
  }

  const unsigned numElem = kl_eigenVectors.dimension(0);
  for (unsigned i=0; i<numElem; i++)
  {
    rf_values(i) = rf_mean_;

    for (unsigned j=0; j<rv_dim_; j++)
      rf_values(i) += rf_coeffs_[j] * kl_eigenVectors(i,j);
  }
}

}
