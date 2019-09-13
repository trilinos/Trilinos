/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#ifndef RFGen_RealizeKL_h
#define RFGen_RealizeKL_h

#include "RFGen_Shards.h"

namespace RFGen
{

class RealizeKL
{
public:
  explicit
  RealizeKL(
    const unsigned rv_dim,
    const double rf_mean,
    const double rf_variance);

  virtual ~RealizeKL()
  {}

  void computeRealization(
    const shards::Array<double,shards::NaturalOrder,Eigen> &rv_coeffs,
    const shards::Array<double,shards::NaturalOrder,Eigen> &kl_eigenValues,
    const shards::Array<double,shards::NaturalOrder,Cell,Eigen> &kl_eigenVectors,
    shards::Array<double,shards::NaturalOrder,Cell> &rf_values);

private:
  const unsigned rv_dim_;

  const double rf_mean_;
  const double rf_variance_;

  std::vector<double> rf_coeffs_;
};

}

#endif
