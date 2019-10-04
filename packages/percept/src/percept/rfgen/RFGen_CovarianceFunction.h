/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef RFGen_CovarianceFunction_h
#define RFGen_CovarianceFunction_h

#include "RFGen_Shards.h"

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

namespace RFGen
{

enum CovarianceType
{
  EXP_L2=0,
  EXP_L1,
  EXP_1D_L1
};

// 
// Abstract class for SCALAR covariance kernel function
//

class CovarianceFunction
{
public:
  explicit
  CovarianceFunction(
    const int spatialDim)
    : 
    spatialDim_(spatialDim)
    {}

  virtual ~CovarianceFunction()
  {}

  virtual void computeValues(
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x1,
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x2,
    shards::Array<double,shards::NaturalOrder,Point,Point> &result) const = 0;
  
protected:
  const int spatialDim_;
};

Teuchos::RCP<CovarianceFunction> 
buildCovarianceFunction(
  const unsigned covar_type,
  const int sdim,
  std::vector<double> cv_scales);

class ExpMagL2CovarianceFunction : public CovarianceFunction
{
public:
  explicit
  ExpMagL2CovarianceFunction(
    const int spatialDim,
    const std::vector<double> lengths)
    : CovarianceFunction(spatialDim),
    lengths_(lengths)
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(lengths[0]<=0, "Length scale should be positive");
    }

  virtual ~ExpMagL2CovarianceFunction()
  {}

  void computeValues(
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x1,
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x2,
    shards::Array<double,shards::NaturalOrder,Point,Point> &result) const
  {
    double mag;
    for (int q1=0; q1<x1.dimension(0); q1++)
    {
      for (int q2=0; q2<x2.dimension(0); q2++)
      {
	mag = 0.0;
	for (int d=0; d<spatialDim_; d++)
	{
	  mag += std::pow( ( x1(q1,d) - x2(q2,d) ) / lengths_[d], 2.0);
	}
	mag = std::sqrt(mag);
	
	result(q1,q2) = std::exp(-mag);
      }
    }
  }

 private:
  const std::vector<double> lengths_;
};

class ExpMagL1CovarianceFunction : public CovarianceFunction
{
public:
  explicit
  ExpMagL1CovarianceFunction(
    const int spatialDim,
    const std::vector<double> lengths)
    : CovarianceFunction(spatialDim),
    lengths_(lengths)
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(lengths[0]<=0, "Length scale should be positive");
    }

  virtual ~ExpMagL1CovarianceFunction()
  {}

  void computeValues(
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x1,
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x2,
    shards::Array<double,shards::NaturalOrder,Point,Point> &result) const
  {
    double temp;
    for (int q1=0; q1<x1.dimension(0); q1++)
    {
      for (int q2=0; q2<x2.dimension(0); q2++)
      {
	temp = 0.0;
	for (int d=0; d<spatialDim_; d++)
	{
	  temp += std::abs(x1(q1,d) - x2(q2,d)) / lengths_[d];
	}
	
	result(q1,q2) = std::exp(-temp);
      }
    }
  }

 private:
  const std::vector<double> lengths_;
};

class ExpMag1DL1CovarianceFunction : public CovarianceFunction
{
public:
  explicit
  ExpMag1DL1CovarianceFunction(
    const int spatialDim,
    const std::vector<double> lengths)
    : CovarianceFunction(spatialDim),
    lengths_(lengths)
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(lengths[0]<=0, "Length scale should be positive");
    }

  virtual ~ExpMag1DL1CovarianceFunction()
  {}

  void computeValues(
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x1,
    const shards::Array<double,shards::NaturalOrder,Point,Dim> &x2,
    shards::Array<double,shards::NaturalOrder,Point,Point> &result) const
  {
    double temp;
    for (int q1=0; q1<x1.dimension(0); q1++)
    {
      for (int q2=0; q2<x2.dimension(0); q2++)
      {
	temp = std::fabs(x1(q1,0) - x2(q2,0)) / lengths_[0];
	
	result(q1,q2) = std::exp(-temp);
      }
    }
  }

 private:
  const std::vector<double> lengths_;
};

}

#endif // RFGen_CovarianceFunction_h
