/*
 * Akri_DistributedVector.cpp
 *
 *  Created on: Dec 21, 2025
 *      Author: drnoble
 */
#include <Akri_DistributedVector.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_util/util/ReportHandler.hpp>

namespace krino {

DistributedVector xpby(const DistributedVector & x, const double b, const DistributedVector & y)
{
  STK_ThrowAssert(x.size() == y.size() && x.local_size() == y.local_size());
  DistributedVector result = x;
  for (size_t i=0; i<x.size(); ++i)
        result[i] += b*y[i];
  return result;
}

DistributedVector scalar_times_vector(const double a, const DistributedVector & x)
{
  DistributedVector result = x;
  for (auto & entry : result)
    entry *= a;
  return result;
}

double Dot(const DistributedVector & x, const DistributedVector & y)
{
  STK_ThrowAssert(x.size() == y.size() && x.local_size() == y.local_size());
  double dot = 0.;
  for (size_t i=0; i<x.local_size(); ++i)
    dot += float(x[i]*y[i]); // reduced order addend to produce parallel consistency

  const double localDot = dot;
  if (x.comm() != stk::parallel_machine_null())
    stk::all_reduce_sum( x.comm(), &localDot, &dot, 1 );
  return dot;
}

DistributedVector vectorSubtract(const DistributedVector& x, const DistributedVector& y)
{
  STK_ThrowAssert(x.size() == y.size() && x.local_size() == y.local_size());
  DistributedVector result = x;
  for (size_t i=0; i<x.size(); ++i)
      result[i] -= y[i];
  return result;
}

}


