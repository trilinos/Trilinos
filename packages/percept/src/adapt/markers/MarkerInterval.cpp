// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/markers/MarkerInterval.hpp>

#include <percept/PerceptMesh.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <adapt/AdaptedMeshVerifier.hpp>

#include <algorithm>
#include <functional>

namespace percept{

MarkerInterval::MarkerInterval(stk::mesh::BulkData & bulkData, MarkerInfo & markerInfo) 
  : Marker(bulkData,markerInfo)
{
  if (markerInfo_.intervalLowerFraction_ > markerInfo_.intervalUpperFraction_ ||
      markerInfo_.intervalLowerFraction_ < 0.0 || markerInfo_.intervalLowerFraction_ > 1.0 ||
      markerInfo_.intervalUpperFraction_ < 0.0 || markerInfo_.intervalUpperFraction_ > 1.0)
    {
      throw std::runtime_error("MarkerInterval: upper and lower fractions should be in the range [0,1] with lower less than the upper fraction.");
    }
}

MarkerInterval::~MarkerInterval()
{}

void 
MarkerInterval::getThresholdMinMax(double& min, double& max)
{
  // refine when error indicator is within interval [a,b]
  //  a = min_error + lower_fraction * (max_error - min_error)
  //  b = min_error + upper_fraction * (max_error - min_error)

  min = markerInfo_.minErrorIndicator_ + (markerInfo_.intervalLowerFraction_)*(markerInfo_.maxErrorIndicator_ - markerInfo_.minErrorIndicator_);
  max = markerInfo_.minErrorIndicator_ + (markerInfo_.intervalUpperFraction_)*(markerInfo_.maxErrorIndicator_ - markerInfo_.minErrorIndicator_);
}

bool 
MarkerInterval::refine_element(const double error) {
  static double min,max;
  getThresholdMinMax(min, max);

  return (min <= error && error <= max);
}

bool 
MarkerInterval::unrefine_element(const double error) {
  static double min,max;
  getThresholdMinMax(min, max);

  return (error < min || max < error);
}

void 
MarkerInterval::markUsing(double /*errIndRefineThreshold*/, std::vector<ErrIndInfoTuple>& /*errIndRefFieldVec*/, bool /*do_refine*/)
{
  return;
}

size_t MarkerInterval::estimateNewElements(double /*errIndRefineThreshold*/, std::vector<ErrIndInfoTuple>& errIndRefFieldVec)
{
  size_t numRefinedLocal = 0, numRefinedGlobal = 0;
  for (size_t i = 0; i < errIndRefFieldVec.size(); ++i) {
    if (refine_element(*std::get<0>(errIndRefFieldVec[i])) 
        && !std::get<2>(errIndRefFieldVec[i]))
      ++numRefinedLocal;
  }

  stk::all_reduce_sum(bulkData_.parallel(), &numRefinedLocal, &numRefinedGlobal, 1);

  return numRefinedGlobal;
}

} // namespace percept
