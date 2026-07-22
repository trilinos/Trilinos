// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/markers/MarkerUsingErrIndFraction.hpp>

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


#include <algorithm>
#include <functional>

namespace percept{


//==========================================================================
// Class Definition
//==========================================================================
// MarkerUsingErrIndFraction - computes refineField from errorIndicator field
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MarkerUsingErrIndFraction::MarkerUsingErrIndFraction(stk::mesh::BulkData& bulkData,  MarkerInfo & markerInfo) : Marker(bulkData, markerInfo)
{
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MarkerUsingErrIndFraction::~MarkerUsingErrIndFraction()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- getThresholdMinMax -------------------------------------------------------
//--------------------------------------------------------------------------

void MarkerUsingErrIndFraction::getThresholdMinMax(double& min, double& max)
{
  max = 1.0;
  min = markerInfo_.refineFraction_;
}

//--------------------------------------------------------------------------
//-------- markUsing -------------------------------------------------------
//--------------------------------------------------------------------------

void MarkerUsingErrIndFraction::markUsing(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec, bool do_refine)
{
  const double UF = markerInfo_.unrefineFraction_;

  double g_maxErrorIndicator = markerInfo_.maxErrorIndicator_;

  for (unsigned i = 0; i < errIndRefFieldVec.size(); ++i) {
    // reset to 0, then check criterion
    *std::get<1>(errIndRefFieldVec[i]) = static_cast<percept::RefineFieldType::value_type>(0);
    if (*std::get<0>(errIndRefFieldVec[i]) < UF*g_maxErrorIndicator) {
      *std::get<1>(errIndRefFieldVec[i]) = static_cast<percept::RefineFieldType::value_type>(-1);
    }
    if (do_refine) {
      if ((*std::get<0>(errIndRefFieldVec[i]) >  errIndRefineThreshold*g_maxErrorIndicator)
          && !std::get<2>(errIndRefFieldVec[i]))  {
        *std::get<1>(errIndRefFieldVec[i]) = static_cast<percept::RefineFieldType::value_type>(1);
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- estimateNewElements ---------------------------------------------
//--------------------------------------------------------------------------

size_t MarkerUsingErrIndFraction::estimateNewElements(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec)
{
  double g_maxErrorIndicator = markerInfo_.maxErrorIndicator_;

  size_t numRefinedLocal = 0, numRefinedGlobal = 0;
  for (size_t i = 0; i < errIndRefFieldVec.size(); ++i)
    {
      if ((*std::get<0>(errIndRefFieldVec[i]) > errIndRefineThreshold*g_maxErrorIndicator)
          && !std::get<2>(errIndRefFieldVec[i]) )
        ++numRefinedLocal;
    }

  stk::all_reduce_sum(bulkData_.parallel(), &numRefinedLocal, &numRefinedGlobal, 1);
  return numRefinedGlobal;
}

} // namespace percept
