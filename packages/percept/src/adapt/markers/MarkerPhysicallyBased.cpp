// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/markers/MarkerPhysicallyBased.hpp>

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


//==========================================================================
// Class Definition
//==========================================================================
// MarkerPhysicallyBased - computes refineField from errorIndicator field
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
  MarkerPhysicallyBased::MarkerPhysicallyBased(stk::mesh::BulkData & bulkData, MarkerInfo & markerInfo) 
    : Marker(bulkData,markerInfo), m_unrefinement_multiplier(markerInfo.physicalErrIndUnrefCriterionMultipler_)
{
  if (0.0 > m_unrefinement_multiplier || m_unrefinement_multiplier > 1.0)
    {
      throw std::runtime_error("MarkerPhysicallyBased: physicalErrIndUnrefCriterionMultipler should be in the range of [0,1]");
    }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MarkerPhysicallyBased::~MarkerPhysicallyBased()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- getThresholdMinMax-----------------------------------------------
//--------------------------------------------------------------------------

void MarkerPhysicallyBased::getThresholdMinMax(double& min, double& max)
{
  max = markerInfo_.maxErrorIndicator_;
  min = markerInfo_.physicalErrIndCriterion_;

  if (min > max) {
    min = max;  // no refinement done
    if (!bulkData_.parallel_rank()) std::cout << "Adapt: WARNING physical_error_criterion (" << markerInfo_.physicalErrIndCriterion_
                    << " is too large, thus no refinement will occur. "
                    << "\n   Consider setting it less than this current max value: " << max
                    << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- markUsing -------------------------------------------------------
//--------------------------------------------------------------------------

void MarkerPhysicallyBased::markUsing(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec, bool do_refine)
{
  for (unsigned i = 0; i < errIndRefFieldVec.size(); ++i) {
    // reset to 0, then check criterion
    *std::get<1>(errIndRefFieldVec[i]) = static_cast<percept::RefineFieldType::value_type>(0);
    if (*std::get<0>(errIndRefFieldVec[i]) < m_unrefinement_multiplier*errIndRefineThreshold) {
      *std::get<1>(errIndRefFieldVec[i]) = static_cast<percept::RefineFieldType::value_type>(-1);
    }
    if (do_refine && (*std::get<0>(errIndRefFieldVec[i]) > errIndRefineThreshold)
        && !std::get<2>(errIndRefFieldVec[i]))  {
      *std::get<1>(errIndRefFieldVec[i]) = static_cast<percept::RefineFieldType::value_type>(1);
    }
  }
}

//--------------------------------------------------------------------------
//-------- estimateNewElements ---------------------------------------------
//--------------------------------------------------------------------------

size_t MarkerPhysicallyBased::estimateNewElements(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec)
{
  size_t numRefinedLocal = 0, numRefinedGlobal = 0;
  for (size_t i = 0; i < errIndRefFieldVec.size(); ++i)
    {
      if ((*std::get<0>(errIndRefFieldVec[i]) > errIndRefineThreshold)
          && !std::get<2>(errIndRefFieldVec[i]) )
        ++numRefinedLocal;
    }

  stk::all_reduce_sum(bulkData_.parallel(), &numRefinedLocal, &numRefinedGlobal, 1);
  if (!bulkData_.parallel_rank()) std::cout << "Adapt: estimateNewElements, errIndRefineThreshold= " << errIndRefineThreshold
                          << " numRefinedLocal= " << numRefinedLocal << " numRefinedGlobal= " << numRefinedGlobal
                          << std::endl;
  return numRefinedGlobal;
}

} // namespace percept
