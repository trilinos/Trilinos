// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef MarkerPhysicallyBased_h
#define MarkerPhysicallyBased_h

#include <adapt/markers/Marker.hpp>

namespace percept{

//--------------------------------------------------------------------
/**
 *  @class  MarkerPhysicallyBased
 *  @author srkenno
 *  @brief  marking strategies that use a quantity like vorticity and
 *    a typical value (e.g. free-stream velocity) to choose elements
 *    for refinement - assumes error indicator reduces as elements are
 *    refined, thus the associated error indicator should have a length
 *    scale multiplier (eg. SimpleErrorIndicatorElemAlgorithm::VORTICITY_DX)
 */
//--------------------------------------------------------------------

class MarkerPhysicallyBased : public Marker
{
public:

  MarkerPhysicallyBased(stk::mesh::BulkData& bulkData,  MarkerInfo & markerInfo);
  ~MarkerPhysicallyBased();

  // given a value between thresholdMin,Max, return an estimate of
  //   how many new elements the mesh will have after refining based
  //   on the value of @param errIndRefineThreshold
  virtual size_t estimateNewElements(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec) override;
  virtual void getThresholdMinMax(double& min, double& max) override;

  virtual bool refine_element(const double error) override {return (error > markerInfo_.physicalErrIndCriterion_);}
  virtual bool unrefine_element(const double error) override {return (error < m_unrefinement_multiplier*(markerInfo_.physicalErrIndCriterion_));}

  // given the errIndRefFieldVec data, set its values based on the given
  //   errIndRefineThreshold
  virtual void markUsing(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec, bool do_refine) override;

  // multiplier for unrefinement 
  double m_unrefinement_multiplier;
};


} // namespace percept

#endif
