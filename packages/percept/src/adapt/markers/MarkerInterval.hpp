// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef MarkerInterval_h
#define MarkerInterval_h

#include <adapt/markers/Marker.hpp>

namespace percept{

//--------------------------------------------------------------------
/**
 *  @class  MarkerInterval
 *  @author bcarnes
 *  @brief  mark using some interval between min/max element values of the error indicator.
 */
//--------------------------------------------------------------------

class MarkerInterval : public Marker
{
public:

  MarkerInterval(stk::mesh::BulkData& bulkData,  MarkerInfo & markerInfo);
  ~MarkerInterval();

  // given a value between thresholdMin,Max, return an estimate of
  //   how many new elements the mesh will have after refining based
  //   on the value of @param errIndRefineThreshold
  virtual size_t estimateNewElements(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec) override;
  virtual void getThresholdMinMax(double& min, double& max) override;

  virtual bool refine_element(const double error) override;
  virtual bool unrefine_element(const double error) override;

  // given the errIndRefFieldVec data, set its values based on the given
  //   errIndRefineThreshold
  virtual void markUsing(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec, bool do_refine) override;
};


} // namespace percept

#endif
