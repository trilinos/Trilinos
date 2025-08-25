// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef MarkerUsingErrIndFraction_h
#define MarkerUsingErrIndFraction_h

#include <adapt/markers/Marker.hpp>

namespace percept{

//--------------------------------------------------------------------
/**
 *  @class  MarkerUsingErrIndFraction
 *  @author srkenno
 *  @brief  marking strategies using a fraction of the error indicator
 *    to be refined or unrefined, coupled with element budget
 *
 */
//--------------------------------------------------------------------

class MarkerUsingErrIndFraction : public Marker
{
public:

  MarkerUsingErrIndFraction(stk::mesh::BulkData& bulkData,  MarkerInfo & markerInfo);
  ~MarkerUsingErrIndFraction();

  // given a value between thresholdMin,Max, return an estimate of 
  //   how many new elements the mesh will have after refining based
  //   on the value of @param errIndRefineThreshold
  virtual size_t estimateNewElements(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec) override;
  virtual void getThresholdMinMax(double& min, double& max) override;

  virtual bool refine_element  (const double error) override {return (error > (markerInfo_.refineFraction_)  *(markerInfo_.maxErrorIndicator_));}
  virtual bool unrefine_element(const double error) override {return (error < (markerInfo_.unrefineFraction_)*(markerInfo_.maxErrorIndicator_));}

  // given the errIndRefFieldVec data, set its values based on the given
  //   errIndRefineThreshold
  virtual void markUsing(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec, bool do_refine) override;

};


} // namespace percept

#endif
