// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef Marker_h
#define Marker_h

#include <vector>
#include <tuple>

#include <stk_mesh/base/FieldBase.hpp>
#include <percept/FieldTypes.hpp>

#include "Teuchos_RCP.hpp"

namespace percept{

class BoundingRegion;

//--------------------------------------------------------------------
/**
 *  @class  Marker
 *  @author srkenno
 *  @brief  marking strategies coupled with element budget
 *
 */
//--------------------------------------------------------------------
typedef std::tuple<double *, percept::RefineFieldType::value_type *, bool > ErrIndInfoTuple;

struct MarkerInfo
{
  MarkerInfo()
    :
    maxErrorIndicator_(0.0),
    minErrorIndicator_(0.0),
    errorIndicator_(0),
    refineField_(0),
    refineFieldOrig_(0),
    refineLevelField_(0),
    transitionElementField_(0),
    transitionElementField2d_(0),
    numInitialElements_(0),
    maxRefinementLevel_(0),
    useMarker_(false),
    maxMarkerIterations_(100),
    maxRefinementNumberOfElementsFraction_(1.0),
    physicalErrIndCriterion_(0.0),
    physicalErrIndUnrefCriterionMultipler_(1.0),
    refineFraction_(0.0),
    unrefineFraction_(0.0),
    intervalUpperFraction_(1.0),
    intervalLowerFraction_(0.0),
    debug_(false)
  {}

  double maxErrorIndicator_;
  double minErrorIndicator_;

  percept::ErrorFieldType *errorIndicator_;
  percept::RefineFieldType *refineField_;
  percept::RefineFieldType *refineFieldOrig_;
  percept::RefineLevelType *refineLevelField_;
  percept::TransitionElementType *transitionElementField_;
  percept::TransitionElementType *transitionElementField2d_;

  int numInitialElements_;
  int maxRefinementLevel_;
  bool useMarker_;
  unsigned maxMarkerIterations_;

  double maxRefinementNumberOfElementsFraction_;

  double physicalErrIndCriterion_;
  double physicalErrIndUnrefCriterionMultipler_;

  double refineFraction_;
  double unrefineFraction_;

  double intervalUpperFraction_;
  double intervalLowerFraction_;

  Teuchos::RCP<BoundingRegion> boundingRegion_;

  bool debug_;
};

class Marker
{
public:

  Marker(stk::mesh::BulkData& bulkData,  MarkerInfo & markerInfo, stk::mesh::Selector *globalSelector=0);
  virtual ~Marker();

  // selector used to possibly exclude blocks - if not set, this defaults to universal_part
  void setSelector(stk::mesh::Selector *sel);
  stk::mesh::Selector *getSelector();

  // given a value between thresholdMin,Max, return an estimate of
  //   how many new elements the mesh will have after refining based
  //   on the value of @param errIndRefineThreshold
  virtual size_t estimateNewElements(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec)=0;
  virtual void getThresholdMinMax(double& min, double& max)=0;

  // performs the full marking operations
  virtual void mark();

  virtual bool refine_element(const double error) = 0;
  virtual bool unrefine_element(const double error) = 0;

  // given the errIndRefFieldVec data, set its values based on the given
  //   errIndRefineThreshold
  virtual void markUsing(double errIndRefineThreshold, std::vector<ErrIndInfoTuple>& errIndRefFieldVec, bool do_refine)=0;

  void computeMaxErrorIndicator();
  void checkElementBudget(std::vector<ErrIndInfoTuple >& errIndRefFieldVec);
  stk::mesh::BucketVector const& get_active_buckets(stk::mesh::EntityRank rank,
                                                    const stk::mesh::Selector & selector ,
                                                    bool get_all) const;
  void zero_refine_field();
  size_t count_elements(stk::mesh::Selector sel, bool count_all_including_parents = false);

protected:

  stk::mesh::BulkData & bulkData_;
  MarkerInfo markerInfo_;
  bool doAdapt_;
  stk::mesh::Selector adapterSelector_[4];
  stk::mesh::Selector *m_globalSelector;  // used to choose which blocks to mark (or exclude)
};


} // namespace percept

#endif
