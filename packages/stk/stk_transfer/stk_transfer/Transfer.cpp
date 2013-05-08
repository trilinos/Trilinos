/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <stk_transfer/Transfer.hpp>
#include <stk_transfer/TransferP2P.hpp>

namespace stk {
namespace transfer {

using namespace STK_TransferP2P;

Transfer::Transfer(const Hints &hints) :
  Tolerance      (hints.Tolerance) ,
  ExpansionFactor(hints.ExpansionFactor) 
   {}

void Transfer::PointToPoint(      MDArray &ToValues,
                            const MDArray &ToPoints,
                            const MDArray &FromValues,
                            const MDArray &FromPoints,
                            const stk::ParallelMachine  comm) {
  const int Dim = ToPoints.dimension(1);
  
  if (FromPoints.dimension(0) < Dim+1) 
    std::cerr <<__FILE__<<":"<<__LINE__<< " Not enough points in domain. "<<std::endl;
  if (FromPoints.dimension(1) != Dim) 
    std::cerr <<__FILE__<<":"<<__LINE__<< " Inconsistant spatial dimension. "<<std::endl;
  if (FromPoints.dimension(0) != FromValues.dimension(0)) 
    std::cerr <<__FILE__<<":"<<__LINE__<< " Inconsistant domain data. "<<std::endl;
  if (ToPoints.dimension(0) != ToValues.dimension(0)) 
    std::cerr <<__FILE__<<":"<<__LINE__<< " Inconsistant range data. "<<std::endl;
  if (ToValues.dimension(1) != FromValues.dimension(1)) 
    std::cerr <<__FILE__<<":"<<__LINE__<< " Inconsistant value data. "<<std::endl;

  BoundingBox::Data radius = Tolerance;
  if (3==Dim) {
    IdentProcRelation TotalRangeToDomain;
    PointMap to_points;
    convert_to_map<3>(to_points, ToPoints);
    while (to_points.size()) { // Keep going until all range points are processed.
      IdentProcRelation RangeToDomain;
      point_to_point_coarse_search<3>(RangeToDomain, to_points, FromPoints, radius, comm);
      filter_with_fine_search     <3>(RangeToDomain, to_points, FromPoints, comm);
      delete_range_points_found   <3>(to_points, RangeToDomain, comm);
      TotalRangeToDomain.insert(TotalRangeToDomain.end(),
                                     RangeToDomain.begin(),
                                     RangeToDomain.end());
      radius *= ExpansionFactor; // If points were missed, increase search radius.
    }
    sort(TotalRangeToDomain.begin(), TotalRangeToDomain.end());
    linear_interpolation <3>(ToValues,FromValues,TotalRangeToDomain,ToPoints,FromPoints,comm);
  } else {
    std::cerr <<__FILE__<<":"<<__LINE__<< " Only 3D is supported at this time. "<<std::endl;
  }
}

}}
