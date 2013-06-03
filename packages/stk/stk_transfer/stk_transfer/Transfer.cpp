/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <iostream>

#include <boost/shared_ptr.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearch.hpp>

#include <stk_transfer/Transfer.hpp>
#include <stk_transfer/TransferP2P.hpp>

namespace stk {
namespace transfer {

namespace {
stk::diag::TimerSet &HintsTimerSet()
{
  static stk::diag::TimerSet s_HintTimerSet(0);
  return s_HintTimerSet;
}
}

Transfer::Hints::Hints(stk::diag::Timer *timer) : 
      Tolerance      (0.1), 
      ExpansionFactor(1.5),
      Timer          (timer ? stk::diag::Timer("Transfer",*timer) : 
                              stk::diag::createRootTimer("Transfer", HintsTimerSet())){}


using namespace STK_TransferP2P;

Transfer::Transfer(const std::string &name,
                   const Hints &hints) :
  Name           (name),
  Tolerance      (hints.Tolerance) ,
  ExpansionFactor(hints.ExpansionFactor),
  RootTimer      ("Transfer Root",hints.Timer),
  GhostingTimer  ("Ghosting",     RootTimer)
   {}

template <unsigned DIM>
void Transfer::PToP(      MDArray &ToValues,
                    const MDArray &ToPoints,
                    const MDArray &FromValues,
                    const MDArray &FromPoints,
                    const stk::ParallelMachine  comm) {
  BoundingBox::Data radius = Tolerance;
  IdentProcRelation TotalRangeToDomain;
  PointMap to_points;
  convert_to_map<DIM>(to_points, ToPoints);
  while (to_points.size()) { // Keep going until all range points are processed.
    IdentProcRelation RangeToDomain;
    point_to_point_coarse_search<DIM>(RangeToDomain, to_points, FromPoints, radius, comm);
    filter_with_fine_search     <DIM>(RangeToDomain, to_points, FromPoints, comm);
    delete_range_points_found   <DIM>(to_points, RangeToDomain, comm);
    TotalRangeToDomain.insert(TotalRangeToDomain.end(),
                                   RangeToDomain.begin(),
                                   RangeToDomain.end());
    radius *= ExpansionFactor; // If points were missed, increase search radius.
  }
  sort(TotalRangeToDomain.begin(), TotalRangeToDomain.end());
  linear_interpolation <DIM>(ToValues,
                             FromValues,
                             TotalRangeToDomain,
                             ToPoints,
                             FromPoints,
                             comm);
}

template <unsigned DIM>
void Transfer::PPToP(      MDArray &ToValues,
                     const MDArray &ToPoints,
                     const MDArray &FromValues,
                     const MDArray &FromPoints,
                     const stk::ParallelMachine  comm) {

  BoundingBox::Data radius = Tolerance;
  IdentProcRelation TotalRangeToDomain;
  PointMap to_points;
  convert_to_map<DIM>(to_points, ToPoints);
  STKMesh FromMesh = convert_points_to_mesh(FromPoints, FromValues, comm);
  
  while (to_points.size()) { // Keep going until all range points are processed.
    IdentProcRelation RangeToDomain;
    point_to_point_coarse_search<DIM>(RangeToDomain, to_points, FromMesh, radius, comm);
    
    {
      stk::diag::TimeBlock __timer_ghosting(GhostingTimer);
      copy_domain_to_range_processors(FromMesh, RangeToDomain, Name, comm);
      
      const unsigned my_rank = parallel_machine_rank(comm);
      const IdentProcRelation::const_iterator end=RangeToDomain.end();
      for (IdentProcRelation::iterator i=RangeToDomain.begin(); i!=end; ++i) { 
        unsigned &domain_owning_proc = i->first.proc;
        const unsigned range_owning_rank  = i->second.proc;
        if (domain_owning_proc != my_rank && range_owning_rank == my_rank) {
          domain_owning_proc = my_rank;
        }   
      }
    }
    filter_with_fine_search     <DIM>(RangeToDomain, to_points, FromMesh, comm);
    delete_range_points_found<DIM>(to_points, RangeToDomain, comm);
    TotalRangeToDomain.insert(TotalRangeToDomain.end(),
                              RangeToDomain.begin(),
                              RangeToDomain.end());
    radius *= ExpansionFactor; // If points were missed, increase search radius.
  }
  sort(TotalRangeToDomain.begin(), TotalRangeToDomain.end());
  linear_interpolation <DIM>(ToValues,
                           FromMesh,
                           TotalRangeToDomain,
                           ToPoints,
                           FromMesh,
                           comm);
}


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

  if (3==Dim) {
    if(1==stk::parallel_machine_size(comm))
      PToP<3>(ToValues,ToPoints,FromValues,FromPoints,comm);
    else {
      PPToP<3>(ToValues,ToPoints,FromValues,FromPoints,comm);
    }
  } else {
    std::cerr <<__FILE__<<":"<<__LINE__<< " Only 3D is supported at this time. "<<std::endl;
  }
}




}}
