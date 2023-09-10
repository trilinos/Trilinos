// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_BoundingBox.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <algorithm>

namespace krino{

template<class REAL, unsigned DIM>
void
BoundingBox_T<REAL,DIM>::gather_bboxes( const BoundingBox_T<REAL,DIM> & local_bbox,
                            std::vector < BoundingBox_T<REAL,DIM> > & all_bboxes )
{ /* %TRACE% */  /* %TRACE% */

  //
  // globally communicate the bounding box sizes for all procs
  //

  all_bboxes.clear();
  all_bboxes.reserve( stk::EnvData::parallel_size() );

  // put bbox data in Real vector for communication

  std::vector<Real> bbox_data( 2*DIM );
  for (unsigned i = 0, j = 0; i < DIM; i++ )
  {
    bbox_data[j++] = local_bbox.min[i];
    bbox_data[j++] = local_bbox.max[i];
  }

  std::vector<Real> bboxes_data( 2*DIM * stk::EnvData::parallel_size() );
  stk::parallel_vector_concat(stk::EnvData::parallel_comm(), bbox_data, bboxes_data);

  VecType remote_min, remote_max;
  for ( int i = 0; i < stk::EnvData::parallel_size(); i++ )
  {
    unsigned index = 2*DIM * i;
    for (unsigned j = 0; j < DIM; j++ )
    {
      remote_min[j] = bboxes_data[index++];
      remote_max[j] = bboxes_data[index++];
    }
    all_bboxes.emplace_back( remote_min, remote_max );
  }
}

template<class REAL, unsigned DIM>
void
BoundingBox_T<REAL,DIM>::global_reduce()
{
  const VecType local_min = min;
  const VecType local_max = max;

  stk::all_reduce_min(stk::EnvData::parallel_comm(), local_min.data(), min.data(), 3);
  stk::all_reduce_max(stk::EnvData::parallel_comm(), local_max.data(), max.data(), 3);
}

// Explicit template instantiation
template class BoundingBox_T<float,3>;
template class BoundingBox_T<float,2>;

} // namespace krino
