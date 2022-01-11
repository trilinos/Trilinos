// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_BoundingBox.hpp>
#include <Akri_DiagWriter.hpp>

#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>

#include <algorithm>

namespace krino{

template<class REAL, int DIM>
void
BoundingBox_T<REAL,DIM>::pad_epsilon()
{
  if (!valid()) return;

  const double eps = std::numeric_limits<Real>::epsilon();

  for (int i = 0; i < DIM; i++ )
  {
    min[i] -= std::abs(min[i])*eps;
    max[i] += std::abs(max[i])*eps;
  }
}

template<class REAL, int DIM>
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
  for (int i = 0, j = 0; i < DIM; i++ )
  {
    bbox_data[j++] = local_bbox.min[i];
    bbox_data[j++] = local_bbox.max[i];
  }

  std::vector<Real> bboxes_data( 2*DIM * stk::EnvData::parallel_size() );
  stk::parallel_vector_concat(stk::EnvData::parallel_comm(), bbox_data, bboxes_data);

  VecType remote_min, remote_max;
  for ( int i = 0; i < stk::EnvData::parallel_size(); i++ )
  {
    int index = 2*DIM * i;
    for (int j = 0; j < DIM; j++ )
    {
      remote_min[j] = bboxes_data[index++];
      remote_max[j] = bboxes_data[index++];
    }
    all_bboxes.emplace_back( remote_min, remote_max );
  }
}

template<class REAL, int DIM>
void
BoundingBox_T<REAL,DIM>::global_reduce()
{
  const VecType local_min = min;
  const VecType local_max = max;

  stk::all_reduce_min(stk::EnvData::parallel_comm(), local_min.data(), min.data(), 3);
  stk::all_reduce_max(stk::EnvData::parallel_comm(), local_max.data(), max.data(), 3);
}

template<class REAL, int DIM>
REAL
BoundingBox_T<REAL,DIM>::SqrDistLowerBnd( const VecType & pt ) const
{
  // make sure bbox is valid
  ThrowAssert( valid() );

  Real delta, SqrDist = 0.;

  for ( int i = 0; i < DIM; i++ )
  {
    if ( pt[i] < min[i] )
    {
      delta = min[i] - pt[i];
      SqrDist += delta * delta;
    }
    else if ( pt[i] > max[i] )
    {
      delta = pt[i] - max[i];
      SqrDist += delta * delta;
    }
  }
  return ( SqrDist );
}

template<class REAL, int DIM>
REAL
BoundingBox_T<REAL,DIM>::SqrDistUpperBnd( const VecType & pt ) const
{ /* %TRACE% */  /* %TRACE% */
// make sure bbox is valid
  ThrowAssert( valid() );

  // We are guaranteed that there is a point on the surface on each face of the
  // bounding box.  So we know that the upper bound for the distance to a face is
  // the distance to the farthest point on that face.  So the upper bound for this
  // bounding box is the minimum of the upper bounds for each face.  In other words,
  // the upper bound is the minimum distance to the farthest point on each face.

  VecType close_pt;
  VecType far_pt;

  for ( int i = 0; i < DIM; i++ )
  {
    if ( pt[i] < min[i] )
    {
      close_pt[i] = min[i];
      far_pt[i] = max[i];
    }
    else if ( pt[i] > max[i] )
    {
      close_pt[i] = max[i];
      far_pt[i] = min[i];
    }
    else
    {
      if (pt[i]-min[i] < max[i]-pt[i])
      {
        close_pt[i] = min[i];
        far_pt[i] = max[i];
      }
      else
      {
        close_pt[i] = max[i];
        far_pt[i] = min[i];
      }
    }
  }

  Real SqrDistMin;
  if (3 == DIM)
  {
    SqrDistMin = (pt-VecType(close_pt[0],far_pt[1],far_pt[2])).length_squared();
    SqrDistMin = std::min(SqrDistMin,(pt-VecType(far_pt[0],close_pt[1],far_pt[2])).length_squared());
    SqrDistMin = std::min(SqrDistMin,(pt-VecType(far_pt[0],far_pt[1],close_pt[2])).length_squared());
  }
  else
  {
    ThrowAssert(2 == DIM);
    const Real zero = 0.0;
    SqrDistMin = (pt-VecType(close_pt[0],far_pt[1],zero)).length_squared();
    SqrDistMin = std::min(SqrDistMin,(pt-VecType(far_pt[0],close_pt[1],zero)).length_squared());
  }
  return ( SqrDistMin );
}

template<class REAL, int DIM>
REAL
BoundingBox_T<REAL,DIM>::SqrDistUpperBnd( const BoundingBox_T<REAL,DIM> & pt_box ) const
{ /* %TRACE% */  /* %TRACE% */
  // This is somewhat conservative.  More aggressive methods might be able to reduce this estimate
  // while still always being an upper bound.
  // Here we just estimate the upper bound of this distance from any point in pt_bbox to the
  // surface contained in *this by the following:
  //   Loop the faces of *this.
  //   Find the maximum distance from any point in pt_bbox to any point of the face.
  //   Take the minimum of these maximum face distances.

  if( !pt_box.valid() )
    return std::numeric_limits<REAL>::max();

  // make sure bbox is valid
  ThrowAssert( valid() );

  Real SqrDistMin = 0.0;
  for ( int j = 0; j < DIM; j++ )
  {
    for ( int side = 0; side < 2; ++side)
    {
      Real SqrDistSideMin = 0.0;
      for ( int i = 0; i < DIM; i++ )
      {
        const Real delta = (i!=j) ? std::max(pt_box.max[i]-min[i],max[i]-pt_box.min[i]) :
            ((side==0) ? std::max(min[i]-pt_box.min[i], pt_box.max[i]-min[i]) : std::max(max[i]-pt_box.min[i], pt_box.max[i]-max[i]));
        ThrowAssert(delta >= 0.0);
        SqrDistSideMin += delta*delta;
      }
      if (SqrDistMin==0.0 || SqrDistSideMin<SqrDistMin)
      {
        SqrDistMin = SqrDistSideMin;
      }
    }
  }

  return ( SqrDistMin );
}

// Explicit template instantiation
template class BoundingBox_T<float,3>;
template class BoundingBox_T<float,2>;

} // namespace krino
