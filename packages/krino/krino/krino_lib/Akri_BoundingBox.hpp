// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_BoundingBox_h
#define Akri_BoundingBox_h

#include <Akri_Vec.hpp>

#include <stk_util/util/ReportHandler.hpp>

#include <limits>
#include <vector>

namespace krino {

template<class REAL, int DIM>
class BoundingBox_T {

public:
  typedef REAL Real;
  typedef Vec<Real,3> VecType;

private:
  VecType min;
  VecType max;

public:
  static void gather_bboxes( const BoundingBox_T<REAL,DIM> & local_bbox,
      std::vector< BoundingBox_T<REAL,DIM> > & all_bboxes );

  void global_reduce();

  template<class VECTYPE>
  void accommodate( const VECTYPE & pt );

  template<class REALARG>
  void accommodate( const BoundingBox_T<REALARG,DIM> & bbox );

  const VecType & get_min() const { return min; }
  const VecType & get_max() const { return max; }

  template<class VECTYPE>
  bool contains( const VECTYPE & pt ) const;

  // Lower bound on the square of the distance from the point pt to surfaces contained in *this
  Real SqrDistLowerBnd( const VecType & pt ) const;

  // Upper bound on the square of the distance from the point pt to surfaces contained in *this
  Real SqrDistUpperBnd( const VecType & pt ) const;

  // Upper bound on the square of the distance from any point in pt_bbox to surfaces contained in *this
  Real SqrDistUpperBnd( const BoundingBox_T<REAL,DIM> & pt_box ) const;

  void pad( const Real & dist );
  void pad_epsilon();

  void scale( const Real & fraction );

  Real SqrSize() const;

  bool valid() const
  {
    bool is_not_valid = (min[0]>max[0]);
    return !is_not_valid;
  }
  
  bool intersects( const BoundingBox_T<REAL,DIM> & bbox ) const;
  bool contains( const BoundingBox_T<REAL,DIM> & bbox ) const;
  
  int max_span_direction() const;

  VecType center() const { return Real(0.5)*(min+max); }
  
  void clear() {
    min = VecType(std::numeric_limits<Real>::max(),std::numeric_limits<Real>::max(),std::numeric_limits<Real>::max());
    max = VecType(-std::numeric_limits<Real>::max(),-std::numeric_limits<Real>::max(),-std::numeric_limits<Real>::max());
  }

  template<class VECTYPE>
  BoundingBox_T( const VECTYPE & pt_min,
	       const VECTYPE & pt_max )
  : min(pt_min[0], pt_min[1], pt_min[2]),
    max(pt_max[0], pt_max[1], pt_max[2]) {}

  BoundingBox_T()
  : min(std::numeric_limits<Real>::max(),std::numeric_limits<Real>::max(),std::numeric_limits<Real>::max()),
    max(-std::numeric_limits<Real>::max(),-std::numeric_limits<Real>::max(),-std::numeric_limits<Real>::max()) {}
};

template<class REAL, int DIM>
template<class VECTYPE>
inline void
BoundingBox_T<REAL,DIM>::accommodate( const VECTYPE & pt )
{
  for ( int i = 0; i < DIM; ++i )
  {
    if ( pt[i] < min[i] ) min[i] = pt[i];
    if ( pt[i] > max[i] ) max[i] = pt[i];
  }
}

template<class REAL, int DIM>
template<class REALARG>
inline void
BoundingBox_T<REAL,DIM>::accommodate( const BoundingBox_T<REALARG,DIM> & bbox )
{
  for ( int i = 0; i < DIM; ++i )
  {
    if ( bbox.get_min()[i] < min[i] ) min[i] = bbox.get_min()[i];
    if ( bbox.get_max()[i] > max[i] ) max[i] = bbox.get_max()[i];
  }
}

template<class REAL, int DIM>
template<class VECTYPE>
inline bool
BoundingBox_T<REAL,DIM>::contains( const VECTYPE & pt ) const
{
  for ( int i = 0; i < DIM; ++i )
  {
    if ( pt[i] < min[i] || pt[i] > max[i] ) return false;
  }
  return true;
}

template<class REAL, int DIM>
inline bool
BoundingBox_T<REAL,DIM>::intersects( const BoundingBox_T<REAL,DIM> & bbox ) const
{
  for ( int i = 0; i < DIM; ++i )
  {
    if ( bbox.max[i] < min[i] || bbox.min[i] > max[i] ) return false;
  }
  return true;
}

template<class REAL, int DIM>
inline bool
BoundingBox_T<REAL,DIM>::contains( const BoundingBox_T<REAL,DIM> & bbox ) const
{
  for ( int i = 0; i < DIM; ++i )
  {
    if ( bbox.min[i] < min[i] || bbox.max[i] > max[i] ) return false;
  }
  return true;
}

template<class REAL, int DIM>
inline int
BoundingBox_T<REAL,DIM>::max_span_direction() const
{
  int max_dir = 0;
  for ( int i = 1; i < DIM; ++i )
  {
    if ( max[i] - min[i] > max[max_dir] - min[max_dir] ) max_dir = i;
  }
  return max_dir;
}

template<class REAL, int DIM>
inline REAL
BoundingBox_T<REAL,DIM>::SqrSize() const
{
  Real sqrSize = 0.0;
  if (valid())
  {
    for ( int i = 0; i < DIM; ++i )
    {
      sqrSize += (max[i]-min[i])*(max[i]-min[i]);
    }
  }
  return sqrSize;
}

template<class REAL, int DIM>
inline void
BoundingBox_T<REAL,DIM>::pad( const Real & dist )
{
  if (!valid()) return;

  VecType extension;
  for (int i = 0; i < DIM; i++ )
    extension[i] = dist;

  min -= extension;
  max += extension;
}

template<class REAL, int DIM>
inline void
BoundingBox_T<REAL,DIM>::scale( const Real & scale_factor )
{
  if (!valid()) return;

  Real fraction = 0.5 * ( scale_factor - 1. );
  VecType extension;
  for (int i = 0; i < DIM; i++ )
    extension[i] = fraction * ( max[i] - min[i] );
  min -= extension;
  max += extension;
}

typedef BoundingBox_T<float,3> BoundingBox;

} // namespace krino

#endif // Akri_BoundingBox_h
