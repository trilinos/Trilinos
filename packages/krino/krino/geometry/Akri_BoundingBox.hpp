// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_BoundingBox_h
#define Akri_BoundingBox_h

#include <stk_math/StkVector.hpp>

#include <limits>
#include <vector>

namespace krino {

template<class REAL, unsigned DIM>
class BoundingBox_T {

public:
  typedef REAL Real;
  typedef stk::math::Vec<Real,3> VecType;

private:
  VecType min;
  VecType max;

public:
  static const BoundingBox_T<REAL,DIM> ENTIRE_DOMAIN;

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

  void pad( const Real & dist );
  void pad_epsilon();

  template<class VECTYPE>
  void shift( const VECTYPE & shiftVec );

  void scale( const Real & fraction );

  Real size_squared() const;

  bool valid() const
  {
    bool is_not_valid = (min[0]>max[0]);
    return !is_not_valid;
  }
  
  bool intersects( const BoundingBox_T<REAL,DIM> & bbox ) const;
  bool contains( const BoundingBox_T<REAL,DIM> & bbox ) const;
  
  unsigned max_span_direction() const;

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

template<class REAL, unsigned DIM>
template<class VECTYPE>
inline void
BoundingBox_T<REAL,DIM>::accommodate( const VECTYPE & pt )
{
  for ( unsigned i = 0; i < DIM; ++i )
  {
    if ( pt[i] < min[i] ) min[i] = pt[i];
    if ( pt[i] > max[i] ) max[i] = pt[i];
  }
}

template<class REAL, unsigned DIM>
template<class REALARG>
inline void
BoundingBox_T<REAL,DIM>::accommodate( const BoundingBox_T<REALARG,DIM> & bbox )
{
  for ( unsigned i = 0; i < DIM; ++i )
  {
    if ( bbox.get_min()[i] < min[i] ) min[i] = bbox.get_min()[i];
    if ( bbox.get_max()[i] > max[i] ) max[i] = bbox.get_max()[i];
  }
}

template<class REAL, unsigned DIM>
template<class VECTYPE>
inline bool
BoundingBox_T<REAL,DIM>::contains( const VECTYPE & pt ) const
{
  for ( unsigned i = 0; i < DIM; ++i )
  {
    if ( pt[i] < min[i] || pt[i] > max[i] ) return false;
  }
  return true;
}

template<class REAL, unsigned DIM>
inline bool
BoundingBox_T<REAL,DIM>::intersects( const BoundingBox_T<REAL,DIM> & bbox ) const
{
  for ( unsigned i = 0; i < DIM; ++i )
  {
    if ( bbox.max[i] < min[i] || bbox.min[i] > max[i] ) return false;
  }
  return true;
}

template<class REAL, unsigned DIM>
inline bool
BoundingBox_T<REAL,DIM>::contains( const BoundingBox_T<REAL,DIM> & bbox ) const
{
  for ( unsigned i = 0; i < DIM; ++i )
  {
    if ( bbox.min[i] < min[i] || bbox.max[i] > max[i] ) return false;
  }
  return true;
}

template<class REAL, unsigned DIM>
inline unsigned
BoundingBox_T<REAL,DIM>::max_span_direction() const
{
  unsigned max_dir = 0;
  for ( unsigned i = 1; i < DIM; ++i )
  {
    if ( max[i] - min[i] > max[max_dir] - min[max_dir] ) max_dir = i;
  }
  return max_dir;
}

template<class REAL, unsigned DIM>
inline REAL
BoundingBox_T<REAL,DIM>::size_squared() const
{
  Real sqrSize = 0.0;
  if (valid())
  {
    for ( unsigned i = 0; i < DIM; ++i )
    {
      sqrSize += (max[i]-min[i])*(max[i]-min[i]);
    }
  }
  return sqrSize;
}

template<class REAL, unsigned DIM>
inline void
BoundingBox_T<REAL,DIM>::pad( const Real & dist )
{
  if (!valid()) return;

  VecType extension;
  for (unsigned i = 0; i < DIM; i++ )
    extension[i] = dist;

  min -= extension;
  max += extension;
}

template<class REAL, unsigned DIM>
inline void
BoundingBox_T<REAL,DIM>::scale( const Real & scale_factor )
{
  if (!valid()) return;

  Real fraction = 0.5 * ( scale_factor - 1. );
  VecType extension;
  for (unsigned i = 0; i < DIM; i++ )
    extension[i] = fraction * ( max[i] - min[i] );
  min -= extension;
  max += extension;
}

template<class REAL, unsigned DIM>
inline void
BoundingBox_T<REAL,DIM>::pad_epsilon()
{
  if (!valid()) return;
  pad(std::sqrt(size_squared())*std::numeric_limits<Real>::epsilon());
}

template<class REAL, unsigned DIM>
template<class VECTYPE>
inline void
BoundingBox_T<REAL,DIM>::shift( const VECTYPE & shiftVec )
{
  for ( unsigned i = 0; i < DIM; ++i )
  {
    min[i] += shiftVec[i];
    max[i] += shiftVec[i];
  }
}

template<class REAL, unsigned DIM>
const BoundingBox_T<REAL,DIM> BoundingBox_T<REAL,DIM>::ENTIRE_DOMAIN(VecType(-std::numeric_limits<REAL>::max(), -std::numeric_limits<REAL>::max(), -std::numeric_limits<REAL>::max()),
      VecType(std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max(), std::numeric_limits<REAL>::max()));

typedef BoundingBox_T<float,3> BoundingBox;

} // namespace krino

#endif // Akri_BoundingBox_h
