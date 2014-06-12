/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_BoundingBox_hpp
#define stk_search_BoundingBox_hpp

#include <iosfwd>
#include <stk_search/BoundingBoxCompare.hpp>
#include <cmath>
#include <stdint.h>

namespace stk_classic {
namespace search {
namespace box {

template <class K, class T, int Dim>
  struct PointBoundingBox;

template <class K, class T, int Dim>
  struct SphereBoundingBox;

template <class K, class T, int Dim>
  struct AxisAlignedBoundingBox;

template <class K = uint64_t, class T = float, int Dim = 3>
struct PointBoundingBox
{
  typedef K Key;
  typedef T Data;

  const static int DIMENSION = Dim;

  inline Data lower(int axis) const { return center[axis]; }
  inline Data middle(int axis) const { return center[axis]; }
  inline Data upper(int axis) const { return center[axis]; }

  inline Data length(int /* axis */ = 0) const { return 0; }

  PointBoundingBox()
    : center(), key()
  {}

  PointBoundingBox(const Data center_[], const Key & key_)
    : center(), key(key_)
  {
    set_center(center_);
  }

  PointBoundingBox(const PointBoundingBox &point_)
    : center(), key(point_.key)
  {
    set_center(point_.center);
  }

  PointBoundingBox & operator = ( const PointBoundingBox &point_) {
    if (this != &point_) {
      set_center(point_.center);
      key = point_.key;
    }
    return *this;
  }

  inline
  void set_center(const Data center_[]) {
    for (int i = 0; i<DIMENSION; ++i)
      center[i] = center_[i];
  }

  bool intersect( const PointBoundingBox<K,T,Dim> & point_) const {
    Data dist = 0;
    for (int i=0; i<DIMENSION; ++i) {
      dist += std::abs(middle(i) - point_.middle(i));
    }
    return dist == 0;
  }

  bool intersect( const SphereBoundingBox<K,T,Dim> & sphere_) const {
    Data tmp = 0;
    Data dist = 0;
    for (int i=0; i<DIMENSION; ++i) {
      tmp = sphere_.middle(i) - middle(i);
      dist += tmp * tmp;
    }
    return dist < sphere_.radius * sphere_.radius;
  }

  bool intersect( const AxisAlignedBoundingBox<K,T,Dim> & box_) const {
    for (int i=0; i<DIMENSION; ++i) {
      if( middle(i) < box_.lower(i) || middle(i) > box_.upper(i)) {
        return false;
      }
    }
    return true;
  }

  Data center[DIMENSION];
  Key  key;

};

template <class K = uint64_t, class T = float, int Dim = 3>
struct SphereBoundingBox
{
  typedef K Key;
  typedef T Data;

  const static int DIMENSION = Dim;

  inline Data lower(int axis) const { return center[axis] - radius; }
  inline Data middle(int axis) const { return center[axis]; }
  inline Data upper(int axis) const { return center[axis] + radius; }

  inline Data length(int /* axis */ = 0) const { return 2*radius; }

  inline void expand(const Data delta) { radius += delta; }


  SphereBoundingBox()
    : center(), radius(0), key()
  {}

  SphereBoundingBox(const SphereBoundingBox &sphere_)
    : center(), radius(sphere_.radius), key(sphere_.key)
  {
    set_center(sphere_.center);
  }

  SphereBoundingBox & operator = ( const SphereBoundingBox &sphere_) {
    if (this != &sphere_) {
      set_center(sphere_.center);
      radius = sphere_.radius;
      key = sphere_.key;
    }
    return *this;
  }

  SphereBoundingBox( const PointBoundingBox<K,T,Dim> & point_ )
    : center(), radius(0), key(point_.key)
  {
    set_center(point_.center);
  }

  SphereBoundingBox( const AxisAlignedBoundingBox<K,T,Dim> & box_)
    : center(), radius(0), key(box_.key)
  {
    Data dist;
    Data radius_squared = 0;
    for (int i=0; i<Dim; ++i) {
      center[i] = box_.middle(i);
      dist = box_.upper(i) - center[i];
      radius_squared += dist * dist;
    }
    radius = std::sqrt(radius_squared);
  }

  SphereBoundingBox(const Data center_[], const Data radius_, const Key & key_)
    : center(), radius(radius_), key(key_)
  {
    set_center(center_);
  }

  inline
  void set_center(const Data center_[]) {
    for (int i = 0; i<DIMENSION; ++i)
      center[i] = center_[i];
  }

  inline
  void set_radius(const Data radius_) {
    radius = radius_;
  }

  inline
  bool intersect( const PointBoundingBox<K,T,Dim> & point_) const {
    return point_.intersect(*this);
  }

  bool intersect( const SphereBoundingBox<K,T,Dim> & sphere_) const {
    Data tmp = 0;
    Data dist = 0;
    for (int i=0; i<DIMENSION; ++i) {
      tmp = middle(i) - sphere_.middle(i);
      dist += tmp * tmp;
    }
    Data radius_sum = radius + sphere_.radius;
    return dist < radius_sum * radius_sum;
  }

  bool intersect( const AxisAlignedBoundingBox<K,T,Dim> & box_) const {
    Data tmp = 0;
    Data dist = 0;
    for (int i=0; i<DIMENSION; ++i) {
      if (middle(i) < box_.lower(i)) {
        tmp = middle(i) - box_.lower(i);
        dist += tmp*tmp;
      }
      else if (middle(i) > box_.upper(i)) {
        tmp = middle(i) - box_.upper(i);
        dist += tmp*tmp;
      }
    }

    return dist <= radius*radius;;
  }
  Data center[DIMENSION];
  Data radius;
  Key  key;

};


template <class K = uint64_t, class T = float, int Dim = 3>
struct AxisAlignedBoundingBox
{
  typedef K Key;
  typedef T Data;

  const static int DIMENSION = Dim;

  inline Data lower(int axis)  const { return box[axis]; }
  inline Data middle(int axis) const { return lower(axis) + length(axis)/2; }
  inline Data upper(int axis)  const { return box[axis+DIMENSION]; }

  inline Data length(int axis) const { return upper(axis) - lower(axis); }

  inline void expand(const Data delta) {
    for (int i=0; i<DIMENSION; ++i) {
      box[i] -= delta;
      box[i+DIMENSION] += delta;
    }
  }

  AxisAlignedBoundingBox()
    : box(), key()
  {}

  AxisAlignedBoundingBox(const Data box_[], const Key & key_)
    : box(), key(key_)
  {
    set_box(box_);
  }

  AxisAlignedBoundingBox(const AxisAlignedBoundingBox &box_)
    : box(), key(box_.key)
  {
    set_box(box_.box);
  }

  AxisAlignedBoundingBox & operator = ( const AxisAlignedBoundingBox &box_) {
    if (this != &box_) {
      set_box(box_.box);
      key = box_.key;
    }
    return *this;
  }

  AxisAlignedBoundingBox( const PointBoundingBox<K,T,Dim> & point_)
    : box(), key(point_.key)
  {
    for (int i=0; i<Dim; ++i) {
      box[i] = point_.lower(i);
      box[i+Dim] = point_.upper(i);
    }
  }

  AxisAlignedBoundingBox( const SphereBoundingBox<K,T,Dim> & sphere_)
    : box(), key(sphere_.key)
  {
    for (int i=0; i<Dim; ++i) {
      box[i] = sphere_.lower(i);
      box[i+Dim] = sphere_.upper(i);
    }
  }

  inline
  void set_box(const Data box_[]) {
    for (int i = 0; i<2*DIMENSION; ++i)
      box[i] = box_[i];
  }

  inline
  bool intersect( const PointBoundingBox<K,T,Dim> & point_) const {
    return point_.intersect(*this);
  }

  inline
  bool intersect( const SphereBoundingBox<K,T,Dim> & sphere_) const {
    return sphere_.intersect(*this);
  }

  bool intersect( const AxisAlignedBoundingBox<K,T,Dim> & box_) const {
    for (int i=0; i<DIMENSION; ++i) {
      if( upper(i) < box_.lower(i) || lower(i) > box_.upper(i)) return false;
    }
    return true;
  }

  Data box[2*DIMENSION];
  Key  key;

};



} // namespace box
} // namespace search
} // namespace stk_classic

#endif // stk_search_BoundingBox_hpp
