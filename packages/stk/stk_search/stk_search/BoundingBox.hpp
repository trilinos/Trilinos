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

#include <boost/concept/assert.hpp>
#include <boost/geometry/algorithms/convert.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>

#include <boost/geometry/geometries/adapted/c_array.hpp>

BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(boost::geometry::cs::cartesian)

namespace stk {
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

  typedef PointBoundingBox<K,T,Dim> self_type;

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

  template <class KT, class TT>
  bool intersect( const PointBoundingBox<KT,TT,Dim> & point_) const {
    Data dist = 0;
    for (int i=0; i<DIMENSION; ++i) {
      dist += std::abs(middle(i) - point_.middle(i));
    }
    return dist == 0;
  }

  template <class KT, class TT>
  bool intersect( const SphereBoundingBox<KT,TT,Dim> & sphere_) const {
    Data tmp = 0;
    Data dist = 0;
    for (int i=0; i<DIMENSION; ++i) {
      tmp = sphere_.middle(i) - middle(i);
      dist += tmp * tmp;
    }
    return dist < sphere_.radius * sphere_.radius;
  }

  template <class KT, class TT>
  bool intersect( const AxisAlignedBoundingBox<KT,TT,Dim> & box_) const {
    for (int i=0; i<DIMENSION; ++i) {
      if( middle(i) < box_.lower(i) || middle(i) > box_.upper(i)) {
        return false;
      }
    }
    return true;
  }

  bool operator==(const self_type & rhs) const
  {
    bool result = key == rhs.key;
    for (int i=0; result && i<Dim; ++i) {
      result = center[i] == rhs.center[i];
    }
    return result;
  }

  bool operator!=(const self_type & rhs) const
  { return !(*this == rhs); }

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
  inline void scale (const Data delta) { radius *= delta; }

  typedef SphereBoundingBox<K,T,Dim> self_type;

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
    : center(), radius(std::abs(radius_)), key(key_)
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

  template <class KT, class TT>
  bool intersect( const PointBoundingBox<KT,TT,Dim> & point_) const {
    return point_.intersect(*this);
  }

  template <class KT, class TT>
  bool intersect( const SphereBoundingBox<KT,TT,Dim> & sphere_) const {
    Data tmp = 0;
    Data dist = 0;
    for (int i=0; i<DIMENSION; ++i) {
      tmp = middle(i) - static_cast<Data>(sphere_.middle(i));
      dist += tmp * tmp;
    }
    Data radius_sum = radius + static_cast<Data>(sphere_.radius);
    return dist < radius_sum * radius_sum;
  }

  template <class KT, class TT>
  bool intersect( const AxisAlignedBoundingBox<KT,TT,Dim> & box_) const {
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

  bool operator==(const self_type & rhs) const
  {
    bool result = key == rhs.key;
    result = result && (radius == rhs.radius);
    for (int i=0; result && i<Dim; ++i) {
      result = center[i] == rhs.center[i];
    }
    return result;
  }

  bool operator!=(const self_type & rhs) const
  { return !(*this == rhs); }


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

  typedef AxisAlignedBoundingBox<K,T,Dim> self_type;

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
  inline void scale (const Data delta) {
    for (int i=0; i<DIMENSION; ++i) {
      const Data f = (delta-1)/2; // Scale like r *= delta
      const Data d = f*(box[i+DIMENSION] - box[i]);
      box[i]           -= d;
      box[i+DIMENSION] += d;
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

  template <class KT, class TT>
  bool intersect( const PointBoundingBox<KT,TT,Dim> & point_) const {
    return point_.intersect(*this);
  }

  template <class KT, class TT>
  bool intersect( const SphereBoundingBox<KT,TT,Dim> & sphere_) const {
    return sphere_.intersect(*this);
  }

  template <class KT, class TT>
  bool intersect( const AxisAlignedBoundingBox<KT,TT,Dim> & box_) const {
    for (int i=0; i<DIMENSION; ++i) {
      if( upper(i) < box_.lower(i) || lower(i) > box_.upper(i)) return false;
    }
    return true;
  }

  bool operator==(const self_type & rhs) const
  {
    bool result = key == rhs.key;
    for (int i=0; result && i<2*Dim; ++i) {
      result = box[i] == rhs.box[i];
    }
    return result;
  }

  bool operator!=(const self_type & rhs) const
  { return !(*this == rhs); }

  Data box[2*DIMENSION];
  Key  key;
};

} // namespace box
} // namespace search
} // namespace stk

namespace boost { namespace geometry { namespace traits {

// PointBoundingBox
template <class K, class T, int Dim>
struct tag< stk::search::box::PointBoundingBox<K,T,Dim> >
{
    typedef point_tag type;
};

template <class K, class T, int Dim>
struct coordinate_type< stk::search::box::PointBoundingBox<K,T,Dim> >
{
    typedef T type;
};

template <class K, class T, int Dim>
struct coordinate_system< stk::search::box::PointBoundingBox<K,T,Dim> >
{
    typedef boost::geometry::cs::cartesian type;
};

template <class K, class T, int Dim>
struct dimension< stk::search::box::PointBoundingBox<K,T,Dim> > : public boost::mpl::int_<Dim> {};

template <class K, class T, int Dim, size_t IndexDimension>
struct access< stk::search::box::PointBoundingBox<K,T,Dim>, IndexDimension>
{
    typedef T coordinate_type;
    typedef stk::search::box::PointBoundingBox<K,T,Dim> point_type;

    static inline coordinate_type get(point_type const& b)
    {
        return b.center[IndexDimension];
    }

    static inline void set(point_type & b, coordinate_type const& value)
    {
        b.center[IndexDimension] = value;
    }
};

template <class K, class T, int Dim, size_t IndexDimension>
struct indexed_access< stk::search::box::PointBoundingBox<K,T,Dim>, max_corner, IndexDimension>
{
    typedef T coordinate_type;
    typedef stk::search::box::PointBoundingBox<K,T,Dim> box_type;

    static inline coordinate_type get(box_type const& b)
    {
        return b.center[IndexDimension];
    }

    static inline void set(box_type & b, coordinate_type const& value)
    {
        b.center[IndexDimension] = value;
    }
};


// SphereBoundingBox
template <class K, class T, int Dim>
struct tag< stk::search::box::SphereBoundingBox<K,T,Dim> >
{
    typedef box_tag type;
};

template <class K, class T, int Dim>
struct point_type< stk::search::box::SphereBoundingBox<K,T,Dim> >
{
    typedef T type[Dim];
};

template <class K, class T, int Dim, size_t IndexDimension>
struct indexed_access< stk::search::box::SphereBoundingBox<K,T,Dim>, min_corner, IndexDimension>
{
    typedef T coordinate_type;
    typedef stk::search::box::SphereBoundingBox<K,T,Dim> box_type;

    static inline coordinate_type get(box_type const& b)
    {
        return b.center[IndexDimension] - b.radius;
    }
};

template <class K, class T, int Dim, size_t IndexDimension>
struct indexed_access< stk::search::box::SphereBoundingBox<K,T,Dim>, max_corner, IndexDimension>
{
    typedef T coordinate_type;
    typedef stk::search::box::SphereBoundingBox<K,T,Dim> box_type;

    static inline coordinate_type get(box_type const& b)
    {
        return b.center[IndexDimension] + b.radius;
    }
};

// AxisAlignedBoundingBox
template <class K, class T, int Dim>
struct tag< stk::search::box::AxisAlignedBoundingBox<K,T,Dim> >
{
    typedef box_tag type;
};

template <class K, class T, int Dim>
struct point_type< stk::search::box::AxisAlignedBoundingBox<K,T,Dim> >
{
    typedef T type[Dim];
};

template <class K, class T, int Dim, size_t IndexDimension>
struct indexed_access< stk::search::box::AxisAlignedBoundingBox<K,T,Dim>, min_corner, IndexDimension>
{
    typedef T coordinate_type;
    typedef stk::search::box::AxisAlignedBoundingBox<K,T,Dim> box_type;

    static inline coordinate_type get(box_type const& b)
    {
        return b.box[IndexDimension];
    }

    static inline void set(box_type & b, coordinate_type const& value)
    {
        b.box[IndexDimension] = value;
    }
};

template <class K, class T, int Dim, size_t IndexDimension>
struct indexed_access< stk::search::box::AxisAlignedBoundingBox<K,T,Dim>, max_corner, IndexDimension>
{
    typedef T coordinate_type;
    typedef stk::search::box::AxisAlignedBoundingBox<K,T,Dim> box_type;

    static inline coordinate_type get(box_type const& b)
    {
        return b.box[Dim + IndexDimension];
    }

    static inline void set(box_type & b, coordinate_type const& value)
    {
        b.box[Dim + IndexDimension] = value;
    }
};

}}} // namespace boost::geometry::traits

#endif // stk_search_BoundingBox_hpp
