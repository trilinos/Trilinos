/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_BoundingBox2_hpp
#define stk_search_BoundingBox2_hpp

#include <stk_util/environment/ReportHandler.hpp>

#include <boost/concept/assert.hpp>
#include <boost/geometry/algorithms/convert.hpp>
#include <boost/geometry/geometries/concepts/point_concept.hpp>

#include <iosfwd>
#include <cmath>
#include <stdint.h>

namespace stk { namespace search {

template <typename T>
class Point3
{
public:
  static const size_t Dim = 3;
  typedef T value_type;

  Point3(value_type x = value_type(), value_type y = value_type(), value_type z = value_type())
    : m_value()
  {
    m_value[0] = x; m_value[1] = y; m_value[2] = z;
  }

  value_type const& operator[](size_t index) const
  {
    ThrowAssert(index < Dim);
    return m_value[index];
  }

  value_type & operator[](size_t index)
  {
    ThrowAssert(index < Dim);
    return m_value[index];
  }

  bool operator==(Point3<value_type> const& p) const
  {
    return  m_value[0] == p.m_value[0]
         && m_value[1] == p.m_value[1]
         && m_value[2] == p.m_value[2];
  }

  bool operator!=(Point3<value_type> const& p) const
  { return !(*this == p); }

  friend std::ostream& operator<<(std::ostream & out, Point3<value_type> const& p)
  {
    out << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return out;
  }

private:
  value_type m_value[Dim];
};

template <typename T>
class Sphere
{
public:
  typedef T value_type;
  typedef Point3<value_type> point_type;

  Sphere( point_type const& x_center = point_type(), value_type const& x_radius = static_cast<value_type>(-1))
    : m_center(x_center)
    , m_radius(x_radius)
  {}

  point_type const& center() const { return m_center; }
  value_type const& radius() const { return m_radius; }

  void set_center(point_type const& c) { m_center = c; }
  void set_radius(value_type const& r) { m_radius = r; }

  bool operator==(Sphere<value_type> const& s) const
  { return m_radius == s.m_radius && m_center == s.m_center; }

  bool operator!=(Sphere<value_type> const& s) const
  { return !(*this == s); }

  bool is_valid() const
  { return static_cast<value_type>(0) <= m_radius; }

  friend std::ostream& operator<<(std::ostream & out, Sphere<value_type> const& s)
  {
    out << "{" << s.center() << ":" << s.radius() << "}";
    return out;
  }

private:
  point_type m_center;
  value_type m_radius;
};


template <typename T>
class AABox
{
public:
  typedef T value_type;
  typedef Point3<value_type> point_type;

  AABox( point_type const& x_min_corner = point_type(1,1,1), point_type const& x_max_corner = point_type(0,0,0))
    : m_min_corner(x_min_corner)
    , m_max_corner(x_max_corner)
  {}

  point_type const& min_corner() const { return m_min_corner; }
  point_type const& max_corner() const { return m_max_corner; }

  void set_min_corner(point_type const& x_min_corner) { m_min_corner = x_min_corner; }
  void set_max_corner(point_type const& x_max_corner) { m_max_corner = x_max_corner; }

  bool operator==(AABox<value_type> const& b) const
  { return m_min_corner == b.m_min_corner && m_max_corner == b.m_max_corner; }

  bool operator!=(AABox<value_type> const& b) const
  { return !(*this == b); }


  bool is_valid() const
  {
    return    m_min_corner[0] <= m_max_corner[0]
           && m_min_corner[1] <= m_max_corner[1]
           && m_min_corner[2] <= m_max_corner[2];
  }

  friend std::ostream& operator<<(std::ostream & out, AABox<value_type> const& b)
  {
    out << "{" << b.min_corner() << "->" << b.max_corner() << "}";
    return out;
  }

private:
  point_type m_min_corner;
  point_type m_max_corner;
};


template <typename T>
inline Point3<T> const& min_corner(Point3<T> const& p)
{
  return p;
}

template <typename T>
inline Point3<T> const& max_corner(Point3<T> const& p)
{
  return p;
}

template <typename T>
inline Point3<T> const& center(Point3<T> const& p)
{
  return p;
}

template <typename T>
inline Point3<T> min_corner(Sphere<T> const& s)
{
  Point3<T> p(s.center());
  p[0] -= s.radius();
  p[1] -= s.radius();
  p[2] -= s.radius();
  return p;
}

template <typename T>
inline Point3<T> max_corner(Sphere<T> const& s)
{
  Point3<T> p(s.center());
  p[0] += s.radius();
  p[1] += s.radius();
  p[2] += s.radius();
  return p;
}

template <typename T>
inline Point3<T> const& center(Sphere<T> const& s)
{
  return s.center();
}

template <typename T>
inline Point3<T> const& min_corner(AABox<T> const& b)
{
  return b.min_corner();
}

template <typename T>
inline Point3<T> const& max_corner(AABox<T> const& b)
{
  return b.max_corner();
}

template <typename T>
inline Point3<T> center(AABox<T> const& b)
{
  Point3<T> p;
  p[0] = (b.min_corner()[0] + b.max_corner()[0])/2;
  p[1] = (b.min_corner()[1] + b.max_corner()[1])/2;
  p[2] = (b.min_corner()[2] + b.max_corner()[2])/2;
  return p;
}


template <typename T>
inline bool intersects(Point3<T> const& a, Point3<T> const& b)
{
  return (a==b);
}


template <typename T>
inline bool intersects(Point3<T> const& a, Sphere<T> const& b)
{
  T dst2 = (a[0]-b.center()[0])*(a[0]-b.center()[0]);
  dst2 += (a[1]-b.center()[1])*(a[1]-b.center()[1]);
  dst2 += (a[2]-b.center()[2])*(a[2]-b.center()[2]);

  return (dst2 <= b.radius()*b.radius());
}

template <typename T>
inline bool intersects(Sphere<T> const& a, Point3<T> const& b)
{
  return intersects(b,a);
}

template <typename T>
inline bool intersects(Point3<T> const& a, AABox<T> const& b)
{
  return b.min_corner()[0] <= a[0] && a[0] <= b.max_corner()[0]
      && b.min_corner()[1] <= a[1] && a[1] <= b.max_corner()[1]
      && b.min_corner()[2] <= a[2] && a[2] <= b.max_corner()[2];
}

template <typename T>
inline bool intersects(AABox<T> const& a, Point3<T> const& b)
{
  return intersects(b,a);
}

}} //namespace stk::search

namespace boost { namespace geometry { namespace traits {

// traits for stk::search::Point3<T>
template <typename T> struct tag< stk::search::Point3<T> > { typedef point_tag type; };
template <typename T> struct coordinate_type< stk::search::Point3<T> > { typedef T type; };
template <typename T> struct coordinate_system< stk::search::Point3<T> > { typedef cs::cartesian type; };
template <typename T> struct dimension< stk::search::Point3<T> > : public boost::mpl::int_<3> {};

template <typename T, int Index>
struct access< stk::search::Point3<T>, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Point3<T> const& p) { return p[Index]; }
  static inline void set( stk::search::Point3<T> const& p, T const& v) { p[Index] = v; }
};

// traits for stk::search::Sphere<T>
template <typename T> struct tag< stk::search::Sphere<T> > { typedef box_tag type; };
template <typename T> struct point_type< stk::search::Sphere<T> > { typedef stk::search::Point3<T> type; };

template <typename T, int Index>
struct indexed_access< stk::search::Sphere<T>, min_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Sphere<T> const& s) { return s.center()[Index] - s.radius(); }
};

template <typename T, int Index>
struct indexed_access< stk::search::Sphere<T>, max_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Sphere<T> const& s) { return s.center()[Index] + s.radius(); }
};


// traits for stk::search::AABox<T>
template <typename T> struct tag< stk::search::AABox<T> > { typedef box_tag type; };
template <typename T> struct point_type< stk::search::AABox<T> > { typedef stk::search::Point3<T> type; };

template <typename T, int Index>
struct indexed_access< stk::search::AABox<T>, min_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::AABox<T> const& s) { return s.min_corner()[Index]; }
};

template <typename T, int Index>
struct indexed_access< stk::search::AABox<T>, max_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::AABox<T> const& s) { return s.max_corner()[Index]; }
};

}}} // namespace boost::geometry::traits

#endif //stk_search_BoundingBox2_hpp
