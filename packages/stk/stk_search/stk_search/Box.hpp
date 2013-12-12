#ifndef STK_SEARCH_BOX_HPP
#define STK_SEARCH_BOX_HPP

#include <stk_search/Point.hpp>
#include <Geom_AxisAlignedBB.h>

namespace stk { namespace search {

template <typename T>
class Box
{
public:
  typedef T value_type;
  typedef Point<value_type> point_type;
  static const int Dim = 3;

  Box( point_type const& x_min_corner = point_type(1,1,1), point_type const& x_max_corner = point_type(0,0,0))
    : m_min_corner(x_min_corner)
    , m_max_corner(x_max_corner)
  {}

  point_type const& min_corner() const { return m_min_corner; }
  point_type      & min_corner()       { return m_min_corner; }
  point_type const& max_corner() const { return m_max_corner; }
  point_type      & max_corner()       { return m_max_corner; }

  value_type get_x_min() const { return m_min_corner[0]; }
  value_type get_y_min() const { return m_min_corner[1]; }
  value_type get_z_min() const { return m_min_corner[2]; }
  value_type get_x_max() const { return m_max_corner[0]; }
  value_type get_y_max() const { return m_max_corner[1]; }
  value_type get_z_max() const { return m_max_corner[2]; }

  void set_min_corner(point_type const& x_min_corner) { m_min_corner = x_min_corner; }
  void set_max_corner(point_type const& x_max_corner) { m_max_corner = x_max_corner; }

  bool operator==(Box<value_type> const& b) const
  { return m_min_corner == b.m_min_corner && m_max_corner == b.m_max_corner; }

  bool operator!=(Box<value_type> const& b) const
  { return !(*this == b); }

  friend std::ostream& operator<<(std::ostream & out, Box<value_type> const& b)
  {
    out << "{" << b.min_corner() << "->" << b.max_corner() << "}";
    return out;
  }

private:
  point_type m_min_corner;
  point_type m_max_corner;
};

}} //namespace stk::search


namespace boost {
namespace geometry {
namespace traits {

// traits for stk::search::Box<T>
template <typename T> struct tag< stk::search::Box<T> > { typedef box_tag type; };
template <typename T> struct point_type< stk::search::Box<T> > { typedef stk::search::Point<T> type; };

template <typename T, size_t Index>
struct indexed_access< stk::search::Box<T>, min_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Box<T> const& s) { return s.min_corner()[Index]; }
};

template <typename T, size_t Index>
struct indexed_access< stk::search::Box<T>, max_corner, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Box<T> const& s) { return s.max_corner()[Index]; }
};


// traits for ::geometry::AxisAlignedBB
template <> struct tag< ::geometry::AxisAlignedBB > { typedef box_tag type; };
template <> struct point_type< ::geometry::AxisAlignedBB > { typedef stk::search::Point<Real> type; };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, min_corner, 0 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_x_min(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, min_corner, 1 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_y_min(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, min_corner, 2 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_z_min(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, max_corner, 0 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_x_max(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, max_corner, 1 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_y_max(); } };

template <>
struct indexed_access< ::geometry::AxisAlignedBB, max_corner, 2 >
{  static inline Real get( ::geometry::AxisAlignedBB const& s) { return s.get_z_max(); } };

}}} // namespace boost::geometry::traits

#endif //STK_SEARCH_BOX_HPP

