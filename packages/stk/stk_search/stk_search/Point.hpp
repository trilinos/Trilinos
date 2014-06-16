#ifndef STK_SEARCH_POINT_HPP
#define STK_SEARCH_POINT_HPP

#include <stk_util/environment/ReportHandler.hpp>
#include <boost/geometry/geometry.hpp>
#include <iosfwd>

namespace stk { namespace search {

template <typename T>
class Point
{
public:
  static const unsigned Dim = 3;
  typedef T value_type;

  Point(value_type x = value_type(), value_type y = value_type(), value_type z = value_type())
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

  bool operator==(Point<value_type> const& p) const
  {
    return  m_value[0] == p.m_value[0]
         && m_value[1] == p.m_value[1]
         && m_value[2] == p.m_value[2];
  }

  bool operator!=(Point<value_type> const& p) const
  { return !(*this == p); }

  value_type get_x_min() const { return m_value[0]; }
  value_type get_y_min() const { return m_value[1]; }
  value_type get_z_min() const { return m_value[2]; }
  value_type get_x_max() const { return m_value[0]; }
  value_type get_y_max() const { return m_value[1]; }
  value_type get_z_max() const { return m_value[2]; }


  friend std::ostream& operator<<(std::ostream & out, Point<value_type> const& p)
  {
    out << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return out;
  }

private:
  value_type m_value[Dim];
};

}} // stk::search

namespace boost { namespace geometry { namespace traits {

// traits for stk::search::Point<T>
template <typename T> struct tag< stk::search::Point<T> > { typedef point_tag type; };
template <typename T> struct coordinate_type< stk::search::Point<T> > { typedef T type; };
template <typename T> struct coordinate_system< stk::search::Point<T> > { typedef cs::cartesian type; };
template <typename T> struct dimension< stk::search::Point<T> > : public boost::mpl::int_<3> {};

template <typename T, size_t Index>
struct access< stk::search::Point<T>, Index >
{
  BOOST_STATIC_ASSERT((Index < 3));
  static inline T const& get( stk::search::Point<T> const& p) { return p[Index]; }
  static inline void set( stk::search::Point<T> const& p, T const& v) { p[Index] = v; }
};

}}} // namespace boost::geometry::traits

#endif //STK_SEARCH_POINT_HPP
