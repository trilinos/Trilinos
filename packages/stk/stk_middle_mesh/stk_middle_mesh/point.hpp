#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <complex>
#include <iostream>
#include <array>
#include "is_point.hpp"
#include "complex_utils.hpp"


namespace stk {
namespace middle_mesh {
namespace utils {

template <typename T>
struct PointT
{
    using value_type = T;

    constexpr PointT(const T& x_, const T& y_, const T& z_ = 0)
      : x(x_)
      , y(y_)
      , z(z_)
    {}

    constexpr PointT()
      : x(T())
      , y(T())
      , z(T())
    {}

    template <typename T2>
    PointT(const PointT<T2>& other)
      : x(other.x)
      , y(other.y)
      , z(other.z)
    {}

    template <typename T2>
    PointT& operator=(const PointT<T2>& other)
    {
      x = other.x;
      y = other.y;
      z = other.z;

      return *this;
    }

    constexpr T get_x() const { return x; }
    constexpr T get_y() const { return y; }
    constexpr T get_z() const { return z; }

    constexpr T& operator[](int idx)
    { 
      std::array<T*, 3> vals = {&x, &y, &z};
      return *(vals[idx]);
    }
    constexpr const T& operator[](int idx) const
    { 
      std::array<const T*, 3> vals = {&x, &y, &z};
      return *(vals[idx]);
    }

    T x;
    T y;
    T z;

    constexpr PointT operator-() const { return PointT(-x, -y, -z); }

    constexpr PointT& operator+=(const PointT& rhs)
    {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;

      return *this;
    }

    constexpr PointT& operator-=(const PointT& rhs) { return operator+=(-rhs); }
};

using Point = PointT<double>;

template <typename T, typename T2>
constexpr PointT<std::common_type_t<T, T2>> operator+(const PointT<T>& a, const PointT<T2>& b)
{
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

template <typename T, typename T2>
constexpr PointT<std::common_type_t<T, T2>> operator-(const PointT<T>& a, const PointT<T2>& b)
{
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}


template <typename T, typename T2, typename T3 = T2, std::enable_if_t<!IsPoint<T3>::Value, bool> = true>
constexpr PointT<std::common_type_t<T, T2>> operator*(const PointT<T>& a, const T2& b)
{
  return {a.x * b, a.y * b, a.z * b};
}

// Idea: test for T2::value_type
template <typename T, typename T2, typename T3 = T2, std::enable_if_t<!IsPoint<T3>::Value, bool> = true>
constexpr PointT<std::common_type_t<T, T2>> operator*(const T2& b, const PointT<T>& a)
{
  return a * b;
}

template <typename T, typename T2>
constexpr PointT<std::common_type_t<T, T2>> operator/(const PointT<T>& a, const T2& b)
{
  return {a.x / b, a.y / b, a.z / b};
}

template <typename T, typename T2>
constexpr bool operator==(const PointT<T>& lhs, const PointT<T2>& rhs)
{
  return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
}

template <typename T, typename T2>
constexpr bool operator!=(const PointT<T>& lhs, const PointT<T2>& rhs)
{
  return !(lhs == rhs);
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const PointT<T>& pt)
{
  os << "(" << pt.get_x() << ", " << pt.get_y() << ", " << pt.get_z() << ")";
  return os;
}

// Note that this is the simple dot product, not the conjugate transpose product
// when T is complex
template <typename T, typename T2>
constexpr std::common_type_t<T, T2> dot(const PointT<T>& a, const PointT<T2>& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

template <typename T, typename T2>
constexpr void dot_rev(const PointT<T>& a, const PointT<T2>& b, PointT<std::common_type_t<T, T2>>& aBar,
                       PointT<std::common_type_t<T, T2>>& bBar, const std::common_type_t<T, T2> dBar)
{
  aBar = b * dBar;
  bBar = a * dBar;
}

template <typename T, typename T2>
constexpr PointT<std::common_type_t<T, T2>> cross(const PointT<T>& a, const PointT<T2>& b)
{
  auto cx = a.y * b.z - a.z * b.y;
  auto cy = -(a.x * b.z - a.z * b.x);
  auto cz = a.x * b.y - a.y * b.x;

  return {cx, cy, cz};
}

// project v onto n
template <typename T, typename T2>
constexpr PointT<std::common_type_t<T, T2>> project(const PointT<T>& v, const PointT<T2>& n)
{
  auto fac = dot(v, n) / dot(n, n);
  return fac * n;
}

} // namespace utils
} // namespace middle_mesh
} // namespace stk

#endif
