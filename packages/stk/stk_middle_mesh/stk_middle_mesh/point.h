#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <iostream>

namespace stk {
namespace middle_mesh {
namespace utils {

struct Point
{
    Point(const double& x_, const double& y_, const double& z_ = 0)
      : x(x_)
      , y(y_)
      , z(z_)
    {}

    Point()
      : x(0)
      , y(0)
      , z(0)
    {}

    double get_x() const { return x; }
    double get_y() const { return y; }
    double get_z() const { return z; }

    double& operator[](int idx) { return *(&(x) + idx); }
    const double& operator[](int idx) const { return *(&(x) + idx); }

    double x;
    double y;
    double z;

    Point operator-() const { return Point(-x, -y, -z); }

    Point& operator+=(const Point& rhs)
    {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;

      return *this;
    }

    Point& operator-=(const Point& rhs) { return operator+=(-rhs); }
};

inline Point operator+(const Point& a, const Point& b)
{
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

inline Point operator-(const Point& a, const Point& b)
{
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Point operator*(const Point& a, const double& b)
{
  return {a.x * b, a.y * b, a.z * b};
}

inline Point operator*(const double& b, const Point& a)
{
  return a * b;
}

inline Point operator/(const Point& a, const double& b)
{
  return {a.x / b, a.y / b, a.z / b};
}

inline std::ostream& operator<<(std::ostream& os, const Point& pt)
{
  os << "(" << pt.get_x() << ", " << pt.get_y() << ", " << pt.get_z() << ")";
  return os;
}

inline double dot(const Point& a, const Point& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline void dot_rev(const Point& a, const Point& b, Point& aBar, Point& bBar, const double dBar)
{
  aBar = b * dBar;
  bBar = a * dBar;
}

inline Point cross(const Point& a, const Point& b)
{
  auto cx = a.y * b.z - a.z * b.y;
  auto cy = -(a.x * b.z - a.z * b.x);
  auto cz = a.x * b.y - a.y * b.x;

  return {cx, cy, cz};
}

// project v onto n
inline Point project(const Point& v, const Point& n)
{
  auto fac = dot(v, n) / dot(n, n);
  return fac * n;
}

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
