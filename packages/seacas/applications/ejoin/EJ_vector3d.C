// Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "EJ_vector3d.h"
#include <cmath>

//----------------------------------------------------------------------------
vector3d::vector3d() = default;

//----------------------------------------------------------------------------
vector3d::vector3d(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

//----------------------------------------------------------------------------
vector3d::vector3d(double location[3]) : x(location[0]), y(location[1]), z(location[2]) {}

//----------------------------------------------------------------------------
vector3d::vector3d(const vector3d &from) = default;

void vector3d::set(double X, double Y, double Z)
{
  x = X;
  y = Y;
  z = Z;
}

void vector3d::set(const double location[3])
{
  x = location[0];
  y = location[1];
  z = location[2];
}

vector3d &vector3d::operator=(const vector3d &from) = default;

vector3d &vector3d::reverse()
{
  x = -x;
  y = -y;
  z = -z;
  return *this;
}

bool vector3d::operator==(const vector3d &from) const
{
  return (x == from.x && y == from.y && z == from.z);
}

bool vector3d::operator!=(const vector3d &from) const
{
  return (x != from.x || y != from.y || z != from.z);
}

vector3d operator+(const vector3d &lhs, const vector3d &rhs)
{
  vector3d tmp(lhs);
  return tmp += rhs;
}

vector3d operator-(const vector3d &lhs, const vector3d &rhs)
{
  vector3d tmp(lhs);
  return tmp -= rhs;
}

vector3d operator*(const vector3d &lhs, double scalar)
{
  vector3d tmp(lhs);
  return tmp *= scalar;
}

vector3d vector3d::operator-() const
{
  vector3d tmp(x, y, z);
  return tmp *= -1.0;
}

vector3d operator*(double scalar, const vector3d &from)
{
  vector3d tmp(from);
  return tmp *= scalar;
}

vector3d operator/(const vector3d &lhs, double scalar)
{
  if (scalar != 0.0) {
    vector3d tmp(lhs);
    return tmp /= scalar;
  }

  return {HUGE_VAL, HUGE_VAL, HUGE_VAL};
}

vector3d &vector3d::operator/=(double scalar)
{
  if (scalar != 0.0) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
  }
  else {
    x = HUGE_VAL;
    y = HUGE_VAL;
    z = HUGE_VAL;
  }

  return *this;
}

double vector3d::length() const { return sqrt(x * x + y * y + z * z); }
