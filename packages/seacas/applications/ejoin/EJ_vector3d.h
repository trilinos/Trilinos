#ifndef VECTOR3D_H
#define VECTOR3D_H

// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

class vector3d
{
public:
  // construction
  vector3d();
  vector3d(double X, double Y, double Z);
  explicit vector3d(double location[3]);
  vector3d(const vector3d &from);

  double x{}, y{}, z{};

  vector3d &operator=(const vector3d &from);
  bool      operator==(const vector3d &from) const;
  bool      operator!=(const vector3d &from) const;
  void      set(double X, double Y, double Z);
  void      set(const double location[3]);
  vector3d &reverse();

  vector3d operator-() const;

  vector3d &operator+=(const vector3d &from);
  vector3d &operator-=(const vector3d &from);

  vector3d &operator+=(double scalar);
  vector3d &operator-=(double scalar);
  vector3d &operator*=(double scalar);
  vector3d &operator/=(double scalar);

  double   length() const;
  vector3d cross(const vector3d &from) const;
};

vector3d operator*(double scalar, const vector3d &from);
vector3d operator*(const vector3d &lhs, double scalar);
vector3d operator/(const vector3d &lhs, double scalar);

vector3d operator+(const vector3d &lhs, const vector3d &rhs);
vector3d operator-(const vector3d &lhs, const vector3d &rhs);

//----------------------------------------------------------------------------
inline vector3d vector3d::cross(const vector3d &from) const
{
  return vector3d(y * from.z - z * from.y, z * from.x - x * from.z, x * from.y - y * from.x);
}
//----------------------------------------------------------------------------
inline vector3d &vector3d::operator+=(const vector3d &from)
{
  x += from.x;
  y += from.y;
  z += from.z;
  return *this;
}
//----------------------------------------------------------------------------
inline vector3d &vector3d::operator-=(const vector3d &from)
{
  x -= from.x;
  y -= from.y;
  z -= from.z;
  return *this;
}
//----------------------------------------------------------------------------
inline vector3d &vector3d::operator*=(double scalar)
{
  x *= scalar;
  y *= scalar;
  z *= scalar;
  return *this;
}

//----------------------------------------------------------------------------
inline vector3d &vector3d::operator+=(double scalar)
{
  x += scalar;
  y += scalar;
  z += scalar;
  return *this;
}

//----------------------------------------------------------------------------
inline vector3d &vector3d::operator-=(double scalar)
{
  x -= scalar;
  y -= scalar;
  z -= scalar;
  return *this;
}

#endif
