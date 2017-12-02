/*
 * Copyright(C) 1999-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef VECTOR3D
#define VECTOR3D

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
  bool operator==(const vector3d &from) const;
  bool operator!=(const vector3d &from) const;
  void set(double X, double Y, double Z);
  void set(double location[3]);
  vector3d &reverse();

  vector3d operator-() const;

  vector3d &operator+=(const vector3d &from);
  vector3d &operator-=(const vector3d &from);
  vector3d &operator*=(double scalar);
  vector3d &operator/=(double scalar);

  double length() const;
  double normalize(double tolerance = 1e-06);
  vector3d cross(const vector3d &from) const;
  static vector3d plane_normal(const vector3d &v1, const vector3d &v2, const vector3d &v3);
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

#endif
