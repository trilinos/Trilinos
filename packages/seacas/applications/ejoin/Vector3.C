// Magic Software, Inc.
// http://www.magic-software.com
// Copyright(C) 2010 Sandia Corporation.
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Vector3.h"
#include <cmath>

//----------------------------------------------------------------------------
Vector3::Vector3 ()
  : x(0.0), y(0.0), z(0.0)
{
}

//----------------------------------------------------------------------------
Vector3::Vector3 (double fX, double fY, double fZ)
  : x(fX), y(fY), z(fZ)
{}

//----------------------------------------------------------------------------
Vector3::Vector3 (double Coordinate[3])
  : x(Coordinate[0]), y(Coordinate[1]), z(Coordinate[2])
{}

//----------------------------------------------------------------------------
Vector3::Vector3 (const Vector3& from)
  : x(from.x), y(from.y), z(from.z)
{}

void Vector3::set(double fX, double fY, double fZ)
{
    x = fX;
    y = fY;
    z = fZ;
}

void Vector3::set(double Coordinate[3])
{
    x = Coordinate[0];
    y = Coordinate[1];
    z = Coordinate[2];
}

Vector3& Vector3::operator= (const Vector3& from)
{
    x = from.x;
    y = from.y;
    z = from.z;
    return *this;
}

Vector3& Vector3::reverse()
{
  x = -x;
  y = -y;
  z = -z;
  return *this;
}


bool Vector3::operator== (const Vector3& from) const
{
    return ( x == from.x && y == from.y && z == from.z );
}

bool Vector3::operator!= (const Vector3& from) const
{
    return ( x != from.x || y != from.y || z != from.z );
}

const Vector3 operator+ (const Vector3& lhs, const Vector3& rhs)
{
  Vector3 Sum(lhs);
  return Sum += rhs;
}

const Vector3 operator- (const Vector3& lhs, const Vector3& rhs)
{
  Vector3 Diff(lhs);
  return Diff -= rhs;
}

const Vector3 operator* (const Vector3& lhs, double scalar)
{
  Vector3 Prod(lhs);
  return Prod *= scalar;
}

Vector3 Vector3::operator- () const
{
    Vector3 Neg;
    return Neg *= -1.0;
}

const Vector3 operator* (double scalar, const Vector3& from)
{
  Vector3 Prod(from);
  return Prod *= scalar;
}

double Vector3::squared_length () const
{
    return x*x + y*y + z*z;
}

double Vector3::dot (const Vector3& from) const
{
    return x*from.x + y*from.y + z*from.z;
}

const Vector3 operator/ (const Vector3& lhs, double scalar)
{
  Vector3 Quot(lhs);
  return Quot /= scalar;
}

Vector3& Vector3::operator/= (double scalar)
{
  if ( scalar != 0.0 ) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
  } else {
    x = HUGE_VAL;
    y = HUGE_VAL;
    z = HUGE_VAL;
  }

  return *this;
}

double Vector3::length () const
{
  return std::sqrt(x*x + y*y + z*z);
}

double Vector3::normalize (double tolerance)
{
  double mylength = length();

  if ( mylength > tolerance ) {
    x /= mylength;
    y /= mylength;
    z /= mylength;
  } else {
    mylength = 0.0;
  }

  return mylength;
}

Vector3 Vector3::plane_normal(const Vector3 &v1,
			      const Vector3 &v2,
			      const Vector3 &v3)
{
  Vector3 v32 = v3;
  v32 -= v2;
  Vector3 v12 = v1;
  v12 -= v2;
  return v32.cross(v12);
}
