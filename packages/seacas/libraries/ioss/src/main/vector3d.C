// Copyright (c) 2014, Sandia Corporation.
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
// 
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
// 

#include <math.h>

#include "vector3d.h"

//----------------------------------------------------------------------------
vector3d::vector3d ()
  : x(0.0), y(0.0), z(0.0)
{
}

//----------------------------------------------------------------------------
vector3d::vector3d (double X, double Y, double Z)
  : x(X), y(Y), z(Z)
{}

//----------------------------------------------------------------------------
vector3d::vector3d (double location[3])
  : x(location[0]), y(location[1]), z(location[2])
{}

//----------------------------------------------------------------------------
vector3d::vector3d (const vector3d& from)
  : x(from.x), y(from.y), z(from.z)
{}

void vector3d::set(double X, double Y, double Z)
{
    x = X;
    y = Y;
    z = Z;
}

void vector3d::set(double location[3])
{
    x = location[0];
    y = location[1];
    z = location[2];
}

vector3d& vector3d::operator= (const vector3d& from)
{
    x = from.x;
    y = from.y;
    z = from.z;
    return *this;
}

vector3d& vector3d::reverse()
{
  x = -x;
  y = -y;
  z = -z;
  return *this;
}


bool vector3d::operator== (const vector3d& from) const
{
    return ( x == from.x && y == from.y && z == from.z );
}

bool vector3d::operator!= (const vector3d& from) const
{
    return ( x != from.x || y != from.y || z != from.z );
}

vector3d operator+ (const vector3d& lhs, const vector3d& rhs)
{
  vector3d tmp(lhs);
  return tmp += rhs;
}

vector3d operator- (const vector3d& lhs, const vector3d& rhs)
{
  vector3d tmp(lhs);
  return tmp -= rhs;
}

vector3d operator* (const vector3d& lhs, double scalar)
{
  vector3d tmp(lhs);
  return tmp *= scalar;
}

vector3d vector3d::operator- () const
{
  vector3d tmp(x, y, z);
  return tmp *= -1.0;
}

vector3d operator* (double scalar, const vector3d& from)
{
  vector3d tmp(from);
  return tmp *= scalar;
}

double vector3d::squared_length () const
{
    return x*x + y*y + z*z;
}

double vector3d::dot (const vector3d& from) const
{
    return x*from.x + y*from.y + z*from.z;
}

vector3d operator/ (const vector3d& lhs, double scalar)
{
  if ( scalar != 0.0 ) {
    vector3d tmp(lhs);
    return tmp /= scalar;
  } else {
    return vector3d(HUGE_VAL, HUGE_VAL, HUGE_VAL);
  }
}

vector3d& vector3d::operator/= (double scalar)
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

double vector3d::length () const
{
  return sqrt(x*x + y*y + z*z);
}

double vector3d::normalize (double tolerance)
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

vector3d vector3d::plane_normal(const vector3d &v1,
				const vector3d &v2,
				const vector3d &v3)
{
  vector3d v32 = v3;   v32 -= v2;
  vector3d v12 = v1;   v12 -= v2;
  return v32.cross(v12);
}
