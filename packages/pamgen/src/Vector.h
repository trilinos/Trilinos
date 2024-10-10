// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SNL_VLIB
#define SNL_VLIB

#include <math.h>
#include <iostream>

#include "pamgen_code_types.h"
namespace PAMGEN_NEVADA {

/*****************************************************************************/
struct Vector 
/*****************************************************************************/
// Represents a Cartesian vector of the appropriate dimensionality.
{
  Real x;
  Real y;
  Real z;
                /**** Constructors, destructors and assigns ****/
    Vector();
    Vector(Real xx, Real yy, Real zz = 0.0);

    Vector(const Vector &src);

    Vector& operator=(const Vector &src);

                /**** Access functions ****/
    const Real * operator() (void) const { return &x; }
          Real * operator() (void)       { return &x; }
    Real * operator() (long long i) { return (   (Real *) & x)+i;}
    Real& X() { return x;} 
    Real& Y() { return y;}
    const Real& X() const { return x;} 
    const Real& Y() const { return y;}


  Real & Z() { return z;}
  const Real & Z() const { return z;}


    void XYZ(Real& xx, Real& yy, Real& zz) const;

                /**** Definition functions ****/

    void X(Real x);
    void Y(Real y);
    void Z(Real z);

    void XYZ(Real xx, Real yy, Real zz = 0.0);
 
                /**** Operator overloads ****/

// Unary operations

    Vector  operator-() const;
    Vector  operator+() const;

// Operations with scalars

    Vector& operator*=(Real);
    Vector& operator/=(Real);

// Operations with other vectors

    Vector& operator+=(const Vector&);
    Vector& operator-=(const Vector&);
};

/*****************************************************************************/
/**************************** Vector *****************************************/

                /**** Constructors, destructor, assigns ****/

/*****************************************************************************/
inline Vector::Vector()
/*****************************************************************************/
  : x(0.0)
  , y(0.0)
  , z(0.0)
{}

/*****************************************************************************/
inline Vector::Vector( Real xx, Real yy, Real zz)
/*****************************************************************************/
  : x(xx)
  , y(yy)
  , z(zz)
{ }

/*****************************************************************************/
inline Vector::Vector(const Vector& src) 
/*****************************************************************************/
:
  x(src.x)
, y(src.y)
, z(src.z)
{}

/*****************************************************************************/
inline Vector& Vector::operator=(const Vector& src)
/*****************************************************************************/
{
  x = src.x;
  y = src.y;
  z = src.z;
  return *this;
}

/*****************************************************************************/
inline void Vector::XYZ(Real& xx, Real& yy, Real& zz) const 
/*****************************************************************************/
{
  xx = x;
  yy = y;
  zz = z;
}
                /**** Definition functions ****/

/*****************************************************************************/
inline void Vector::X(Real sc)
/*****************************************************************************/
{ 
  x = sc; 
}

/*****************************************************************************/
inline void Vector::Y(Real sc)
/*****************************************************************************/
{ 
  y = sc;
}

/*****************************************************************************/
inline void Vector::Z(Real sc)
/*****************************************************************************/
{ 
  z = sc;
}

/*****************************************************************************/
inline void Vector::XYZ(Real xx, Real yy, Real zz)
/*****************************************************************************/
{
  x = xx;
  y = yy;
  z = zz;
}

                /**** Operator overloads ****/

// Unary operators

/*****************************************************************************/
inline Vector Vector::operator-() const
/*****************************************************************************/
{
  return Vector(-x
              , -y
              , -z
                );
}

/*****************************************************************************/
inline Vector Vector::operator+() const 
/*****************************************************************************/
{
  return *this;
}

// operators with scalars

/*****************************************************************************/
inline Vector operator*(const Vector& v, Real sc)
/*****************************************************************************/
{
  return Vector(v.x*sc
              , v.y*sc
              , v.z*sc
                );
}

/*****************************************************************************/
inline Vector  operator*(Real sc, const Vector& v)
/*****************************************************************************/
{
  return v*sc;
}

inline Vector& Vector::operator*=(Real sc){
  x *= sc;
  y *= sc;
  z *= sc;
  return *this;
}

/*****************************************************************************/
inline Vector operator/(const Vector& v, Real sc)
/*****************************************************************************/
{
  return Vector(v.x/sc
              , v.y/sc
              , v.z/sc
                );
}

/*****************************************************************************/
inline Vector& Vector::operator/=(Real sc)
/*****************************************************************************/
{
  x /= sc;
  y /= sc;
  z /= sc;
  return *this;
}

// Operators with vectors

/*****************************************************************************/
inline Real operator*(const Vector& va, const Vector& vb)
/*****************************************************************************/
{
  return va.x*vb.x
       + va.y*vb.y
       + va.z*vb.z
        ;
}

/*****************************************************************************/
inline Vector operator+(const Vector& va, const Vector& vb)
/*****************************************************************************/
{
  return Vector(va.x + vb.x
              , va.y + vb.y
              , va.z + vb.z
                );
}

/*****************************************************************************/
inline Vector& Vector::operator+=(const Vector& vb)
/*****************************************************************************/
{
  x += vb.x;
  y += vb.y;
  z += vb.z;
  return *this;
}

/*****************************************************************************/
inline Vector operator-(const Vector& va, const Vector& vb)
/*****************************************************************************/
{
  return Vector(va.x - vb.x
		, va.y - vb.y
		, va.z - vb.z
                );
}

/*****************************************************************************/
inline Vector& Vector::operator-=(const Vector& vb)
/*****************************************************************************/
{
  x -= vb.x;
  y -= vb.y;
  z -= vb.z;
  return *this;
}

/*****************************************************************************/
inline bool operator==(const Vector& va, const Vector &vb)
/*****************************************************************************/
{
  return va.x == vb.x
      && va.y == vb.y
      && va.z == vb.z
  ;
}

/*****************************************************************************/
inline bool operator!=(const Vector& va, const Vector &vb)
/*****************************************************************************/
{
  return va.x != vb.x
      || va.y != vb.y
      || va.z != vb.z
  ;
}
                /**** Methods ****/
// Stream i/o

/*****************************************************************************/
inline std::ostream& operator<<(std::ostream &str, const Vector &v)
/*****************************************************************************/
{
  str << '(' << v.x
      << ',' << v.y
      << ',' << v.z
      << ')';
  return str;
}
} // end namespace PAMGEN_NEVADA {
#endif
