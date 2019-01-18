/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef Stk_Math_Vector_hpp
#define Stk_Math_Vector_hpp

#include <array>
#include <cmath>
#include <type_traits>
#include <ostream>
#include <iterator>
#include <cassert>
#include <string>
#include <limits>

namespace stk {
namespace math {

enum class MemberInit { NONE }; // Only entry should be MemberInit::NONE

template<class REAL, unsigned DIM>
class Vec {

typedef typename std::array<REAL,DIM>::const_iterator const_iterator;
typedef typename std::array<REAL,DIM>::iterator iterator;

 public:
  typedef REAL Real;
  static const Vec<REAL,DIM> ZERO;


  Vec(const Vec<REAL,DIM>& rhs) { vec = rhs.vec; }
  Vec(Vec<REAL,DIM>&&) = default;

  explicit Vec(const double * rhs) { assert(rhs); for (unsigned i=0; i<DIM; ++i) vec[i] = rhs[i]; }
  Vec(const double * rhs, const unsigned len) { assert(len == 0 || rhs); assert(DIM >= len); for (unsigned i=0; i<len; ++i) vec[i] = rhs[i]; for (unsigned i=len; i<DIM; ++i) vec[i] = 0; }
  Vec(const double x, const double y, const double z) { static_assert(DIM==3, "Invalid dimension"); vec[0] = x; vec[1] = y; vec[2] = z; }
  Vec(const double x, const double y) { static_assert(DIM==2, "Invalid dimension"); vec[0] = x; vec[1] = y; }
  explicit Vec(const double x) { static_assert(DIM==1, "Invalid dimension"); vec[0] = x; }
  constexpr Vec() : vec{} {}
  Vec(MemberInit) {/* MemberInit::NONE */}

  Vec<REAL,DIM>& operator = (const Vec<REAL,DIM>& rhs) { vec = rhs.vec; return *this; }
  Vec<REAL,DIM>& operator = (Vec<REAL,DIM>&&) = default;
  Vec<REAL,DIM>& operator -= (const Vec<REAL,DIM>& rhs) { for (unsigned i=0; i<DIM; ++i) vec[i] -= rhs.vec[i]; return *this; }
  Vec<REAL,DIM>& operator += (const Vec<REAL,DIM>& rhs) { for (unsigned i=0; i<DIM; ++i) vec[i] += rhs.vec[i]; return *this; }
  Vec<REAL,DIM>& operator *= (const REAL rhs) { for (unsigned i=0; i<DIM; ++i) vec[i] *= rhs; return *this; }
  Vec<REAL,DIM>& operator /= (const REAL rhs) { for (unsigned i=0; i<DIM; ++i) vec[i] /= rhs; return *this; }

  Vec<REAL,DIM> operator - (void) const { Vec<REAL,DIM> result(MemberInit::NONE); for (unsigned i=0; i<DIM; ++i) result.vec[i] = -vec[i]; return result; }
  Vec<REAL,DIM> operator + (const Vec<REAL,DIM>& rhs) const { Vec<REAL,DIM> result(MemberInit::NONE); for (unsigned i=0; i<DIM; ++i) result.vec[i] = vec[i]+rhs.vec[i]; return result; }
  Vec<REAL,DIM> operator - (const Vec<REAL,DIM>& rhs) const { Vec<REAL,DIM> result(MemberInit::NONE); for (unsigned i=0; i<DIM; ++i) result.vec[i] = vec[i]-rhs.vec[i]; return result; }

  bool zero_length() const; ///< Return true if vector is zero length

  bool operator == (const Vec<REAL,DIM>& rhs) const;
  bool operator != (const Vec<REAL,DIM>& rhs) const;
  bool operator < (const Vec<REAL,DIM>& rhs) const;

  REAL unitize(); ///< Turn vector into a unit length vector
  REAL length() const; ///< Vector magnitude
  REAL length_squared() const; ///< Magnitude squared (faster, avoids square root)
  Vec<REAL,DIM> unit_vector() const; ///< Return unit length vector leaving original unchanged

  const REAL & operator[](const unsigned i) const {return vec[i];}
  REAL & operator[](const unsigned i) {return vec[i];}
  const REAL * data() const {return vec.data(); }
  REAL * data() {return vec.data(); }

  friend std::ostream& operator<<( std::ostream& out, const Vec<REAL,DIM>& rhs )
  {
    out << "Vec" << DIM << "d: ";
    for (unsigned i=0; i<DIM; ++i) out << rhs.vec[i] << " ";
    return out;
  }

  std::string to_string() const
  {
      std::string output;
      for ( size_t i=0; i<DIM; i++)
      {
          output += std::to_string(vec[i]);
          if ( i != DIM-1 )
              output += " ";
      }
      return output;
  }

  unsigned dimension() const { return DIM; }
  unsigned size() const { return DIM; }

  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;

  void set_invalid() { for (unsigned i=0; i<DIM; ++i) vec[i] = std::nan(""); }
  bool is_valid() const { for ( REAL item : vec ) { if ( std::isnan(item) ) return false; } return true; }

 private:
  std::array<REAL, DIM> vec;
};

template<class REAL, unsigned DIM>
const Vec<REAL,DIM> Vec<REAL,DIM>::ZERO = Vec<REAL,DIM>();

template<class REAL, unsigned DIM>
inline bool Vec<REAL,DIM>::zero_length() const
{
  // This is a strict equality!
  for (unsigned i=0; i<DIM; ++i)
  {
    if (vec[i] != 0.0) return false;
  }
  return true;
}

template<class REAL, unsigned DIM>
inline bool Vec<REAL,DIM>::operator == (const Vec<REAL,DIM>& rhs) const
{
  // This is a strict equality!
  for (unsigned i=0; i<DIM; ++i)
  {
    if (vec[i] != rhs.vec[i]) return false;
  }
  return true;
}

template<class REAL, unsigned DIM>
inline bool Vec<REAL,DIM>::operator != (const Vec<REAL,DIM>& rhs) const
{
  for (unsigned i=0; i<DIM; ++i)
  {
    if (vec[i] != rhs.vec[i]) return true;
  }
  return false;
}

template<class REAL, unsigned DIM>
inline bool Vec<REAL,DIM>::operator < (const Vec<REAL,DIM>& rhs) const
{
  for (unsigned i=0; i<DIM; ++i)
  {
    if (vec[i] < rhs.vec[i]) return true;
    if (vec[i] > rhs.vec[i]) return false;
  }
  return false;
}

template<class REAL, unsigned DIM>
inline REAL Vec<REAL,DIM>::unitize()
{
  const REAL len = length();
  if (len > 0.0)
  {
    const REAL inv_length = 1.0/len;
    for (unsigned i=0; i<DIM; ++i) vec[i] *= inv_length;
  }
  else
  {
    set_invalid();
  }
  return len;
}

template<class REAL, unsigned DIM>
inline REAL Vec<REAL,DIM>::length_squared() const
{
  REAL lengthSquared = 0.0;
  for (unsigned i=0; i<DIM; ++i) lengthSquared += vec[i]*vec[i];
  return lengthSquared;
}

template<class REAL, unsigned DIM>
inline REAL Vec<REAL,DIM>::length() const
{
  REAL lengthSquared = 0.0;
  for (unsigned i=0; i<DIM; ++i) lengthSquared += vec[i]*vec[i];
  return std::sqrt(lengthSquared);
}

template<class REAL, unsigned DIM>
inline Vec<REAL,DIM> Vec<REAL,DIM>::unit_vector() const
{
  Vec<REAL,DIM> result(MemberInit::NONE);
  for (unsigned i=0; i<DIM; ++i) result.vec[i] = vec[i];
  result.unitize();
  return result;
}

template<class REAL, unsigned DIM, class T>
inline Vec<REAL,DIM> operator * (const T scalar, const Vec<REAL,DIM>& rhs) {
  Vec<REAL,DIM> result(MemberInit::NONE);
  for (unsigned i=0; i<DIM; ++i) result[i] = scalar*rhs[i];
  return result;
}

template<class REAL, unsigned DIM, class T>
inline Vec<REAL,DIM> operator * (const Vec<REAL,DIM>& rhs, const T scalar) {
  Vec<REAL,DIM> result(MemberInit::NONE);
  for (unsigned i=0; i<DIM; ++i) result[i] = scalar*rhs[i];
  return result;
}

template<class REAL, unsigned DIM, class T>
inline Vec<REAL,DIM> operator / (const Vec<REAL,DIM>& rhs, const T scalar) {
  Vec<REAL,DIM> result(MemberInit::NONE);
  for (unsigned i=0; i<DIM; ++i) result[i] = rhs[i]/scalar;
  return result;
}

template<class REAL, unsigned DIM>
inline REAL Dot(const Vec<REAL,DIM>& a, const Vec<REAL,DIM>& b) {
  REAL dot = 0.0;
  for (unsigned i=0; i<DIM; ++i) dot += a[i]*b[i];
  return dot;
}

template<class REAL, unsigned DIM>
inline Vec<REAL,DIM> Cross(const Vec<REAL,DIM>& a, const Vec<REAL,DIM>& b) {
  static_assert(DIM==3, "Invalid dimension");
  return Vec<REAL,DIM>(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

/// Optimized Cross Products to Unit Vectors
template<class REAL, unsigned DIM>
inline Vec<REAL,DIM> crossX(const Vec<REAL,DIM> &v) {
  static_assert(DIM==3, "Invalid dimension");
  return Vec<REAL,DIM>(0.0, v[2], -v[1]);
}

template<class REAL, unsigned DIM>
inline Vec<REAL,DIM> crossY(const Vec<REAL,DIM> &v) {
  static_assert(DIM==3, "Invalid dimension");
  return Vec<REAL,DIM>(-v[2], 0.0, v[0]);
}

template<class REAL, unsigned DIM>
inline Vec<REAL,DIM> crossZ(const Vec<REAL,DIM> &v) {
  static_assert(DIM==3, "Invalid dimension");
  return Vec<REAL,DIM>(v[1], -v[0], 0.0);
}

template<class REAL, unsigned DIM>
inline typename Vec<REAL,DIM>::iterator Vec<REAL,DIM>::begin()
{
  return this->vec.begin();
}

template<class REAL, unsigned DIM>
inline typename Vec<REAL,DIM>::const_iterator Vec<REAL,DIM>::begin() const
{
  return this->vec.begin();
}

template<class REAL, unsigned DIM>
inline typename Vec<REAL,DIM>::iterator Vec<REAL,DIM>::end()
{
  return this->vec.end();
}

template<class REAL, unsigned DIM>
inline typename Vec<REAL,DIM>::const_iterator Vec<REAL,DIM>::end() const
{
  return this->vec.end();
}



typedef Vec<double,3> Vector3d;
typedef Vec<double,2> Vector2d;
typedef Vec<double,1> Vector1d;
typedef Vec<float,3> Float3d;

}} 

#endif 
