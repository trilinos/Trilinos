#ifndef MGCVECTOR3_H
#define MGCVECTOR3_H

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

  class Vector3
  {
  public:
    // construction
    Vector3 ();
    Vector3 (double fX, double fY, double fZ);
    explicit Vector3 (double Coordinate[3]);
    Vector3 (const Vector3& from);

    double x, y, z;

    // assignment and comparison
    Vector3& operator= (const Vector3& from);
    bool operator== (const Vector3& from) const;
    bool operator!= (const Vector3& from) const;
    void set(double fX, double fY, double fZ);
    void set(double Coordinate[3]);
    Vector3& reverse();

    // arithmetic operations
    Vector3 operator- () const;

    // arithmetic updates
    Vector3& operator+= (const Vector3& from);
    Vector3& operator+= (double scalar);
    Vector3& operator-= (const Vector3& from);
    Vector3& operator-= (double scalar);
    Vector3& operator*= (double scalar);
    Vector3& operator/= (double scalar);

    // vector operations
    double length () const;
    double squared_length () const;
    double dot (const Vector3& from) const;
    double normalize (double tolerance = 1e-06);
    Vector3 cross (const Vector3& from) const;
    static Vector3 plane_normal(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);
  };

  const Vector3 operator* (double scalar, const Vector3& vec);
  const Vector3 operator+ (const Vector3& vec, const Vector3& vec2);
  const Vector3 operator- (const Vector3& vec, const Vector3& vec2);
  const Vector3 operator* (double scalar, const Vector3& vec);
  const Vector3 operator/ (const Vector3& vec, double scalar);

  //----------------------------------------------------------------------------
  inline Vector3 Vector3::cross (const Vector3& from) const
  {
    return Vector3(y*from.z - z*from.y,
		   z*from.x - x*from.z,
		   x*from.y - y*from.x);
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator+= (const Vector3& from)
  {
    x += from.x;
    y += from.y;
    z += from.z;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator+= (double scalar)
  {
    x += scalar;
    y += scalar;
    z += scalar;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator-= (const Vector3& from)
  {
    x -= from.x;
    y -= from.y;
    z -= from.z;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator-= (double scalar)
  {
    x -= scalar;
    y -= scalar;
    z -= scalar;
    return *this;
  }
  //----------------------------------------------------------------------------
  inline Vector3& Vector3::operator*= (double scalar)
  {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
  }
#endif

