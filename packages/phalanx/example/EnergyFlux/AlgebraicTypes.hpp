// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_ENTITY_TEMPLATES
#define PHX_ENTITY_TEMPLATES

#include <iostream>
#include "Phalanx_ConfigDefs.hpp"
#include "Sacado.hpp"

// *********************************************************************
// Vector
// *********************************************************************

template <typename T> class MyVector {
  T val[3];
public:
  MyVector() { }
  MyVector(const T& value)
  {
    this->init(value);
  }
  MyVector(const T& v0, const T& v1, const T& v2)
  {
    val[0] = v0;
    val[1] = v1;
    val[2] = v2;
  }
  void init(const T& value) 
  {
    for (std::size_t i=0; i < 3; ++i)
      val[i] = value;
  }
  T& operator[](std::size_t i)
  {
    return val[i];
  }
  const T& operator[](std::size_t i) const
  {
    return val[i];
  }
  MyVector<T> operator+()
  {
    return *this;
  }
  MyVector<T> operator-()
  {
    return MyVector<T>(-val[0], -val[1], -val[2]);
  }
  MyVector<T>& operator+=(const MyVector<T>& a)
  {
    for (std::size_t i=0; i < 3; ++i)
      val[i] += a[i];
    return *this;
  }
  MyVector<T>& operator-=(const MyVector<T>& a)
  {
    for (std::size_t i=0; i < 3; ++i)
      val[i] -= a[i];
    return *this;
  }
  MyVector<T>& operator*=(const MyVector<T>& a)
  {
    for (std::size_t i=0; i < 3; ++i)
      val[i] *= a.val[i];
    return *this;
  }
  MyVector<T>& operator/=(const MyVector<T>& a)
  {
    for (std::size_t i=0; i < 3; ++i)
      val[i] /= a.val[i];
    return *this;
  }
  void print(std::ostream& os) const
  {
    for (std::size_t i=0; i < 3; ++i)
      os << "MyVector[" << i << "] = " << val[i] << std::endl;
  }

};

// MyVector/Scalar operations
template <class T, class U>
inline MyVector< typename Sacado::Promote<T,U>::type>
operator*(const MyVector<T>& v, const U& sc)
{
  typedef typename Sacado::Promote<T,U>::type ValueT;
  return MyVector<ValueT>(v[0] * sc, v[1] * sc, v[2] * sc);
}

template <class T, class U>
inline MyVector< typename Sacado::Promote<T,U>::type>
operator*(const T& sc, const MyVector<U>& u)
{
  typedef typename Sacado::Promote<T,U>::type ValueT;
  return MyVector<ValueT>(u[0] * sc, u[1] * sc, u[2] * sc);
}

// MyVector/MyVector operations
template <class T, class U>
inline MyVector< typename Sacado::Promote<T,U>::type>
operator*(const MyVector<T>& v , const MyVector<U>& u)
{
  typedef typename Sacado::Promote<T,U>::type ValueT;
  return MyVector<ValueT>(v[0] * u[0], v[1] * u[1], v[2] * u[2]);
}

template <class T, class U>
inline MyVector< typename Sacado::Promote<T,U>::type>
operator/(const MyVector<T>& v , const MyVector<U>& u)
{
  typedef typename Sacado::Promote<T,U>::type ValueT;
  return MyVector<ValueT>(v[0] / u[0], v[1] / u[1], v[2] / u[2]);
}

template <class T, class U>
inline MyVector< typename Sacado::Promote<T,U>::type>
operator+(const MyVector<T>& v , const MyVector<U>& u)
{
  typedef typename Sacado::Promote<T,U>::type ValueT;
  return MyVector<ValueT>(v[0] + u[0], v[1] + u[1], v[2] + u[2]);
}

template <class T, class U>
inline MyVector< typename Sacado::Promote<T,U>::type>
operator-(const MyVector<T>& v , const MyVector<U>& u)
{
  typedef typename Sacado::Promote<T,U>::type ValueT;
  return MyVector<ValueT>(v[0] - u[0], v[1] - u[1], v[2] - u[2]);
}

// Printing
template <typename T>
std::ostream& operator<<(std::ostream& os, const MyVector<T>& v)
{
  v.print(os);
  return os;
}

// *********************************************************************
// MyTensor
// *********************************************************************
template <typename T> class MyTensor {
  T val[9];
public:
  void init(const T& value) 
  {
    for (std::size_t i=0; i < 9; ++i)
      val[i] = value;
  }
  T& operator[](std::size_t i)
  {
    return val[i];
  }
  MyTensor<T>& operator*(MyVector<T>& a)
  {
    for (std::size_t i=0; i < 9; ++i)
      val[i] = val[i] * a[i];
    return *this;
  }
  void print() 
  {
    for (std::size_t i=0; i < 9; ++i)
      std::cout << "MyTensor[" << i << "] = " << val[i] << std::endl;
  }
};

// *********************************************************************
// *********************************************************************
						 
#endif
