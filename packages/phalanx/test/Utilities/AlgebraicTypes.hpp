#ifndef PHX_ENTITY_TEMPLATES
#define PHX_ENTITY_TEMPLATES

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
  MyVector<T>& operator*(const T& a)
  {
    for (std::size_t i=0; i < 3; ++i)
      val[i] = val[i] * a;
    return *this;
  }
  MyVector<T>& operator/(const T& a)
  {
    for (std::size_t i=0; i < 3; ++i)
      val[i] = val[i] / a;
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
  v.template print(os);
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
