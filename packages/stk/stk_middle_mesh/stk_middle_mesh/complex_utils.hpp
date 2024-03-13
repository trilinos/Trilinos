#ifndef STK_MIDDLE_MESH_COMPLEX_UTILS_H
#define STK_MIDDLE_MESH_COMPLEX_UTILS_H

#include <complex>

namespace std
{
//-----------------------------------------------------------------------------
// The standard library doesn't have overloads for std::complex and integer operations

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator+(const std::complex<T>& a, T2 b)
{
  return {a.real() + b, a.imag()};
}

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator+(T2 b, const std::complex<T>& a)
{
  return {a.real() + b, a.imag()};
}

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator-(const std::complex<T>& a, T2 b)
{
  return {a.real() - b, a.imag()};
}

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator-(T2 b, const std::complex<T>& a)
{
  return {b - a.real(), -a.imag()};
}

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator*(const std::complex<T>& a, T2 b)
{
  return {a.real() * b, a.imag() * b};
}

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator*(T2 b, const std::complex<T>& a)
{
  return {a.real() * b, a.imag() * b};
}

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator/(const std::complex<T>& a, T2 b)
{
  return {a.real() / b, a.imag() / b};
}

template <typename T, typename T2, std::enable_if_t<std::is_integral_v<T2>, bool> = true>
std::complex<T> operator/(T2 b, const std::complex<T>& a)
{
  return std::complex<T>(b, 0) / a;
}


}
#endif
