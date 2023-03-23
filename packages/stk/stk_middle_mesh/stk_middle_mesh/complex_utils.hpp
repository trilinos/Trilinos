#ifndef COMPLEX_UTILS_H
#define COMPLEX_UTILS_H

#include <complex>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// the standard library does not support mixed integer-std::complex
// arithmetic, add that support here

template <typename T>
std::complex<T> operator+(const std::complex<T>& v1, const int& v2)
{
  return std::complex<T>(v1.real() + v2, v1.imag());
}

template <typename T>
std::complex<T> operator+(const int& v1, const std::complex<T>& v2)
{
  return v2 + v1;
}

template <typename T>
std::complex<T> operator-(const std::complex<T>& v1, const int& v2)
{
  return std::complex<T>(v1.real() - v2, v1.imag());
}

template <typename T>
std::complex<T> operator-(const int& v1, const std::complex<T>& v2)
{
  return std::complex<T>(v1 - v2.real(), -v2.imag());
}

template <typename T>
std::complex<T> operator*(const std::complex<T>& v1, const int& v2)
{
  return std::complex<T>(v1.real() * v2, v1.imag() * v2);
}

template <typename T>
std::complex<T> operator*(const int& v1, const std::complex<T>& v2)
{
  return v2 * v1;
}

template <typename T>
std::complex<T> operator/(const std::complex<T>& v1, const int& v2)
{
  return std::complex<T>(v1.real() / v2, v1.imag() / v2);
}

template <typename T>
std::complex<T> operator/(const int& v1, const std::complex<T>& v2)
{
  return std::complex<T>(v1, 0) / v2;
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
