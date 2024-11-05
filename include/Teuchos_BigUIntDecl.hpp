// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_BIG_UINT_DECL_HPP
#define TEUCHOS_BIG_UINT_DECL_HPP

#include <iosfwd>
#include <cstdint>

/*! \file Teuchos_BigUIntDecl.hpp
    \brief Arbitrary-precision unsigned integer declaration.
*/

namespace Teuchos {

/** \brief Arbitrary-precision unsigned integer class.

This class implements an unsigned integer type of arbitrary precision.
The precision is chosen at compile time via a template parameter.
The template parameter specifies how many 32-bit "digits" will compose
the full integer. Thus, the number has (32*n) bits of precision.

This class was primarily created to serve the Teuchos::print_double
function for printing floating-point values in a lossless and minimal way.
*/

template <int n>
class BigUInt {
 private:
  std::uint32_t x[n];
 public:
  BigUInt();
  BigUInt(std::uint64_t v);
  BigUInt(std::string const& s);
  explicit operator bool() const noexcept;
  std::uint32_t& operator[](int i);
  std::uint32_t const& operator[](int i) const;
  BigUInt& operator+=(std::uint32_t b);
  BigUInt& operator+=(BigUInt const& b);
  BigUInt& operator-=(std::uint32_t b);
  BigUInt& operator-=(BigUInt const& b);
  BigUInt& operator*=(std::uint32_t b);
  BigUInt& operator<<=(std::uint32_t b);
  BigUInt& operator>>=(std::uint32_t b);
};

template <int n>
std::ostream& operator<<(std::ostream& os, BigUInt<n> a);
template <int n>
BigUInt<n> operator+(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
BigUInt<n> operator-(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
BigUInt<n> operator*(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
BigUInt<n> operator/(BigUInt<n> const& a, std::uint32_t const& b);
template <int n>
BigUInt<n> operator/(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
int comp(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
bool operator>=(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
bool operator<=(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
bool operator<(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
bool operator>(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
bool operator==(BigUInt<n> const& a, BigUInt<n> const& b);
template <int n>
void divmod(BigUInt<n>& quotient, BigUInt<n>& x, std::uint32_t const& b);

}

#endif
