// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_BIG_UINT_HPP
#define TEUCHOS_BIG_UINT_HPP

#include <iostream>
#include <Teuchos_BigUIntDecl.hpp>

/*! \file Teuchos_BigUInt.hpp
    \brief Arbitrary-precision unsigned integer definition.
*/

namespace Teuchos {

template <int n>
BigUInt<n>::BigUInt() {}

template <int n>
BigUInt<n>::BigUInt(std::uint64_t v) {
  for (int i = 2; i < n; ++i) x[i] = 0;
  if (n > 1) x[1] = std::uint32_t(v >> 32);
  x[0] = std::uint32_t(v);
}

template <int n>
BigUInt<n>::BigUInt(std::string const& s) : BigUInt(std::uint32_t(0)) {
  for (auto c : s) {
    operator*=(10);
    operator+=(c - '0');
  }
}

template <int n>
BigUInt<n>::operator bool() const noexcept {
  for (int i = 0; i < n; ++i) if (x[i]) return true;
  return false;
}

template <int n>
std::uint32_t& BigUInt<n>::operator[](int i) { return x[i]; }

template <int n>
std::uint32_t const& BigUInt<n>::operator[](int i) const { return x[i]; }

template <int n>
BigUInt<n>& BigUInt<n>::operator+=(std::uint32_t b) {
  std::uint32_t carry = b;
  for (int i = 0; i < n; ++i) {
    std::uint64_t ax = x[i];
    auto cx = ax + std::uint64_t(carry);
    carry = std::uint32_t(cx >> 32);
    x[i] = std::uint32_t(cx);
  }
  return *this;
}

template <int n>
BigUInt<n>& BigUInt<n>::operator+=(BigUInt<n> const& b) {
  std::uint32_t carry = 0;
  for (int i = 0; i < n; ++i) {
    std::uint64_t ax = x[i];
    std::uint64_t bx = b[i];
    auto cx = ax + bx + std::uint64_t(carry);
    carry = std::uint32_t(cx >> 32);
    x[i] = std::uint32_t(cx);
  }
  return *this;
}

template <int n>
BigUInt<n>& BigUInt<n>::operator-=(std::uint32_t b) {
  std::int64_t carry = b;
  for (int i = 0; i < n; ++i) {
    std::int64_t ax = x[i];
    auto cx = ax - carry;
    if (cx < 0) {
      carry = 1;
      cx += std::int64_t(1) << 32;
    } else {
      carry = 0;
    }
    x[i] = std::uint32_t(cx);
  }
  return *this;
}

template <int n>
BigUInt<n>& BigUInt<n>::operator-=(BigUInt<n> const& b) {
  std::int64_t carry = 0;
  for (int i = 0; i < n; ++i) {
    std::int64_t ax = x[i];
    std::int64_t bx = b[i];
    auto cx = ax - bx - carry;
    if (cx < 0) {
      carry = 1;
      cx += std::int64_t(1) << 32;
    } else {
      carry = 0;
    }
    x[i] = std::uint32_t(cx);
  }
  return *this;
}

template <int n>
BigUInt<n>& BigUInt<n>::operator*=(std::uint32_t b) {
  std::uint32_t carry = 0;
  for (int i = 0; i < n; ++i) {
    std::uint64_t ax = x[i];
    auto cx = (ax * std::uint64_t(b)) + std::uint64_t(carry);
    carry = std::uint32_t(cx >> 32);
    x[i] = std::uint32_t(cx);
  }
  return *this;
}

template <int n>
BigUInt<n>& BigUInt<n>::operator<<=(std::uint32_t b) {
  auto ndigits = b / 32;
  auto nbits = b - (ndigits * 32);
  for (int i = n - 1; i >= 0; --i) {
    std::uint32_t xi = 0;
    if (i >= int(ndigits)) {
      xi = x[i - ndigits] << nbits;
    }
    // nbits &&, because apparently shifting a 32-bit value by 32 is not allowed
    if (nbits && (i > int(ndigits))) {
      xi |= x[i - ndigits - 1] >> (32 - nbits);
    }
    x[i] = xi;
  }
  return *this;
}

template <int n>
BigUInt<n>& BigUInt<n>::operator>>=(std::uint32_t b) {
  auto ndigits = b / 32;
  auto nbits = b - (ndigits * 32);
  for (int i = 0; i < n; ++i) {
    std::uint32_t xi = 0;
    if (i + ndigits < n) xi = x[i + ndigits] >> nbits;
    if (nbits && i + ndigits + 1 < n) xi |= x[i + ndigits + 1] << (32 - nbits);
    x[i] = xi;
  }
  return *this;
}

template <int n>
std::ostream& operator<<(std::ostream& os, BigUInt<n> a) {
  char buf[n * 20];
  int i = 0;
  while (a) {
    BigUInt<n> quotient;
    divmod(quotient, a, 10);
    auto remainder = a[0];
    a = quotient;
    buf[i++] = char(remainder) + '0';
  }
  for (int j = 0; j < i / 2; ++j) {
    auto tmp = buf[i - j - 1];
    buf[i - j - 1] = buf[j];
    buf[j] = tmp;
  }
  if (i == 0) buf[i++] = '0';
  buf[i] = '\0';
  os << buf;
  return os;
}

template <int n>
BigUInt<n> operator+(BigUInt<n> const& a, BigUInt<n> const& b) {
  auto c = a;
  c += b;
  return c;
}

template <int n>
BigUInt<n> operator-(BigUInt<n> const& a, BigUInt<n> const& b) {
  auto c = a;
  c -= b;
  return c;
}

template <int n>
BigUInt<n> operator*(BigUInt<n> const& a, BigUInt<n> const& b) {
  BigUInt<n> a_times_b_i;
  BigUInt<n> c{0};
  for (int i = n - 1; i >= 0; --i) {
    a_times_b_i = a;
    a_times_b_i *= b[i];
    c <<= 32;
    c += a_times_b_i;
  }
  return c;
}

template <int n>
BigUInt<n> operator/(BigUInt<n> const& a, std::uint32_t const& b) {
  BigUInt<n> quotient;
  auto x = a;
  divmod(quotient, x, b);
  return quotient;
}

template <int n>
BigUInt<n> operator/(BigUInt<n> const& a, BigUInt<n> const& b) {
  if (b > a) return BigUInt<n>(0);
  BigUInt<n> quotient(1);
  auto c = b;
  while (c < a) {
    c <<= 1;
    quotient <<= 1;
  }
  auto factor = quotient;
  factor >>= 1;
  while (factor) {
    int d = comp(a, c);
    if (d == 0) break;
    if (d == -1) {
      c -= b * factor;
      quotient -= factor;
    } else {
      c += b * factor;
      quotient += factor;
    }
    factor >>= 1;
  }
  if (c > a) {
    c -= b;
    quotient -= 1;
  }
  return quotient;
}

template <int n>
int comp(BigUInt<n> const& a, BigUInt<n> const& b) {
  for (int i = n - 1; i >= 0; --i) {
    if (a[i] != b[i]) {
      if (a[i] > b[i]) return 1;
      else return -1;
    }
  }
  return 0;
}

template <int n>
bool operator>=(BigUInt<n> const& a, BigUInt<n> const& b) {
  return comp(a, b) > -1;
}

template <int n>
bool operator<=(BigUInt<n> const& a, BigUInt<n> const& b) {
  return comp(a, b) < 1;
}

template <int n>
bool operator<(BigUInt<n> const& a, BigUInt<n> const& b) {
  return comp(a, b) == -1;
}

template <int n>
bool operator>(BigUInt<n> const& a, BigUInt<n> const& b) {
  return comp(a, b) == 1;
}

template <int n>
bool operator==(BigUInt<n> const& a, BigUInt<n> const& b) {
  for (int i = 0; i < n; ++i) if (a[i] != b[i]) return false;
  return true;
}

template <int n>
void divmod(BigUInt<n>& quotient, BigUInt<n>& x, std::uint32_t const& b) {
  quotient = BigUInt<n>(std::uint32_t(0));
  for (int i = n - 1; i >= 0;) {
    if (x[i]) {
      if (x[i] >= b) {
        auto dividend = x[i];
        auto quotient2 = dividend / b;
        auto remainder = dividend - (quotient2 * b);
        quotient[i] = quotient2;
        x[i] = remainder;
      } else if (i > 0) {
        auto dividend = std::uint64_t(x[i]);
        dividend <<= 32;
        dividend |= x[i - 1];
        auto quotient2 = dividend / std::uint64_t(b);
        auto remainder = dividend - (quotient2 * std::uint64_t(b));
        quotient[i - 1] = std::uint32_t(quotient2);
        x[i] = 0;
        x[i - 1] = std::uint32_t(remainder);
        i = i - 1;
      } else {
        break;
      }
    } else {
      i = i - 1;
    }
  }
}

}

#endif

