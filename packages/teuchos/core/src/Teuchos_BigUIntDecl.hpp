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

#ifndef TEUCHOS_BIG_UINT_DECL_HPP
#define TEUCHOS_BIG_UINT_DECL_HPP

#include <iosfwd>

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
