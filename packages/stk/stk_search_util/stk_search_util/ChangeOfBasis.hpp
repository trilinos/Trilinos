// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
//

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CHANGE_OF_BASIS_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_CHANGE_OF_BASIS_HPP_

#include <array>

namespace stk::search::impl {

inline std::array<double, 9> matmat(const std::array<double, 9>& a, const std::array<double, 9>& b)
{
  std::array<double, 9> c = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (unsigned i=0; i < 3; ++i)
  {
    for (unsigned j=0; j < 3; ++j)
    {
      for (unsigned k=0; k < 3; ++k)
      {
        c[3*i + j] += a[3*i + k] * b[3*k + j];
      }
    }
  }

  return c;
}

inline std::array<double, 3> matvec(const std::array<double, 9>& a, const std::array<double, 3>& x)
{
  std::array<double, 3> b = {0, 0, 0};
  for (unsigned i=0; i < 3; ++i)
  {
    for (unsigned j=0; j < 3; ++j)
    {
      b[i] += a[3*i + j] * x[j];
    }
  }

  return b;
}

inline std::array<double, 3> matvec_transposed(const std::array<double, 9>& a, const std::array<double, 3>& x)
{
  std::array<double, 3> b = {0, 0, 0};
  for (unsigned i=0; i < 3; ++i)
  {
    for (unsigned j=0; j < 3; ++j)
    {
      b[i] += a[3*j + i] * x[j];
    }
  }

  return b;
}

// computes a compound rotation matrix such that the z axis is rotated to align with the
// given vector.
// This type of transformation is not unique, the particular transformation done here is to
// rotate about the X axis, then about the Y axis
std::array<double, 9> compute_rotation_matrix(const std::array<double, 3>& new_z);

class ChangeOfBasis
{
  public:
    ChangeOfBasis(const std::array<double, 9>& basis) :
      m_basis(basis)
    {}

    std::array<double, 3> change_basis_forward(const std::array<double, 3>& pt) const
    {
      return matvec_transposed(m_basis, pt);
    }

    std::array<double, 3> change_basis_reverse(const std::array<double, 3>& pt) const
    {
      return matvec(m_basis, pt);
    }

  private:
    std::array<double, 9> m_basis;  // each column is a basis vector
};

}

#endif