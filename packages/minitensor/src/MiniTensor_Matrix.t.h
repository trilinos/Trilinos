// @HEADER
// ************************************************************************
//
//                           MiniTensor Package
//                 Copyright (2016) Sandia Corporation
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(MiniTensor_Matrix_t_h)
#define MiniTensor_Matrix_t_h

namespace minitensor {

//
// Matrix input
//
template<typename T, Index M, Index N>
std::istream &
operator>>(std::istream & is, Matrix<T, M, N> & A)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  for (Index i = 0; i < num_rows; ++i) {
    for (Index j = 0; j < num_cols; ++j) {
      is >> A(i,j);
    }
  }

  return is;
}

//
// Matrix output
//
template<typename T, Index M, Index N>
std::ostream &
operator<<(std::ostream & os, Matrix<T, M, N> const & A)
{
  Index const
  num_rows = A.get_num_rows();

  Index const
  num_cols = A.get_num_cols();

  Index const
  dimension = num_rows * num_cols;

  if (dimension == 0) {
    return os;
  }

  os << std::scientific << std::setprecision(17);

  for (Index i = 0; i < num_rows; ++i) {

    os << std::setw(24) << A(i,0);

    for (Index j = 1; j < num_cols; ++j) {
      os << "," << std::setw(24) << A(i,j);
    }

    os << std::endl;
  }

  return os;
}

} // namespace minitensor

#endif // MiniTensor_Matrix_t_h
