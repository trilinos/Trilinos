// @HEADER
// ************************************************************************
//
//                    Intrepid MiniTensor Subpackage
//                 Copyright (2013) Sandia Corporation
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

#if !defined(Intrepid_MiniTensor_Tensor_t_cc)
#define Intrepid_MiniTensor_Tensor_t_cc

namespace Intrepid {

  //
  // R^N tensor input
  // \param A tensor
  // \param is input stream
  // \return is input stream
  //
  template<typename T>
  std::istream &
  operator>>(std::istream & is, Tensor<T> & A)
  {

    Index const
    N = A.get_dimension();

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        is >> A(i,j);
      }
    }

    return is;
  }

  //
  // R^N tensor output
  // \param A tensor
  // \param os output stream
  // \return os output stream
  //
  template<typename T>
  std::ostream &
  operator<<(std::ostream & os, Tensor<T> const & A)
  {
    Index const
    N = A.get_dimension();

    if (N == 0) {
      return os;
    }

    for (Index i = 0; i < N; ++i) {

      os << std::scientific << A(i,0);

      for (Index j = 1; j < N; ++j) {
        os << std::scientific << "," << A(i,j);
      }

      os << std::endl;
    }

    return os;
  }

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor_t_cc
