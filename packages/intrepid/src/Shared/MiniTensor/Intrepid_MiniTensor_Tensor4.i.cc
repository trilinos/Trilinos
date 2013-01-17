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

#if !defined(Intrepid_MiniTensor_Tensor4_i_cc)
#define Intrepid_MiniTensor_Tensor4_i_cc

namespace Intrepid
{

  //
  // R^N 4th-order tensor transpose
  // per Holzapfel 1.157
  //
  template<typename T>
  inline Tensor4<T> transpose(Tensor4<T> const & A)
  {
    Index const N = A.get_dimension();

    Tensor4<T> B(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            B(i, j, k, l) = A(k, l, i, j);
          }
        }
      }
    }

    return B;
  }

  //
  // Dimension
  // \return dimension
  //
  template<typename T>
  inline
  Index
  Tensor4<T>::get_dimension() const
  {
    return dimension;
  }

  //
  // R^N indexing for constant 4th order tensor
  // \param i index
  // \param j index
  // \param k index
  // \param l index
  //
  template<typename T>
  inline
  T const &
  Tensor4<T>::operator()(
      Index const i, Index const j, Index const k, Index const l) const
  {
    Index const
    N = get_dimension();

    assert(i < N);
    assert(j < N);
    assert(k < N);
    assert(l < N);

    return e[((i * N + j) * N + k) * N + l];
  }

  //
  // R^N 4th-order tensor indexing
  // \param i index
  // \param j index
  // \param k index
  // \param l index
  //
  template<typename T>
  inline
  T &
  Tensor4<T>::operator()(
      Index const i, Index const j, Index const k, Index const l)
  {
    Index const
    N = get_dimension();

    assert(i < N);
    assert(j < N);
    assert(k < N);
    assert(l < N);

    return e[((i * N + j) * N + k) * N + l];
  }

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor4_i_cc
