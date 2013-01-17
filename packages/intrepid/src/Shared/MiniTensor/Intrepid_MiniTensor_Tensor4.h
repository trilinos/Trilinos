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

#if !defined(Intrepid_MiniTensor_Tensor4_h)
#define Intrepid_MiniTensor_Tensor4_h

#include "Intrepid_MiniTensor_Tensor3.h"

namespace Intrepid {

  ///
  /// Fourth order tensor in R^N.
  ///
  template<typename T>
  class Tensor4
  {
  public:

    ///
    /// Default constructor
    ///
    Tensor4();

    ///
    /// 4th-order tensor constructor with NaNs
    ///
    explicit
    Tensor4(Index const N);

    ///
    /// 4th-order tensor constructor with a scalar
    /// \param s all components set to this scalar
    ///
    Tensor4(Index const N, T const & s);

    ///
    /// Copy constructor
    /// 4th-order tensor constructor with 4th-order tensor
    /// \param A from which components are copied
    ///
    Tensor4(Tensor4<T> const & A);

    ///
    /// 4th-order tensor simple destructor
    ///
    ~Tensor4();

    ///
    /// Indexing for constant 4th-order tensor
    /// \param i index
    /// \param j index
    /// \param k index
    /// \param l index
    ///
    T const &
    operator()(
        Index const i,
        Index const j,
        Index const k,
        Index const l) const;

    ///
    /// 4th-order tensor indexing
    /// \param i index
    /// \param j index
    /// \param k index
    /// \param l index
    ///
    T &
    operator()(
        Index const i,
        Index const j,
        Index const k,
        Index const l);

    ///
    /// \return dimension
    ///
    Index
    get_dimension() const;

    ///
    /// \param N dimension of 4th-order tensor
    ///
    void
    set_dimension(Index const N);

    ///
    /// 4th-order tensor copy assignment
    ///
    Tensor4<T> &
    operator=(Tensor4<T> const & A);

    ///
    /// 4th-order tensor increment
    /// \param A added to this tensor
    ///
    Tensor4<T> &
    operator+=(Tensor4<T> const & A);

    ///
    /// 4th-order tensor decrement
    /// \param A substracted from this tensor
    ///
    Tensor4<T> &
    operator-=(Tensor4<T> const & A);

    ///
    /// Fill 4th-order tensor with zeros
    ///
    void
    clear();

  private:

    ///
    /// Tensor dimension
    ///
    Index
    dimension;

    ///
    /// Tensor components
    ///
    T *
    e;

  };

  ///
  /// 4th-order tensor addition
  /// \param A 4th-order tensor
  /// \param B 4th-order tensor
  /// \return \f$ A + B \f$
  ///
  template<typename T>
  Tensor4<T>
  operator+(Tensor4<T> const & A, Tensor4<T> const & B);

  ///
  /// 4th-order tensor substraction
  /// \param A 4th-order tensor
  /// \param B 4th-order tensor
  /// \return \f$ A - B \f$
  ///
  template<typename T>
  Tensor4<T>
  operator-(Tensor4<T> const & A, Tensor4<T> const & B);

  ///
  /// 4th-order tensor minus
  /// \return \f$ -A \f$
  ///
  template<typename T>
  Tensor4<T>
  operator-(Tensor4<T> const & A);

  ///
  /// 4th-order equality
  /// Tested by components
  ///
  template<typename T>
  bool
  operator==(Tensor4<T> const & A, Tensor4<T> const & B);

  ///
  /// 4th-order inequality
  /// Tested by components
  ///
  template<typename T>
  bool
  operator!=(Tensor4<T> const & A, Tensor4<T> const & B);

  ///
  /// Scalar 4th-order tensor product
  /// \param s scalar
  /// \param A 4th-order tensor
  /// \return \f$ s A \f$
  ///
  template<typename T, typename S>
  Tensor4<T>
  operator*(S const & s, Tensor4<T> const & A);

  ///
  /// 4th-order tensor scalar product
  /// \param A 4th-order tensor
  /// \param s scalar
  /// \return \f$ s A \f$
  ///
  template<typename T, typename S>
  Tensor4<T>
  operator*(Tensor4<T> const & A, S const & s);


  /// Tensor4 Tensor4 double dot product
  /// \param A Tensor4
  /// \param B Tensor4
  /// \return a Tensor4 \f$ C_{ijkl} = A_{ijmn} : B){mnkl} \f$
  template<typename T>
  Tensor4<T>
  dotdot(Tensor4<T> const & A, Tensor4<T> const & B);

  ///
  /// 4th-order tensor transpose
  ///
  template<typename T>
  Tensor4<T>
  transpose(Tensor4<T> const & A);

  ///
  /// 4th-order identity I1
  /// \return \f$ \delta_{ik} \delta_{jl} \f$ such that \f$ A = I_1 A \f$
  ///
  template<typename T>
  const Tensor4<T>
  identity_1(Index const N);

  ///
  /// 4th-order identity I2
  /// \return \f$ \delta_{il} \delta_{jk} \f$ such that \f$ A^T = I_2 A \f$
  ///
  template<typename T>
  const Tensor4<T>
  identity_2(Index const N);

  ///
  /// 4th-order identity I3
  /// \return \f$ \delta_{ij} \delta_{kl} \f$ such that \f$ I_A I = I_3 A \f$
  ///
  template<typename T>
  const Tensor4<T>
  identity_3(Index const N);

  ///
  /// 4th-order tensor vector dot product
  /// \param A 4th-order tensor
  /// \param u vector
  /// \return 3rd-order tensor \f$ A dot u \f$ as \f$ B_{ijk}=A_{ijkl}u_{l} \f$
  ///
  template<typename T>
  Tensor3<T>
  dot(Tensor4<T> const & A, Vector<T> const & u);

  ///
  /// vector 4th-order tensor dot product
  /// \param A 4th-order tensor
  /// \param u vector
  /// \return 3rd-order tensor \f$ u dot A \f$ as \f$ B_{jkl}=u_{i} A_{ijkl} \f$
  ///
  template<typename T>
  Tensor3<T>
  dot(Vector<T> const & u, Tensor4<T> const & A);

  ///
  /// 4th-order tensor vector dot2 product
  /// \param A 4th-order tensor
  /// \param u vector
  /// \return 3rd-order tensor \f$ A dot2 u \f$ as \f$ B_{ijl}=A_{ijkl}u_{k} \f$
  ///
  template<typename T>
  Tensor3<T>
  dot2(Tensor4<T> const & A, Vector<T> const & u);

  ///
  /// vector 4th-order tensor dot2 product
  /// \param A 4th-order tensor
  /// \param u vector
  /// \return 3rd-order tensor \f$ u dot2 A \f$ as \f$ B_{ikl}=u_{j}A_{ijkl} \f$
  ///
  template<typename T>
  Tensor3<T>
  dot2(Vector<T> const & u, Tensor4<T> const & A);

  ///
  /// 4th-order tensor 2nd-order tensor double dot product
  /// \param A 4th-order tensor
  /// \param B 2nd-order tensor
  /// \return 2nd-order tensor \f$ A:B \f$ as \f$ C_{ij}=A_{ijkl}B_{kl} \f$
  ///
  template<typename T>
  Tensor<T>
  dotdot(Tensor4<T> const & A, Tensor<T> const & B);

  ///
  /// 2nd-order tensor 4th-order tensor double dot product
  /// \param B 2nd-order tensor
  /// \param A 4th-order tensor
  /// \return 2nd-order tensor \f$ B:A \f$ as \f$ C_{kl}=A_{ijkl}B_{ij} \f$
  ///
  template<typename T>
  Tensor<T>
  dotdot(Tensor<T> const & B, Tensor4<T> const & A);

  ///
  /// 2nd-order tensor 2nd-order tensor tensor product
  /// \param A 2nd-order tensor
  /// \param B 2nd-order tensor
  /// \return \f$ A \otimes B \f$
  ///
  template<typename T>
  Tensor4<T>
  tensor(Tensor<T> const & A, Tensor<T> const & B);

  ///
  /// odot operator useful for \f$ \frac{\partial A^{-1}}{\partial A} \f$
  /// see Holzapfel eqn 6.165
  /// \param A 2nd-order tensor
  /// \param B 2nd-order tensor
  /// \return \f$ A \odot B \f$ which is
  /// \f$ C_{ijkl} = \frac{1}{2}(A_{ik} B_{jl} + A_{il} B_{jk}) \f$
  ///
  template<typename T>
  Tensor4<T>
  odot(Tensor<T> const & A, Tensor<T> const & B);

  ///
  /// 4th-order input
  /// \param A 4th-order tensor
  /// \param is input stream
  /// \return is input stream
  ///
  template<typename T>
  std::istream &
  operator>>(std::istream & is, Tensor4<T> & A);

  ///
  /// 4th-order output
  /// \param A 4th-order tensor
  /// \param os output stream
  /// \return os output stream
  ///
  template<typename T>
  std::ostream &
  operator<<(std::ostream & os, Tensor4<T> const & A);

} // namespace Intrepid

#include "Intrepid_MiniTensor_Tensor4.i.cc"
#include "Intrepid_MiniTensor_Tensor4.t.cc"

#endif //Intrepid_MiniTensor_Tensor4_h
