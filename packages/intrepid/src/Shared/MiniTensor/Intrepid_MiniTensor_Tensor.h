// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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

#if !defined(Intrepid_MiniTensor_Tensor_h)
#define Intrepid_MiniTensor_Tensor_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include <boost/tuple/tuple.hpp>

#include "Intrepid_MiniTensor_Vector.h"

namespace Intrepid {

  ///
  /// Second order tensor in R^N.
  ///
  template<typename T>
  class Tensor
  {
  public:

    ///
    /// Component type
    ///
    typedef T type;

    ///
    /// Default constructor
    ///
    Tensor();

    ///
    /// Constructor that initializes to NaNs
    /// \param N dimension
    ///
    explicit
    Tensor(Index const N);

    ///
    /// Create tensor from a scalar
    /// \param N dimension
    /// \param s all components are set equal to this value
    ///
    Tensor(Index const N, T const & s);

    ///
    /// Create tensor specifying components
    /// \param N dimension
    /// \param  s00, s01, ... components in the R^2 canonical basis
    ///
    Tensor(T const & s00, T const & s01, T const & s10, T const & s11);

    ///
    /// Create tensor specifying components
    /// \param N dimension
    /// \param  s00, s01, ... components in the R^3 canonical basis
    ///
    Tensor(
        T const & s00, T const & s01, T const & s02,
        T const & s10, T const & s11, T const & s12,
        T const & s20, T const & s21, T const & s22);

    ///
    /// Create tensor from array - const version
    /// \param data_ptr pointer into the array
    ///
    Tensor(Index const N, T const * data_ptr);

    ///
    /// Copy constructor
    /// \param A the values of its components are copied to the new tensor
    ///
    Tensor(Tensor<T> const & A);

    ///
    /// Simple destructor
    ///
    ~Tensor();

    ///
    /// Indexing for constant tensor
    /// \param i index
    /// \param j index
    ///
    T const &
    operator()(Index const i, Index const j) const;

    ///
    /// Tensor indexing
    /// \param i index
    /// \param j index
    ///
    T &
    operator()(Index const i, Index const j);

    ///
    /// \return dimension
    ///
    Index
    get_dimension() const;

    ///
    /// \param N dimension of 2nd-order tensor
    ///
    void
    set_dimension(Index const N);

    ///
    /// Fill components from array defined by pointer.
    /// \param data_ptr pointer into array for filling components
    ///
    void
    fill(T const * data_ptr);

    ///
    /// Copy assignment
    /// \param A the values of its components are copied to this tensor
    ///
    Tensor<T> &
    operator=(Tensor<T> const & A);

    ///
    /// Tensor increment
    /// \param A added to current tensor
    ///
    Tensor<T> &
    operator+=(Tensor<T> const & A);

    ///
    /// Tensor decrement
    /// \param A substracted from current tensor
    ///
    Tensor<T> &
    operator-=(Tensor<T> const & A);

    ///
    /// Fill with zeros
    ///
    void
    clear();

    ///
    /// Tensor order
    ///
    static
    Index
    order() {return 2U;};

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
  /// Tensor addition
  /// \return \f$ A + B \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  operator+(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Tensor substraction
  /// \return \f$ A - B \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  operator-(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Tensor minus
  /// \return \f$ -A \f$
  ///
  template<typename T>
  Tensor<T>
  operator-(Tensor<T> const & A);

  ///
  /// Tensor equality
  /// Tested by components
  /// \return \f$ A \equiv B \f$
  ///
  template<typename T>
  bool
  operator==(Tensor<T> const & A, Tensor<T> const & B);

  ///
  /// Tensor inequality
  /// Tested by components
  /// \return \f$ A \neq B \f$
  ///
  template<typename T>
  bool
  operator!=(Tensor<T> const & A, Tensor<T> const & B);

  ///
  /// Tensor vector product v = A u
  /// \param A tensor
  /// \param u vector
  /// \return \f$ A u \f$
  ///
  template<typename S, typename T>
  Vector<typename Promote<S, T>::type>
  operator*(Tensor<T> const & A, Vector<S> const & u);

  ///
  /// Vector tensor product v = u A
  /// \param A tensor
  /// \param u vector
  /// \return \f$ u A = A^T u \f$
  ///
  template<typename S, typename T>
  Vector<typename Promote<S, T>::type>
  operator*(Vector<S> const & u, Tensor<T> const & A);

  ///
  /// Tensor dot product C = A B
  /// \return \f$ A \cdot B \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  operator*(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Scalar tensor product
  /// \param s scalar
  /// \param A tensor
  /// \return \f$ s A \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  operator*(S const & s, Tensor<T> const & A);

  ///
  /// Tensor scalar product
  /// \param A tensor
  /// \param s scalar
  /// \return \f$ s A \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  operator*(Tensor<T> const & A, S const & s);

  ///
  /// Tensor scalar division
  /// \param A tensor
  /// \param s scalar
  /// \return \f$ A / s \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  operator/(Tensor<T> const & A, S const & s);

  ///
  /// Tensor input
  /// \param A tensor
  /// \param is input stream
  /// \return is input stream
  ///
  template<typename T>
  std::istream &
  operator>>(std::istream & is, Tensor<T> & A);

  ///
  /// Tensor output
  /// \param A tensor
  /// \param os output stream
  /// \return os output stream
  ///
  template<typename T>
  std::ostream &
  operator<<(std::ostream & os, Tensor<T> const & A);

  ///
  /// Tensor vector product v = A u
  /// \param A tensor
  /// \param u vector
  /// \return \f$ A u \f$
  ///
  template<typename S, typename T>
  Vector<typename Promote<S, T>::type>
  dot(Tensor<T> const & A, Vector<S> const & u);

  ///
  /// Vector tensor product v = u A
  /// \param A tensor
  /// \param u vector
  /// \return \f$ u A = A^T u \f$
  ///
  template<typename S, typename T>
  Vector<typename Promote<S, T>::type>
  dot(Vector<S> const & u, Tensor<T> const & A);

  ///
  /// Tensor tensor product C = A B
  /// \param A tensor
  /// \param B tensor
  /// \return a tensor \f$ A \cdot B \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dot(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Tensor tensor product C = A^T B
  /// \param A tensor
  /// \param B tensor
  /// \return a tensor \f$ A^T \cdot B \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  t_dot(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Tensor tensor product C = A B^T
  /// \param A tensor
  /// \param B tensor
  /// \return a tensor \f$ A \cdot B^T \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dot_t(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Tensor tensor product C = A^T B^T
  /// \param A tensor
  /// \param B tensor
  /// \return a tensor \f$ A^T \cdot B^T \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  t_dot_t(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Tensor tensor double dot product (contraction)
  /// \param A tensor
  /// \param B tensor
  /// \return a scalar \f$ A : B \f$
  ///
  template<typename S, typename T>
  typename Promote<S, T>::type
  dotdot(Tensor<S> const & A, Tensor<T> const & B);

  ///
  /// Dyad
  /// \param u vector
  /// \param v vector
  /// \return \f$ u \otimes v \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dyad(Vector<S> const & u, Vector<T> const & v);

  ///
  /// Bun operator, just for Jay
  /// \param u vector
  /// \param v vector
  /// \return \f$ u \otimes v \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  bun(Vector<S> const & u, Vector<T> const & v);

  ///
  /// Tensor product
  /// \param u vector
  /// \param v vector
  /// \return \f$ u \otimes v \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  tensor(Vector<S> const & u, Vector<T> const & v);

  ///
  /// Diagonal tensor from vector
  /// \param v vector
  /// \return A = diag(v)
  ///
  template<typename T>
  Tensor<T>
  diag(Vector<T> const & v);

  ///
  /// Diagonal of tensor in a vector
  /// \param A tensor
  /// \return v = diag(A)
  ///
  template<typename T>
  Vector<T>
  diag(Tensor<T> const & A);

  ///
  /// Zero 2nd-order tensor
  /// All components are zero
  ///
  template<typename T>
  const Tensor<T>
  zero(Index const N);

  ///
  /// 2nd-order identity tensor
  ///
  template<typename T>
  const Tensor<T>
  identity(Index const N);

  ///
  /// 2nd-order identity tensor, à la Matlab
  ///
  template<typename T>
  const Tensor<T>
  eye(Index const N);

  ///
  /// R^N 2nd-order tensor transpose
  ///
  template<typename T>
  Tensor<T>
  transpose(Tensor<T> const & A);

  ///
  /// C^N 2nd-order tensor adjoint
  ///
  template<typename T>
  Tensor<T>
  adjoint(Tensor<T> const & A);

  ///
  /// Symmetric part of 2nd-order tensor
  /// \return \f$ \frac{1}{2}(A + A^T) \f$
  ///
  template<typename T>
  Tensor<T>
  sym(Tensor<T> const & A);

  ///
  /// Skew symmetric part of 2nd-order tensor
  /// \return \f$ \frac{1}{2}(A - A^T) \f$
  ///
  template<typename T>
  Tensor<T>
  skew(Tensor<T> const & A);

  ///
  /// Skew symmetric 2nd-order tensor from vector valid for R^3 only.
  /// R^N with N != 3 will produce an error
  /// \param u vector
  /// \return \f$ {{0, -u_2, u_1}, {u_2, 0, -u_0}, {-u_1, u+0, 0}} \f$
  ///
  template<typename T>
  Tensor<T>
  skew(Vector<T> const & u);

} // namespace Intrepid

#include "Intrepid_MiniTensor_Tensor.i.h"
#include "Intrepid_MiniTensor_Tensor.t.h"

#endif //Intrepid_MiniTensor_Tensor_h
