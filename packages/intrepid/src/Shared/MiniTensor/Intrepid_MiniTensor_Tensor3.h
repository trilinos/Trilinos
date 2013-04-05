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

#if !defined(Intrepid_MiniTensor_Tensor3_h)
#define Intrepid_MiniTensor_Tensor3_h

#include "Intrepid_MiniTensor_Tensor.h"
#include "Teuchos_ArrayRCP.hpp"

namespace Intrepid {

  ///
  /// Third order tensor in R^N.
  ///
  template<typename T>
  class Tensor3
  {
  public:

    ///
    /// Component type
    ///
    typedef T type;

    ///
    /// Default constructor
    ///
    Tensor3();

    ///
    /// 3rd-order tensor constructor with NaNs
    ///
    explicit
    Tensor3(Index const N);

    ///
    /// 3rd-order tensor constructor with a scalar
    /// \param s all components set to this scalar
    ///
    Tensor3(Index const N, T const & s);

    ///
    /// Copy constructor
    /// 3rd-order tensor constructor from 3rd-order tensor
    /// \param A from which components are copied
    ///
    Tensor3(Tensor3<T> const & A);

    ///
    /// 3rd-order tensor simple destructor
    ///
    ~Tensor3();

    ///
    /// Indexing for constant 3rd-order tensor
    /// \param i index
    /// \param j index
    /// \param k index
    ///
    T const &
    operator()(Index const i, Index const j, Index const k) const;

    ///
    /// 3rd-order tensor indexing
    /// \param i index
    /// \param j index
    /// \param k index
    ///
    T &
    operator()(Index const i, Index const j, Index const k);

    ///
    /// \return dimension
    ///
    Index
    get_dimension() const;

    ///
    /// \param N dimension of 3rd-order tensor
    ///
    void
    set_dimension(Index const N);

    ///
    /// 3rd-order tensor copy assignment
    ///
    Tensor3<T> &
    operator=(Tensor3<T> const & A);

    ///
    /// 3rd-order tensor increment
    /// \param A added to this tensor
    ///
    Tensor3<T> &
    operator+=(Tensor3<T> const & A);

    ///
    /// 3rd-order tensor decrement
    /// \param A substracted from this tensor
    ///
    Tensor3<T> &
    operator-=(Tensor3<T> const & A);

    ///
    /// Fill 3rd-order tensor with zeros
    ///
    void
    clear();

    ///
    /// Tensor order
    ///
    static
    Index
    order() {return 3U;};

  private:

    ///
    /// Tensor dimension
    ///
    Index
    dimension;

    ///
    /// Tensor components
    ///
    Teuchos::ArrayRCP<T>
    e;

  };

  ///
  /// 3rd-order tensor addition
  /// \param A 3rd-order tensor
  /// \param B 3rd-order tensor
  /// \return \f$ A + B \f$
  ///
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  operator+(Tensor3<S> const & A, Tensor3<T> const & B);

  ///
  /// 3rd-order tensor substraction
  /// \param A 3rd-order tensor
  /// \param B 3rd-order tensor
  /// \return \f$ A - B \f$
  ///
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  operator-(Tensor3<S> const & A, Tensor3<T> const & B);

  ///
  /// 3rd-order tensor minus
  /// \return \f$ -A \f$
  ///
  template<typename T>
  Tensor3<T>
  operator-(Tensor3<T> const & A);

  ///
  /// 3rd-order tensor equality
  /// Tested by components
  ///
  template<typename T>
  bool
  operator==(Tensor3<T> const & A, Tensor3<T> const & B);

  ///
  /// 3rd-order tensor inequality
  /// Tested by components
  ///
  template<typename T>
  bool
  operator!=(Tensor3<T> const & A, Tensor3<T> const & B);

  ///
  /// Scalar 3rd-order tensor product
  /// \param s scalar
  /// \param A 3rd-order tensor
  /// \return \f$ s A \f$
  ///
  template<typename S, typename T>
  typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T> > >::type
  operator*(S const & s, Tensor3<T> const & A);

  ///
  /// 3rd-order tensor scalar product
  /// \param A 3rd-order tensor
  /// \param s scalar
  /// \return \f$ s A \f$
  ///
  template<typename S, typename T>
  typename lazy_disable_if< order_1234<S>, apply_tensor3< Promote<S,T> > >::type
  operator*(Tensor3<T> const & A, S const & s);

  ///
  /// 3rd-order tensor scalar division
  /// \param A 3rd-order tensor
  /// \param s scalar
  /// \return \f$ A / s \f$
  ///
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  operator/(Tensor3<T> const & A, S const & s);

  ///
  /// 3rd-order tensor vector product
  /// \param A 3rd-order tensor
  /// \param u vector
  /// \return \f$ C = A \cdot u := C_{ij} = A_{ijp} u_p \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dot(Tensor3<T> const & A, Vector<S> const & u);

  ///
  /// vector 3rd-order tensor product
  /// \param A 3rd-order tensor
  /// \param u vector
  /// \return \f$ C = u \cdot A := C_{ij} = u_p A{pij} \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dot(Vector<S> const & u, Tensor3<T> const & A);

  ///
  /// 3rd-order tensor vector product
  /// \param A 3rd-order tensor
  /// \param u vector
  /// \return \f$ C = A \cdot u := C_{ij} = A_{ipj} u_p \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dot2(Tensor3<T> const & A, Vector<S> const & u);

  ///
  /// vector 3rd-order tensor product
  /// \param u vector
  /// \param A 3rd-order tensor
  /// \return \f$ C = u \cdot A := C_{ij} = u_p A_{ipj} \f$
  ///
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dot2(Vector<S> const & u, Tensor3<T> const & A);

  ///
  /// 3rd-order tensor 2nd-order tensor product
  /// \param A 3rd-order tensor
  /// \param B 2nd-order tensor
  /// \return \f$ C = A \cdot B := C_{ijk} = A_{ijp} B_{pk} \f$
  ///
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot(Tensor3<T> const & A, Tensor<S> const & B);

  ///
  /// 2nd-order tensor 3rd-order tensor product
  /// \param A 2nd-order tensor
  /// \param B 3rd-order tensor
  /// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{pjk} \f$
  ///
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot(Tensor<S> const & A, Tensor3<T> const & B);

  ///
  /// 3rd-order tensor 2nd-order tensor product
  /// \param A 3rd-order tensor
  /// \param B 2nd-order tensor
  /// \return \f$ C = A \cdot B := C_{ijk} = A_{ipj} B_{pk} \f$
  ///
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot2(Tensor3<T> const & A, Tensor<S> const & B);

  ///
  /// 2nd-order tensor 3rd-order tensor product
  /// \param A 2nd-order tensor
  /// \param B 3rd-order tensor
  /// \return \f$ C = A \cdot B := C_{ijk} = A_{ip} B_{jpk} \f$
  ///
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot2(Tensor<S> const & A, Tensor3<T> const & B);

  ///
  /// 3rd-order tensor input
  /// \param A 3rd-order tensor
  /// \param is input stream
  /// \return is input stream
  ///
  template<typename T>
  std::istream &
  operator>>(std::istream & is, Tensor3<T> & A);

  ///
  /// 3rd-order tensor output
  /// \param A 3rd-order tensor
  /// \param os output stream
  /// \return os output stream
  ///
  template<typename T>
  std::ostream &
  operator<<(std::ostream & os, Tensor3<T> const & A);

} // namespace Intrepid

#include "Intrepid_MiniTensor_Tensor3.i.h"
#include "Intrepid_MiniTensor_Tensor3.t.h"

#endif //Intrepid_MiniTensor_Tensor3_h
