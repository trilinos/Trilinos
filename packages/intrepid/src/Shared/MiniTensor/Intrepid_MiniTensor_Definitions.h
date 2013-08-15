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

#if !defined(Intrepid_MiniTensor_Definitions_h)
#define Intrepid_MiniTensor_Definitions_h

#include <complex>

#include "Sacado.hpp"
#include "Sacado_mpl_disable_if.hpp"

namespace Intrepid {

///
/// Indexing type
///
typedef unsigned int Index;

///
/// High count type
///
typedef long unsigned int LongCount;

///
/// Floating point type
///
typedef double Real;

///
/// Complex type
///
typedef std::complex<Real> Complex;

///
/// The classes
///
template <typename T, Index N> class Vector;
template <typename T, Index N> class Tensor;
template <typename T, Index N> class Tensor3;
template <typename T, Index N> class Tensor4;

///
/// Indicator for dynamic storage
///
Index const
DYNAMIC = 0;

///
/// For use with type promotion
///
using Sacado::Promote;
using Sacado::mpl::lazy_disable_if;
using Sacado::mpl::disable_if_c;

/// Vector
template <typename T>
struct is_vector {
  static const bool value = false;
};

template <typename T, Index N>
struct is_vector< Vector<T, N> > {
  static const bool value = true;
};

template <typename T, Index N>
struct apply_vector {
  typedef Vector<typename T::type, N> type;
};

/// 2nd-order tensor
template <typename T>
struct is_tensor {
  static const bool value = false;
};

template <typename T, Index N>
struct is_tensor< Tensor<T, N> > {
  static const bool value = true;
};

template <typename T, Index N>
struct apply_tensor {
  typedef Tensor<typename T::type, N> type;
};

/// 3rd-order tensor
template <typename T>
struct is_tensor3 {
  static const bool value = false;
};

template <typename T, Index N>
struct is_tensor3< Tensor3<T, N> > {
  static const bool value = true;
};

template <typename T, Index N>
struct apply_tensor3 {
  typedef Tensor3<typename T::type, N> type;
};

/// 4th-order tensor
template <typename T>
struct is_tensor4 {
  static const bool value = false;
};

template <typename T, Index N>
struct is_tensor4< Tensor4<T, N> > {
  static const bool value = true;
};

template <typename T, Index N>
struct apply_tensor4 {
  typedef Tensor4<typename T::type, N> type;
};

/// Tensors from 1st to 4th order
template <typename T>
struct order_1234 {
  static const bool value = false;
};

template <typename T, Index N>
struct order_1234< Vector<T, N> > {
  static const bool value = true;
};

template <typename T, Index N>
struct order_1234< Tensor<T, N> > {
  static const bool value = true;
};

template <typename T, Index N>
struct order_1234< Tensor3<T, N> > {
  static const bool value = true;
};

template <typename T, Index N>
struct order_1234< Tensor4<T, N> > {
  static const bool value = true;
};

} // namespace Intrepid

namespace Sacado {

///
/// Specialization of Promote for Intrepid::Index
///
template <>
struct Promote<double, Intrepid::Index> {
  typedef double type;
};

template <>
struct Promote<Intrepid::Index, double> {
  typedef double type;
};

template <>
struct Promote<float, Intrepid::Index> {
  typedef float type;
};

template <>
struct Promote<Intrepid::Index, float> {
  typedef float type;
};

template <>
struct Promote<std::complex<double>, Intrepid::Index> {
  typedef std::complex<double> type;
};

template <>
struct Promote<Intrepid::Index, std::complex<double> > {
  typedef std::complex<double> type;
};

template <>
struct Promote<std::complex<float>, Intrepid::Index> {
  typedef std::complex<float> type;
};

template <>
struct Promote<Intrepid::Index, std::complex<float> > {
  typedef std::complex<float> type;
};

///
/// Sacado traits specializations for Vector
///
template <typename T, Intrepid::Index N>
struct ScalarType< Intrepid::Vector<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct ValueType< Intrepid::Vector<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct IsADType< Intrepid::Vector<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Intrepid::Index N>
struct IsScalarType< Intrepid::Vector<T, N> > {
  static bool const value = false;
};

template <typename T, Intrepid::Index N>
struct Value< Intrepid::Vector<T, N> > {
  typedef typename ValueType<T>::type value_type;
  static const Intrepid::Vector<value_type, N> &
  eval(Intrepid::Vector<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Vector<value_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      v(i) = x(i).val();
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct ScalarValue< Intrepid::Vector<T, N> > {
  typedef typename ValueType<T>::type value_type;
  typedef typename ScalarType<T>::type scalar_type;
  static const Intrepid::Vector<scalar_type, N> &
  eval(Intrepid::Vector<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Vector<scalar_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      v(i) = ScalarValue<value_type>::eval(x(i).val());
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct StringName< Intrepid::Vector<T, N> > {
  static std::string
  eval()
  {
    return std::string("Intrepid::Vector< ") +
        StringName<T>::eval() + std::string(" >");
  }
};

template <typename T, Intrepid::Index N>
struct IsEqual< Intrepid::Vector<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Intrepid::Index N>
struct IsStaticallySized< Intrepid::Vector<T, N> > {
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Intrepid::Vector<T, Intrepid::DYNAMIC> >
{
  static const bool value = false;
};

///
/// Sacado traits specializations for Tensor
///
template <typename T, Intrepid::Index N>
struct ScalarType< Intrepid::Tensor<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct ValueType< Intrepid::Tensor<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct IsADType< Intrepid::Tensor<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Intrepid::Index N>
struct IsScalarType< Intrepid::Tensor<T, N> > {
  static bool const value = false;
};

template <typename T, Intrepid::Index N>
struct Value< Intrepid::Tensor<T, N> > {
  typedef typename ValueType<T>::type value_type;
  static const Intrepid::Tensor<value_type, N> &
  eval(Intrepid::Tensor<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Tensor<value_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      for (Intrepid::Index j = 0; j < dimension; ++j) {
        v(i,j) = x(i,j).val();
      }
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct ScalarValue< Intrepid::Tensor<T, N> > {
  typedef typename ValueType<T>::type value_type;
  typedef typename ScalarType<T>::type scalar_type;
  static const Intrepid::Tensor<scalar_type, N> &
  eval(Intrepid::Tensor<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Tensor<scalar_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      for (Intrepid::Index j = 0; j < dimension; ++j) {
        v(i,j) = ScalarValue<value_type>::eval(x(i,j).val());
      }
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct StringName< Intrepid::Tensor<T, N> > {
  static std::string
  eval()
  {
    return std::string("Intrepid::Tensor< ") +
        StringName<T>::eval() + std::string(" >");
  }
};

template <typename T, Intrepid::Index N>
struct IsEqual< Intrepid::Tensor<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Intrepid::Index N>
struct IsStaticallySized< Intrepid::Tensor<T, N> > {
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Intrepid::Tensor<T, Intrepid::DYNAMIC> >
{
  static const bool value = false;
};

///
/// Sacado traits specializations for Tensor3
///
template <typename T, Intrepid::Index N>
struct ScalarType< Intrepid::Tensor3<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct ValueType< Intrepid::Tensor3<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct IsADType< Intrepid::Tensor3<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Intrepid::Index N>
struct IsScalarType< Intrepid::Tensor3<T, N> > {
  static bool const value = false;
};

template <typename T, Intrepid::Index N>
struct Value< Intrepid::Tensor3<T, N> > {
  typedef typename ValueType<T>::type value_type;
  static const Intrepid::Tensor3<value_type, N> &
  eval(Intrepid::Tensor3<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Tensor3<value_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      for (Intrepid::Index j = 0; j < dimension; ++j) {
        for (Intrepid::Index k = 0; k < dimension; ++k) {
          v(i,j,k) = x(i,j,k).val();
        }
      }
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct ScalarValue< Intrepid::Tensor3<T, N> > {
  typedef typename ValueType<T>::type value_type;
  typedef typename ScalarType<T>::type scalar_type;
  static const Intrepid::Tensor3<scalar_type, N> &
  eval(Intrepid::Tensor3<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Tensor3<scalar_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      for (Intrepid::Index j = 0; j < dimension; ++j) {
        for (Intrepid::Index k = 0; k < dimension; ++k) {
          v(i,j,k) = ScalarValue<value_type>::eval(x(i,j,k).val());
        }
      }
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct StringName< Intrepid::Tensor3<T, N> > {
  static std::string
  eval()
  {
    return std::string("Intrepid::Tensor3< ") +
        StringName<T>::eval() + std::string(" >");
  }
};

template <typename T, Intrepid::Index N>
struct IsEqual< Intrepid::Tensor3<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Intrepid::Index N>
struct IsStaticallySized< Intrepid::Tensor3<T, N> >
{
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Intrepid::Tensor3<T, Intrepid::DYNAMIC> >
{
  static const bool value = false;
};

///
/// Sacado traits specializations for Tensor4
///
template <typename T, Intrepid::Index N>
struct ScalarType< Intrepid::Tensor4<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct ValueType< Intrepid::Tensor4<T, N> > {
  typedef T type;
};

template <typename T, Intrepid::Index N>
struct IsADType< Intrepid::Tensor4<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Intrepid::Index N>
struct IsScalarType< Intrepid::Tensor4<T, N> > {
  static bool const value = false;
};

template <typename T, Intrepid::Index N>
struct Value< Intrepid::Tensor4<T, N> > {
  typedef typename ValueType<T>::type value_type;
  static const Intrepid::Tensor4<value_type, N> &
  eval(Intrepid::Tensor4<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Tensor4<value_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      for (Intrepid::Index j = 0; j < dimension; ++j) {
        for (Intrepid::Index k = 0; k < dimension; ++k) {
          for (Intrepid::Index l = 0; l < dimension; ++l) {
            v(i,j,k,l) = x(i,j,k,l).val();
          }
        }
      }
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct ScalarValue< Intrepid::Tensor4<T, N> > {
  typedef typename ValueType<T>::type value_type;
  typedef typename ScalarType<T>::type scalar_type;
  static const Intrepid::Tensor4<scalar_type, N> &
  eval(Intrepid::Tensor4<T, N> const & x)
  {
    Intrepid::Index const
    dimension = x.get_dimension();

    Intrepid::Tensor4<scalar_type, N>
    v(dimension);

    for (Intrepid::Index i = 0; i < dimension; ++i) {
      for (Intrepid::Index j = 0; j < dimension; ++j) {
        for (Intrepid::Index k = 0; k < dimension; ++k) {
          for (Intrepid::Index l = 0; l < dimension; ++l) {
            v(i,j,k,l) = ScalarValue<value_type>::eval(x(i,j,k,l).val());
          }
        }
      }
    }

    return v;
  }
};

template <typename T, Intrepid::Index N>
struct StringName< Intrepid::Tensor4<T, N> > {
  static std::string
  eval()
  {
    return std::string("Intrepid::Tensor4< ") +
        StringName<T>::eval() + std::string(" >");
  }
};

template <typename T, Intrepid::Index N>
struct IsEqual< Intrepid::Tensor4<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Intrepid::Index N>
struct IsStaticallySized< Intrepid::Tensor4<T, N> >
{
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Intrepid::Tensor4<T, Intrepid::DYNAMIC> >
{
  static const bool value = false;
};

} // namespace Sacado

#endif // Intrepid_MiniTensor_Definitions_h
