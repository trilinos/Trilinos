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

#include "Intrepid_ConfigDefs.hpp"
#include "Sacado.hpp"

namespace Intrepid {

/// Static assertion
#define STATIC_ASSERT(condition, name)\
  typedef char static_assertion_failed_ ## name [(condition) ? 1 : -1]

/// Indexing type
typedef unsigned int Index;

/// High count type
typedef long unsigned int LongCount;

/// Floating point type
typedef double Real;

/// Complex type
typedef std::complex<Real> Complex;

/// The classes
template <typename T, Index N> class Vector;
template <typename T, Index N> class Tensor;
template <typename T, Index N> class Tensor3;
template <typename T, Index N> class Tensor4;

/// Indicator for dynamic storage
Index const
DYNAMIC = 0;

/// For use with type promotion
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

/// For Sacado traits

using std::string;

template <Index N>
struct dimension_string {
  static string eval() {return string("INVALID");}
};

template <>
struct dimension_string<DYNAMIC> {
  static string eval() {return string("DYNAMIC");}
};

template <>
struct dimension_string<1> {
  static string eval() {return string("1");}
};

template <>
struct dimension_string<2> {
  static string eval() {return string("2");}
};

template <>
struct dimension_string<3> {
  static string eval() {return string("3");}
};

template <>
struct dimension_string<4> {
  static string eval() {return string("4");}
};

} // namespace Intrepid

namespace Sacado {

using Intrepid::DYNAMIC;
using Intrepid::Index;
using Intrepid::Vector;
using Intrepid::Tensor;
using Intrepid::Tensor3;
using Intrepid::Tensor4;
using Intrepid::dimension_string;
using std::complex;
using std::string;

/// Specialization of Promote for Index
template <>
struct Promote<double, Index> {
  typedef double type;
};

template <>
struct Promote<Index, double> {
  typedef double type;
};

template <>
struct Promote<float, Index> {
  typedef float type;
};

template <>
struct Promote<Index, float> {
  typedef float type;
};

template <>
struct Promote<complex<double>, Index> {
  typedef complex<double> type;
};

template <>
struct Promote<Index, complex<double> > {
  typedef complex<double> type;
};

template <>
struct Promote<complex<float>, Index> {
  typedef complex<float> type;
};

template <>
struct Promote<Index, complex<float> > {
  typedef complex<float> type;
};

/// Sacado traits specializations for Vector
template <typename T, Index N>
struct ScalarType< Vector<T, N> > {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType< Vector<T, N> > {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType< Vector<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType< Vector<T, N> > {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value< Vector<T, N> > {
  typedef typename ValueType< Vector<T, N> >::type value_type;
  static const Vector<value_type, N>
  eval(Vector<T, N> const & x)
  {
    Vector<value_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N>
struct ScalarValue< Vector<T, N> > {
  typedef typename ScalarType< Vector<T, N> >::type scalar_type;
  static const Vector<scalar_type, N>
  eval(Vector<T, N> const & x)
  {
    Vector<scalar_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }
    return v;
  }
};

template <typename T, Index N>
struct StringName< Vector<T, N> > {
  static string
  eval()
  {
    return string("Vector<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual< Vector<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized< Vector<T, N> > {
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Vector<T, DYNAMIC> >
{
  static const bool value = false;
};

/// Sacado traits specializations for Tensor
template <typename T, Index N>
struct ScalarType< Tensor<T, N> > {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType< Tensor<T, N> > {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType< Tensor<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType< Tensor<T, N> > {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value< Tensor<T, N> > {
  typedef typename ValueType< Tensor<T, N> >::type value_type;
  static const Tensor<value_type, N>
  eval(Tensor<T, N> const & x)
  {
    Tensor<value_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N>
struct ScalarValue< Tensor<T, N> > {
  typedef typename ScalarType< Tensor<T, N> >::type scalar_type;
  static const Tensor<scalar_type, N>
  eval(Tensor<T, N> const & x)
  {
    Tensor<scalar_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N>
struct StringName< Tensor<T, N> > {
  static string
  eval()
  {
    return string("Tensor<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual< Tensor<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized< Tensor<T, N> > {
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Tensor<T, DYNAMIC> >
{
  static const bool value = false;
};

/// Sacado traits specializations for Tensor3
template <typename T, Index N>
struct ScalarType< Tensor3<T, N> > {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType< Tensor3<T, N> > {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType< Tensor3<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType< Tensor3<T, N> > {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value< Tensor3<T, N> > {
  typedef typename ValueType< Tensor3<T, N> >::type value_type;
  static const Tensor3<value_type, N>
  eval(Tensor3<T, N> const & x)
  {
    Tensor3<value_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N>
struct ScalarValue< Tensor3<T, N> > {
  typedef typename ScalarType< Tensor3<T, N> >::type scalar_type;
  static const Tensor3<scalar_type, N>
  eval(Tensor3<T, N> const & x)
  {
    Tensor3<scalar_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N>
struct StringName< Tensor3<T, N> > {
  static string
  eval()
  {
    return string("Tensor3<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual< Tensor3<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized< Tensor3<T, N> >
{
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Tensor3<T, DYNAMIC> >
{
  static const bool value = false;
};

/// Sacado traits specializations for Tensor4
template <typename T, Index N>
struct ScalarType< Tensor4<T, N> > {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType< Tensor4<T, N> > {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType< Tensor4<T, N> > {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType< Tensor4<T, N> > {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value< Tensor4<T, N> > {
  typedef typename ValueType< Tensor4<T, N> >::type value_type;
  static const Tensor4<value_type, N>
  eval(Tensor4<T, N> const & x)
  {
    Tensor4<value_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N>
struct ScalarValue< Tensor4<T, N> > {
  typedef typename ScalarType< Tensor4<T, N> >::type scalar_type;
  static const Tensor4<scalar_type, N>
  eval(Tensor4<T, N> const & x)
  {
    Tensor4<scalar_type, N> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N>
struct StringName< Tensor4<T, N> > {
  static string
  eval()
  {
    return string("Tensor4<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual< Tensor4<T, N> > {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized< Tensor4<T, N> >
{
  static const bool value = true;
};

template <typename T>
struct IsStaticallySized< Tensor4<T, DYNAMIC> >
{
  static const bool value = false;
};

} // namespace Sacado

#endif // Intrepid_MiniTensor_Definitions_h
