// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Definitions_h)
#define MiniTensor_Definitions_h

#include <complex>
#include <type_traits>

#include "MiniTensor_config.h"
#include "Sacado.hpp"

#if !defined( KOKKOS_INLINE_FUNCTION )
#define KOKKOS_INLINE_FUNCTION  inline
#endif

#if defined(KOKKOS_ENABLE_CUDA)
#define MT_ERROR_EXIT(...) \
  Kokkos::abort(#__VA_ARGS__)
#else
#define MT_ERROR_EXIT(...) \
  fprintf(stderr, "ERROR in: %s\n", __PRETTY_FUNCTION__); \
  fprintf(stderr, __VA_ARGS__); \
  fprintf(stderr, "\n"); \
  exit(1)
#endif // KOKKOS_ENABLE_CUDA

#if defined(KOKKOS_ENABLE_CUDA)
#define MT_WARNING(...) \
  Kokkos::abort(#__VA_ARGS__)
#else
#define MT_WARNING(...) \
  fprintf(stderr, "WARNING in: %s\n", __PRETTY_FUNCTION__); \
  fprintf(stderr, __VA_ARGS__); \
  fprintf(stderr, "\n")
#endif // KOKKOS_ENABLE_CUDA

namespace minitensor {

/// Indexing type
using Index = uint32_t;

constexpr Index
INDEX_SIZE{32};

/// High count type
using LongIndex = uint64_t;

constexpr Index
LONG_INDEX_SIZE{64};

/// Floating point type
using Real = double;

/// Complex type
using Complex = std::complex<Real>;

/// The classes
template <typename T, Index N> class Vector;
template <typename T, Index N> class Tensor;
template <typename T, Index N> class Tensor3;
template <typename T, Index N> class Tensor4;
template <typename T, Index M, Index N> class Matrix;

/// Indicator for dynamic storage
constexpr Index
DYNAMIC{0};

/// For use with type promotion
using Sacado::Promote;
using Sacado::mpl::lazy_disable_if;
using Sacado::mpl::disable_if_c;

/// Vector
template <typename T>
struct is_vector {
  static bool const value = false;
};

template <typename T, Index N>
struct is_vector<Vector<T, N>> {
  static bool const value = true;
};

template <typename T, Index N>
struct apply_vector {
  typedef Vector<typename T::type, N> type;
};

/// 2nd-order tensor
template <typename T>
struct is_tensor {
  static bool const value = false;
};

template <typename T, Index N>
struct is_tensor<Tensor<T, N>> {
  static bool const value = true;
};

template <typename T, Index N>
struct apply_tensor {
  typedef Tensor<typename T::type, N> type;
};

/// 3rd-order tensor
template <typename T>
struct is_tensor3 {
  static bool const value = false;
};

template <typename T, Index N>
struct is_tensor3<Tensor3<T, N>> {
  static bool const value = true;
};

template <typename T, Index N>
struct apply_tensor3 {
  typedef Tensor3<typename T::type, N> type;
};

/// 4th-order tensor
template <typename T>
struct is_tensor4 {
  static bool const value = false;
};

template <typename T, Index N>
struct is_tensor4<Tensor4<T, N>> {
  static bool const value = true;
};

template <typename T, Index N>
struct apply_tensor4 {
  typedef Tensor4<typename T::type, N> type;
};

/// Matrix
template <typename T>
struct is_matrix {
  static bool const value = false;
};

template <typename T, Index M, Index N>
struct is_matrix<Matrix<T, M, N>> {
  static bool const value = true;
};

template <typename T, Index M, Index N>
struct apply_matrix {
  typedef Matrix<typename T::type, M, N> type;
};

/// Tensors from 1st to 4th order and matrix
template <typename T>
struct order_1234 {
  static bool const value = false;
};

template <typename T, Index N>
struct order_1234<Vector<T, N>> {
  static bool const value = true;
};

template <typename T, Index N>
struct order_1234<Tensor<T, N>> {
  static bool const value = true;
};

template <typename T, Index N>
struct order_1234<Tensor3<T, N>> {
  static bool const value = true;
};

template <typename T, Index N>
struct order_1234<Tensor4<T, N>> {
 static bool const value = true;
  };                        

template<typename T, Index M, Index N>
struct order_1234<Matrix<T, M, N>>{
  static bool const value = true;
};

/// For Sacado traits

using std::string;

template<Index N>
struct dimension_string {
  static string eval() {return string("INVALID");}
};

template<>
struct dimension_string<DYNAMIC> {
  static string eval() {return string("DYNAMIC");}
};

template<>
struct dimension_string<1> {
  static string eval() {return string("1");}
};

template<>
struct dimension_string<2> {
  static string eval() {return string("2");}
};

template<>
struct dimension_string<3> {
  static string eval() {return string("3");}
};

template<>
struct dimension_string<4> {
  static string eval() {return string("4");}
};

} // namespace minitensor

namespace Sacado {

using minitensor::DYNAMIC;
using minitensor::Index;
using minitensor::Vector;
using minitensor::Tensor;
using minitensor::Tensor3;
using minitensor::Tensor4;
using minitensor::Matrix;
using minitensor::dimension_string;
using std::complex;
using std::string;

/// Specialization of Promote for Index
template<>
struct Promote<double, Index> {
  typedef double type;
};

template<>
struct Promote<Index, double> {
  typedef double type;
};

template<>
struct Promote<float, Index> {
  typedef float type;
};

template<>
struct Promote<Index, float> {
  typedef float type;
};

template<>
struct Promote<complex<double>, Index> {
  typedef complex<double> type;
};

template<>
struct Promote<Index, complex<double>> {
  typedef complex<double> type;
};

template<>
struct Promote<complex<float>, Index> {
  typedef complex<float> type;
};

template<>
struct Promote<Index, complex<float>> {
  typedef complex<float> type;
};

/// Sacado traits specializations for Vector
template <typename T, Index N>
struct ScalarType<Vector<T, N>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType<Vector<T, N>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType<Vector<T, N>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType<Vector<T, N>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value<Vector<T, N>> {
  typedef typename ValueType<Vector<T, N>>::type value_type;
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
struct ScalarValue<Vector<T, N>> {
  typedef typename ScalarType<Vector<T, N>>::type scalar_type;
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
struct StringName<Vector<T, N>> {
  static string
  eval()
  {
    return string("Vector<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual<Vector<T, N>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized<Vector<T, N>> {
  static bool const value = true;
};

template <typename T>
struct IsStaticallySized<Vector<T, DYNAMIC>>
{
  static bool const value = false;
};

/// Sacado traits specializations for Tensor
template <typename T, Index N>
struct ScalarType<Tensor<T, N>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType<Tensor<T, N>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType<Tensor<T, N>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType<Tensor<T, N>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value<Tensor<T, N>> {
  typedef typename ValueType<Tensor<T, N>>::type value_type;
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
struct ScalarValue<Tensor<T, N>> {
  typedef typename ScalarType<Tensor<T, N>>::type scalar_type;
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
struct StringName<Tensor<T, N>> {
  static string
  eval()
  {
    return string("Tensor<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual<Tensor<T, N>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized<Tensor<T, N>> {
  static bool const value = true;
};

template <typename T>
struct IsStaticallySized<Tensor<T, DYNAMIC>>
{
  static bool const value = false;
};

/// Sacado traits specializations for Tensor3
template <typename T, Index N>
struct ScalarType<Tensor3<T, N>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType<Tensor3<T, N>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType<Tensor3<T, N>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType<Tensor3<T, N>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value<Tensor3<T, N>> {
  typedef typename ValueType<Tensor3<T, N>>::type value_type;
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
struct ScalarValue<Tensor3<T, N>> {
  typedef typename ScalarType<Tensor3<T, N>>::type scalar_type;
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
struct StringName<Tensor3<T, N>> {
  static string
  eval()
  {
    return string("Tensor3<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual<Tensor3<T, N>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized<Tensor3<T, N>>
{
  static bool const value = true;
};

template <typename T>
struct IsStaticallySized<Tensor3<T, DYNAMIC>>
{
  static bool const value = false;
};

/// Sacado traits specializations for Tensor4
template <typename T, Index N>
struct ScalarType<Tensor4<T, N>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N>
struct ValueType<Tensor4<T, N>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N>
struct IsADType<Tensor4<T, N>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N>
struct IsScalarType<Tensor4<T, N>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N>
struct Value<Tensor4<T, N>> {
  typedef typename ValueType<Tensor4<T, N>>::type value_type;
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
struct ScalarValue<Tensor4<T, N>> {
  typedef typename ScalarType<Tensor4<T, N>>::type scalar_type;
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
struct StringName<Tensor4<T, N>> {
  static string
  eval()
  {
    return string("Tensor4<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N>
struct IsEqual<Tensor4<T, N>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N>
struct IsStaticallySized<Tensor4<T, N>>
{
  static bool const value = true;
};

template <typename T>
struct IsStaticallySized<Tensor4<T, DYNAMIC>>
{
 static bool const value= false;
};

/// Sacado traits specializations for Matrix
template <typename T, Index M, Index N>
struct ScalarType<Matrix<T, M, N>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index M, Index N>
struct ValueType<Matrix<T, M, N>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index M, Index N>
struct IsADType<Matrix<T, M, N>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index M, Index N>
struct IsScalarType<Matrix<T, M, N>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index M, Index N>
struct Value<Matrix<T, M, N>> {
  typedef typename ValueType<Matrix<T, M, N>>::type value_type;
  static const Matrix<value_type, M, N>
  eval(Matrix<T, M, N> const & x)
  {
    Matrix<value_type, M, N> v(x.get_num_rows(), x.get_num_cols());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index M, Index N>
struct ScalarValue<Matrix<T, M, N>> {
  typedef typename ScalarType<Matrix<T, M, N>>::type scalar_type;
  static const Matrix<scalar_type, M, N>
  eval(Matrix<T, M, N> const & x)
  {
    Matrix<scalar_type, M, N> v(x.get_num_rows(), x.get_num_cols());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index M, Index N>
struct StringName<Matrix<T, M, N>> {
  static string
  eval()
  {
    return string("Matrix<") + StringName<T>::eval() + string(", ") +
        dimension_string<M>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index M, Index N>
struct IsEqual<Matrix<T, M, N>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index M, Index N>
struct IsStaticallySized<Matrix<T, M, N>> {
  static bool const value = true;
};

template <typename T, Index M>
struct IsStaticallySized<Matrix<T, M, DYNAMIC>>
{
  static bool const value = false;
};

template <typename T, Index N>
struct IsStaticallySized<Matrix<T, DYNAMIC, N>>
{
  static bool const value = false;
};

template <typename T>
struct IsStaticallySized<Matrix<T, DYNAMIC, DYNAMIC>>
{
  static bool const value = false;
};

} // namespace Sacado

#endif // MiniTensor_Definitions_h
