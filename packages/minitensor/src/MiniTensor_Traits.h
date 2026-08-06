// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Traits_h)
#define MiniTensor_Traits_h

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

/// \addtogroup minitensor_traits
/// @{

/// Indexing type
using Index = uint32_t;

///
/// Number of bits in Index.
///
constexpr Index
INDEX_SIZE{32};

/// High count type
using LongIndex = uint64_t;

///
/// Number of bits in LongIndex.
///
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
  /// Whether T is a Vector.
  static bool const value = false;
};

/// is_vector specialization: T is a Vector.
template <typename T, Index N>
struct is_vector<Vector<T, N>> {
  /// Whether T is a Vector.
  static bool const value = true;
};

///
/// Metafunction that forms a Vector of dimension N whose element
/// type is the nested type of the metafunction T.
///
template <typename T, Index N>
struct apply_vector {
  /// Vector with element type T::type.
  typedef Vector<typename T::type, N> type;
};

/// 2nd-order tensor
template <typename T>
struct is_tensor {
  /// Whether T is a Tensor.
  static bool const value = false;
};

/// is_tensor specialization: T is a Tensor.
template <typename T, Index N>
struct is_tensor<Tensor<T, N>> {
  /// Whether T is a Tensor.
  static bool const value = true;
};

///
/// Metafunction that forms a Tensor of dimension N whose element
/// type is the nested type of the metafunction T.
///
template <typename T, Index N>
struct apply_tensor {
  /// Tensor with element type T::type.
  typedef Tensor<typename T::type, N> type;
};

/// 3rd-order tensor
template <typename T>
struct is_tensor3 {
  /// Whether T is a Tensor3.
  static bool const value = false;
};

/// is_tensor3 specialization: T is a Tensor3.
template <typename T, Index N>
struct is_tensor3<Tensor3<T, N>> {
  /// Whether T is a Tensor3.
  static bool const value = true;
};

///
/// Metafunction that forms a Tensor3 of dimension N whose element
/// type is the nested type of the metafunction T.
///
template <typename T, Index N>
struct apply_tensor3 {
  /// Tensor3 with element type T::type.
  typedef Tensor3<typename T::type, N> type;
};

/// 4th-order tensor
template <typename T>
struct is_tensor4 {
  /// Whether T is a Tensor4.
  static bool const value = false;
};

/// is_tensor4 specialization: T is a Tensor4.
template <typename T, Index N>
struct is_tensor4<Tensor4<T, N>> {
  /// Whether T is a Tensor4.
  static bool const value = true;
};

///
/// Metafunction that forms a Tensor4 of dimension N whose element
/// type is the nested type of the metafunction T.
///
template <typename T, Index N>
struct apply_tensor4 {
  /// Tensor4 with element type T::type.
  typedef Tensor4<typename T::type, N> type;
};

/// Matrix
template <typename T>
struct is_matrix {
  /// Whether T is a Matrix.
  static bool const value = false;
};

/// is_matrix specialization: T is a Matrix.
template <typename T, Index M, Index N>
struct is_matrix<Matrix<T, M, N>> {
  /// Whether T is a Matrix.
  static bool const value = true;
};

///
/// Metafunction that forms an M by N Matrix whose element type
/// is the nested type of the metafunction T.
///
template <typename T, Index M, Index N>
struct apply_matrix {
  /// Matrix with element type T::type.
  typedef Matrix<typename T::type, M, N> type;
};

/// Tensors from 1st to 4th order and matrix
template <typename T>
struct order_1234 {
  /// Whether T is a tensor of order 1 to 4 or a matrix.
  static bool const value = false;
};

/// order_1234 specialization for Vector.
template <typename T, Index N>
struct order_1234<Vector<T, N>> {
  /// Whether T is a tensor of order 1 to 4 or a matrix.
  static bool const value = true;
};

/// order_1234 specialization for Tensor.
template <typename T, Index N>
struct order_1234<Tensor<T, N>> {
  /// Whether T is a tensor of order 1 to 4 or a matrix.
  static bool const value = true;
};

/// order_1234 specialization for Tensor3.
template <typename T, Index N>
struct order_1234<Tensor3<T, N>> {
  /// Whether T is a tensor of order 1 to 4 or a matrix.
  static bool const value = true;
};

/// order_1234 specialization for Tensor4.
template <typename T, Index N>
struct order_1234<Tensor4<T, N>> {
 /// Whether T is a tensor of order 1 to 4 or a matrix.
 static bool const value = true;
  };                        

/// order_1234 specialization for Matrix.
template<typename T, Index M, Index N>
struct order_1234<Matrix<T, M, N>>{
  /// Whether T is a tensor of order 1 to 4 or a matrix.
  static bool const value = true;
};

/// For Sacado traits

using std::string;

///
/// Compile-time conversion of a dimension N to a string, used by
/// the Sacado StringName specializations below.
///
template<Index N>
struct dimension_string {
  /// Return "INVALID" for dimensions without a specialization.
  static string eval() {return string("INVALID");}
};

/// dimension_string specialization for DYNAMIC.
template<>
struct dimension_string<DYNAMIC> {
  /// Return "DYNAMIC".
  static string eval() {return string("DYNAMIC");}
};

/// dimension_string specialization for dimension 1.
template<>
struct dimension_string<1> {
  /// Return "1".
  static string eval() {return string("1");}
};

/// dimension_string specialization for dimension 2.
template<>
struct dimension_string<2> {
  /// Return "2".
  static string eval() {return string("2");}
};

/// dimension_string specialization for dimension 3.
template<>
struct dimension_string<3> {
  /// Return "3".
  static string eval() {return string("3");}
};

/// dimension_string specialization for dimension 4.
template<>
struct dimension_string<4> {
  /// Return "4".
  static string eval() {return string("4");}
};

/// @}
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
  /// Promoted type.
  typedef double type;
};

/// Promote specialization for Index and double.
template<>
struct Promote<Index, double> {
  /// Promoted type.
  typedef double type;
};

/// Promote specialization for float and Index.
template<>
struct Promote<float, Index> {
  /// Promoted type.
  typedef float type;
};

/// Promote specialization for Index and float.
template<>
struct Promote<Index, float> {
  /// Promoted type.
  typedef float type;
};

/// Promote specialization for complex<double> and Index.
template<>
struct Promote<complex<double>, Index> {
  /// Promoted type.
  typedef complex<double> type;
};

/// Promote specialization for Index and complex<double>.
template<>
struct Promote<Index, complex<double>> {
  /// Promoted type.
  typedef complex<double> type;
};

/// Promote specialization for complex<float> and Index.
template<>
struct Promote<complex<float>, Index> {
  /// Promoted type.
  typedef complex<float> type;
};

/// Promote specialization for Index and complex<float>.
template<>
struct Promote<Index, complex<float>> {
  /// Promoted type.
  typedef complex<float> type;
};

/// Sacado traits specializations for Vector
template <typename T, Index N>
struct ScalarType<Vector<T, N>> {
  /// Underlying scalar type of the Vector components.
  typedef typename ScalarType<T>::type type;
};

/// Sacado ValueType specialization for Vector.
template <typename T, Index N>
struct ValueType<Vector<T, N>> {
  /// Value type of the Vector components.
  typedef typename ValueType<T>::type type;
};

/// Sacado IsADType specialization for Vector.
template <typename T, Index N>
struct IsADType<Vector<T, N>> {
  /// Whether the component type is an AD type.
  static bool const value = IsADType<T>::value;
};

/// Sacado IsScalarType specialization for Vector.
template <typename T, Index N>
struct IsScalarType<Vector<T, N>> {
  /// Whether the component type is a scalar type.
  static bool const value = IsScalarType<T>::value;
};

/// Sacado Value specialization for Vector.
template <typename T, Index N>
struct Value<Vector<T, N>> {
  /// Value type of the Vector components.
  typedef typename ValueType<Vector<T, N>>::type value_type;
  /// Extract a Vector of the values of the components of x.
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

/// Sacado ScalarValue specialization for Vector.
template <typename T, Index N>
struct ScalarValue<Vector<T, N>> {
  /// Scalar type of the Vector components.
  typedef typename ScalarType<Vector<T, N>>::type scalar_type;
  /// Extract a Vector of the scalar values of the components
  /// of x.
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

/// Sacado StringName specialization for Vector.
template <typename T, Index N>
struct StringName<Vector<T, N>> {
  /// Return the name of the Vector type as a string.
  static string
  eval()
  {
    return string("Vector<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

/// Sacado IsEqual specialization for Vector.
template <typename T, Index N>
struct IsEqual<Vector<T, N>> {
  /// Compare two components for equality.
  static bool eval(T const & x, T const & y) { return x == y; }
};

/// Sacado IsStaticallySized specialization for Vector.
template <typename T, Index N>
struct IsStaticallySized<Vector<T, N>> {
  /// True: the dimension is fixed at compile time.
  static bool const value = true;
};

/// Sacado IsStaticallySized specialization for dynamic Vector.
template <typename T>
struct IsStaticallySized<Vector<T, DYNAMIC>>
{
  /// False: the dimension is determined at run time.
  static bool const value = false;
};

/// Sacado traits specializations for Tensor
template <typename T, Index N>
struct ScalarType<Tensor<T, N>> {
  /// Underlying scalar type of the Tensor components.
  typedef typename ScalarType<T>::type type;
};

/// Sacado ValueType specialization for Tensor.
template <typename T, Index N>
struct ValueType<Tensor<T, N>> {
  /// Value type of the Tensor components.
  typedef typename ValueType<T>::type type;
};

/// Sacado IsADType specialization for Tensor.
template <typename T, Index N>
struct IsADType<Tensor<T, N>> {
  /// Whether the component type is an AD type.
  static bool const value = IsADType<T>::value;
};

/// Sacado IsScalarType specialization for Tensor.
template <typename T, Index N>
struct IsScalarType<Tensor<T, N>> {
  /// Whether the component type is a scalar type.
  static bool const value = IsScalarType<T>::value;
};

/// Sacado Value specialization for Tensor.
template <typename T, Index N>
struct Value<Tensor<T, N>> {
  /// Value type of the Tensor components.
  typedef typename ValueType<Tensor<T, N>>::type value_type;
  /// Extract a Tensor of the values of the components of x.
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

/// Sacado ScalarValue specialization for Tensor.
template <typename T, Index N>
struct ScalarValue<Tensor<T, N>> {
  /// Scalar type of the Tensor components.
  typedef typename ScalarType<Tensor<T, N>>::type scalar_type;
  /// Extract a Tensor of the scalar values of the components
  /// of x.
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

/// Sacado StringName specialization for Tensor.
template <typename T, Index N>
struct StringName<Tensor<T, N>> {
  /// Return the name of the Tensor type as a string.
  static string
  eval()
  {
    return string("Tensor<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

/// Sacado IsEqual specialization for Tensor.
template <typename T, Index N>
struct IsEqual<Tensor<T, N>> {
  /// Compare two components for equality.
  static bool eval(T const & x, T const & y) { return x == y; }
};

/// Sacado IsStaticallySized specialization for Tensor.
template <typename T, Index N>
struct IsStaticallySized<Tensor<T, N>> {
  /// True: the dimension is fixed at compile time.
  static bool const value = true;
};

/// Sacado IsStaticallySized specialization for dynamic Tensor.
template <typename T>
struct IsStaticallySized<Tensor<T, DYNAMIC>>
{
  /// False: the dimension is determined at run time.
  static bool const value = false;
};

/// Sacado traits specializations for Tensor3
template <typename T, Index N>
struct ScalarType<Tensor3<T, N>> {
  /// Underlying scalar type of the Tensor3 components.
  typedef typename ScalarType<T>::type type;
};

/// Sacado ValueType specialization for Tensor3.
template <typename T, Index N>
struct ValueType<Tensor3<T, N>> {
  /// Value type of the Tensor3 components.
  typedef typename ValueType<T>::type type;
};

/// Sacado IsADType specialization for Tensor3.
template <typename T, Index N>
struct IsADType<Tensor3<T, N>> {
  /// Whether the component type is an AD type.
  static bool const value = IsADType<T>::value;
};

/// Sacado IsScalarType specialization for Tensor3.
template <typename T, Index N>
struct IsScalarType<Tensor3<T, N>> {
  /// Whether the component type is a scalar type.
  static bool const value = IsScalarType<T>::value;
};

/// Sacado Value specialization for Tensor3.
template <typename T, Index N>
struct Value<Tensor3<T, N>> {
  /// Value type of the Tensor3 components.
  typedef typename ValueType<Tensor3<T, N>>::type value_type;
  /// Extract a Tensor3 of the values of the components of x.
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

/// Sacado ScalarValue specialization for Tensor3.
template <typename T, Index N>
struct ScalarValue<Tensor3<T, N>> {
  /// Scalar type of the Tensor3 components.
  typedef typename ScalarType<Tensor3<T, N>>::type scalar_type;
  /// Extract a Tensor3 of the scalar values of the components
  /// of x.
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

/// Sacado StringName specialization for Tensor3.
template <typename T, Index N>
struct StringName<Tensor3<T, N>> {
  /// Return the name of the Tensor3 type as a string.
  static string
  eval()
  {
    return string("Tensor3<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

/// Sacado IsEqual specialization for Tensor3.
template <typename T, Index N>
struct IsEqual<Tensor3<T, N>> {
  /// Compare two components for equality.
  static bool eval(T const & x, T const & y) { return x == y; }
};

/// Sacado IsStaticallySized specialization for Tensor3.
template <typename T, Index N>
struct IsStaticallySized<Tensor3<T, N>>
{
  /// True: the dimension is fixed at compile time.
  static bool const value = true;
};

/// Sacado IsStaticallySized specialization for dynamic Tensor3.
template <typename T>
struct IsStaticallySized<Tensor3<T, DYNAMIC>>
{
  /// False: the dimension is determined at run time.
  static bool const value = false;
};

/// Sacado traits specializations for Tensor4
template <typename T, Index N>
struct ScalarType<Tensor4<T, N>> {
  /// Underlying scalar type of the Tensor4 components.
  typedef typename ScalarType<T>::type type;
};

/// Sacado ValueType specialization for Tensor4.
template <typename T, Index N>
struct ValueType<Tensor4<T, N>> {
  /// Value type of the Tensor4 components.
  typedef typename ValueType<T>::type type;
};

/// Sacado IsADType specialization for Tensor4.
template <typename T, Index N>
struct IsADType<Tensor4<T, N>> {
  /// Whether the component type is an AD type.
  static bool const value = IsADType<T>::value;
};

/// Sacado IsScalarType specialization for Tensor4.
template <typename T, Index N>
struct IsScalarType<Tensor4<T, N>> {
  /// Whether the component type is a scalar type.
  static bool const value = IsScalarType<T>::value;
};

/// Sacado Value specialization for Tensor4.
template <typename T, Index N>
struct Value<Tensor4<T, N>> {
  /// Value type of the Tensor4 components.
  typedef typename ValueType<Tensor4<T, N>>::type value_type;
  /// Extract a Tensor4 of the values of the components of x.
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

/// Sacado ScalarValue specialization for Tensor4.
template <typename T, Index N>
struct ScalarValue<Tensor4<T, N>> {
  /// Scalar type of the Tensor4 components.
  typedef typename ScalarType<Tensor4<T, N>>::type scalar_type;
  /// Extract a Tensor4 of the scalar values of the components
  /// of x.
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

/// Sacado StringName specialization for Tensor4.
template <typename T, Index N>
struct StringName<Tensor4<T, N>> {
  /// Return the name of the Tensor4 type as a string.
  static string
  eval()
  {
    return string("Tensor4<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

/// Sacado IsEqual specialization for Tensor4.
template <typename T, Index N>
struct IsEqual<Tensor4<T, N>> {
  /// Compare two components for equality.
  static bool eval(T const & x, T const & y) { return x == y; }
};

/// Sacado IsStaticallySized specialization for Tensor4.
template <typename T, Index N>
struct IsStaticallySized<Tensor4<T, N>>
{
  /// True: the dimension is fixed at compile time.
  static bool const value = true;
};

/// Sacado IsStaticallySized specialization for dynamic Tensor4.
template <typename T>
struct IsStaticallySized<Tensor4<T, DYNAMIC>>
{
 /// False: the dimension is determined at run time.
 static bool const value= false;
};

/// Sacado traits specializations for Matrix
template <typename T, Index M, Index N>
struct ScalarType<Matrix<T, M, N>> {
  /// Underlying scalar type of the Matrix components.
  typedef typename ScalarType<T>::type type;
};

/// Sacado ValueType specialization for Matrix.
template <typename T, Index M, Index N>
struct ValueType<Matrix<T, M, N>> {
  /// Value type of the Matrix components.
  typedef typename ValueType<T>::type type;
};

/// Sacado IsADType specialization for Matrix.
template <typename T, Index M, Index N>
struct IsADType<Matrix<T, M, N>> {
  /// Whether the component type is an AD type.
  static bool const value = IsADType<T>::value;
};

/// Sacado IsScalarType specialization for Matrix.
template <typename T, Index M, Index N>
struct IsScalarType<Matrix<T, M, N>> {
  /// Whether the component type is a scalar type.
  static bool const value = IsScalarType<T>::value;
};

/// Sacado Value specialization for Matrix.
template <typename T, Index M, Index N>
struct Value<Matrix<T, M, N>> {
  /// Value type of the Matrix components.
  typedef typename ValueType<Matrix<T, M, N>>::type value_type;
  /// Extract a Matrix of the values of the components of x.
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

/// Sacado ScalarValue specialization for Matrix.
template <typename T, Index M, Index N>
struct ScalarValue<Matrix<T, M, N>> {
  /// Scalar type of the Matrix components.
  typedef typename ScalarType<Matrix<T, M, N>>::type scalar_type;
  /// Extract a Matrix of the scalar values of the components
  /// of x.
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

/// Sacado StringName specialization for Matrix.
template <typename T, Index M, Index N>
struct StringName<Matrix<T, M, N>> {
  /// Return the name of the Matrix type as a string.
  static string
  eval()
  {
    return string("Matrix<") + StringName<T>::eval() + string(", ") +
        dimension_string<M>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

/// Sacado IsEqual specialization for Matrix.
template <typename T, Index M, Index N>
struct IsEqual<Matrix<T, M, N>> {
  /// Compare two components for equality.
  static bool eval(T const & x, T const & y) { return x == y; }
};

/// Sacado IsStaticallySized specialization for Matrix.
template <typename T, Index M, Index N>
struct IsStaticallySized<Matrix<T, M, N>> {
  /// True: the dimensions are fixed at compile time.
  static bool const value = true;
};

/// Sacado IsStaticallySized specialization for Matrix with a
/// dynamic number of columns.
template <typename T, Index M>
struct IsStaticallySized<Matrix<T, M, DYNAMIC>>
{
  /// False: the number of columns is determined at run time.
  static bool const value = false;
};

/// Sacado IsStaticallySized specialization for Matrix with a
/// dynamic number of rows.
template <typename T, Index N>
struct IsStaticallySized<Matrix<T, DYNAMIC, N>>
{
  /// False: the number of rows is determined at run time.
  static bool const value = false;
};

/// Sacado IsStaticallySized specialization for Matrix with
/// dynamic numbers of rows and columns.
template <typename T>
struct IsStaticallySized<Matrix<T, DYNAMIC, DYNAMIC>>
{
  /// False: the dimensions are determined at run time.
  static bool const value = false;
};

} // namespace Sacado

#endif // MiniTensor_Traits_h
