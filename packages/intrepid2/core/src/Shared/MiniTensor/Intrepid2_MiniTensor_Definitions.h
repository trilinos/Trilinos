// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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

#if !defined(Intrepid2_MiniTensor_Definitions_h)
#define Intrepid2_MiniTensor_Definitions_h

#include <complex>
#include <type_traits>

#include "Intrepid2_ConfigDefs.hpp"
#include "Sacado.hpp"


#define         NPP_MAX_32U   (32 )

namespace Intrepid2 {

/// Indexing type
using Index = unsigned int;

constexpr Index
INDEX_SIZE{32};

/// High count type
using LongIndex = unsigned long int;

constexpr Index
LONG_INDEX_SIZE{64};

/// Floating point type
using Real = double;

/// Complex type
using Complex = std::complex<Real>;

/// The classes
template <typename T, Index N,  typename ES> class Vector;
template <typename T, Index N,  typename ES> class Tensor;
template <typename T, Index N,  typename ES> class Tensor3;
template <typename T, Index N,  typename ES> class Tensor4;
template <typename T, Index M, Index N,  typename ES> class Matrix;

/// Indicator for dynamic storage
constexpr Index
DYNAMIC{0};

// Default execution space
class NOKOKKOS{};

/// For use with type promotion
using Sacado::Promote;
using Sacado::mpl::lazy_disable_if;
using Sacado::mpl::disable_if_c;

/// Vector
template <typename T>
struct is_vector {
  static const bool value = false;
};

template <typename T, Index N,  typename ES>
struct is_vector< Vector<T, N, ES>> {
  static const bool value = true;
};

template <typename T, Index N,  typename ES>
struct apply_vector {
  typedef Vector<typename T::type, N, ES> type;
};

/// 2nd-order tensor
template <typename T>
struct is_tensor {
  static const bool value = false;
};

template <typename T, Index N,  typename ES>
struct is_tensor< Tensor<T, N, ES>> {
  static const bool value = true;
};

template <typename T, Index N,  typename ES>
struct apply_tensor {
  typedef Tensor<typename T::type, N, ES> type;
};

/// 3rd-order tensor
template <typename T>
struct is_tensor3 {
  static const bool value = false;
};

template <typename T, Index N,  typename ES>
struct is_tensor3< Tensor3<T, N, ES>> {
  static const bool value = true;
};

template <typename T, Index N,  typename ES>
struct apply_tensor3 {
  typedef Tensor3<typename T::type, N, ES> type;
};

/// 4th-order tensor
template <typename T>
struct is_tensor4 {
  static const bool value = false;
};

template <typename T, Index N,  typename ES>
struct is_tensor4<Tensor4<T, N, ES>> {
  static const bool value = true;
};

template <typename T, Index N,  typename ES>
struct apply_tensor4 {
  typedef Tensor4<typename T::type, N, ES> type;
};

/// Matrix
template <typename T>
struct is_matrix {
  static const bool value = false;
};

template <typename T, Index M, Index N,  typename ES>
struct is_matrix<Matrix<T, M, N, ES>> {
  static const bool value = true;
};

template <typename T, Index M, Index N,  typename ES>
struct apply_matrix {
  typedef Matrix<typename T::type, M, N, ES> type;
};

/// Tensors from 1st to 4th order and matrix
template <typename T>
struct order_1234 {
  static const bool value = false;
};

template <typename T, Index N,  typename ES>
struct order_1234< Vector<T, N, ES>> {
  static const bool value = true;
};

template <typename T, Index N,  typename ES>
struct order_1234< Tensor<T, N, ES>> {
  static const bool value = true;
};

template <typename T, Index N,  typename ES>
struct order_1234< Tensor3<T, N, ES>> {
  static const bool value = true;
};

template <typename T, Index N,  typename ES>
struct order_1234< Tensor4<T, N, ES>> {
 static const bool value = true;          
  };                        

template<typename T, Index M, Index N,  typename ES>
struct order_1234<Matrix<T, M, N, ES>>{
  static const bool value = true;
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

} // namespace Intrepid

namespace Sacado {

using Intrepid2::DYNAMIC;
using Intrepid2::NOKOKKOS;
using Intrepid2::Index;
using Intrepid2::Vector;
using Intrepid2::Tensor;
using Intrepid2::Tensor3;
using Intrepid2::Tensor4;
using Intrepid2::Matrix;
using Intrepid2::dimension_string;
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
template <typename T, Index N,  typename ES>
struct ScalarType< Vector<T, N, ES>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct ValueType< Vector<T, N, ES>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct IsADType< Vector<T, N, ES>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N,  typename ES>
struct IsScalarType< Vector<T, N, ES>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N,  typename ES>
struct Value< Vector<T, N, ES>> {
  typedef typename ValueType< Vector<T, N, ES>>::type value_type;
  static const Vector<value_type, N, ES>
  eval(Vector<T, N, ES> const & x)
  {
    Vector<value_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N,  typename ES>
struct ScalarValue< Vector<T, N, ES>> {
  typedef typename ScalarType< Vector<T, N, ES>>::type scalar_type;
  static const Vector<scalar_type, N, ES>
  eval(Vector<T, N, ES> const & x)
  {
    Vector<scalar_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }
    return v;
  }
};

template <typename T, Index N,  typename ES>
struct StringName<Vector<T, N, ES>> {
  static string
  eval()
  {
    return string("Vector<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N,  typename ES>
struct IsEqual< Vector<T, N, ES>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N,  typename ES>
struct IsStaticallySized< Vector<T, N, ES>> {
  static const bool value = true;
};

template <typename T,  typename ES>
struct IsStaticallySized< Vector<T, DYNAMIC,ES>>
{
  static const bool value = false;
};

/// Sacado traits specializations for Tensor
template <typename T, Index N,  typename ES>
struct ScalarType< Tensor<T, N, ES>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct ValueType< Tensor<T, N, ES>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct IsADType< Tensor<T, N, ES>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N,  typename ES>
struct IsScalarType< Tensor<T, N, ES>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N,  typename ES>
struct Value< Tensor<T, N, ES>> {
  typedef typename ValueType< Tensor<T, N, ES>>::type value_type;
  static const Tensor<value_type, N, ES>
  eval(Tensor<T, N, ES> const & x)
  {
    Tensor<value_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N,  typename ES>
struct ScalarValue< Tensor<T, N, ES>> {
  typedef typename ScalarType< Tensor<T, N, ES>>::type scalar_type;
  static const Tensor<scalar_type, N, ES>
  eval(Tensor<T, N, ES> const & x)
  {
    Tensor<scalar_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N,  typename ES>
struct StringName< Tensor<T, N, ES>> {
  static string
  eval()
  {
    return string("Tensor<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N,  typename ES>
struct IsEqual< Tensor<T, N, ES>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N,  typename ES>
struct IsStaticallySized< Tensor<T, N, ES>> {
  static const bool value = true;
};

template <typename T,  typename ES>
struct IsStaticallySized< Tensor<T, DYNAMIC, ES>>
{
  static const bool value = false;
};

/// Sacado traits specializations for Tensor3
template <typename T, Index N,  typename ES>
struct ScalarType< Tensor3<T, N, ES>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct ValueType< Tensor3<T, N, ES>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct IsADType< Tensor3<T, N, ES>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N,  typename ES>
struct IsScalarType< Tensor3<T, N, ES>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N,  typename ES>
struct Value< Tensor3<T, N, ES>> {
  typedef typename ValueType< Tensor3<T, N, ES>>::type value_type;
  static const Tensor3<value_type, N, ES>
  eval(Tensor3<T, N, ES> const & x)
  {
    Tensor3<value_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N,  typename ES>
struct ScalarValue< Tensor3<T, N, ES>> {
  typedef typename ScalarType< Tensor3<T, N, ES>>::type scalar_type;
  static const Tensor3<scalar_type, N, ES>
  eval(Tensor3<T, N, ES> const & x)
  {
    Tensor3<scalar_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N,  typename ES>
struct StringName< Tensor3<T, N, ES>> {
  static string
  eval()
  {
    return string("Tensor3<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N,  typename ES>
struct IsEqual< Tensor3<T, N, ES>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N,  typename ES>
struct IsStaticallySized< Tensor3<T, N, ES>>
{
  static const bool value = true;
};

template <typename T,  typename ES>
struct IsStaticallySized< Tensor3<T, DYNAMIC,ES>>
{
  static const bool value = false;
};

/// Sacado traits specializations for Tensor4
template <typename T, Index N,  typename ES>
struct ScalarType< Tensor4<T, N, ES>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct ValueType< Tensor4<T, N, ES>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index N,  typename ES>
struct IsADType< Tensor4<T, N, ES>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index N,  typename ES>
struct IsScalarType< Tensor4<T, N, ES>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index N,  typename ES>
struct Value< Tensor4<T, N, ES>> {
  typedef typename ValueType< Tensor4<T, N, ES>>::type value_type;
  static const Tensor4<value_type, N, ES>
  eval(Tensor4<T, N, ES> const & x)
  {
    Tensor4<value_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N,  typename ES>
struct ScalarValue< Tensor4<T, N, ES>> {
  typedef typename ScalarType< Tensor4<T, N, ES>>::type scalar_type;
  static const Tensor4<scalar_type, N, ES>
  eval(Tensor4<T, N, ES> const & x)
  {
    Tensor4<scalar_type, N, ES> v(x.get_dimension());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index N,  typename ES>
struct StringName< Tensor4<T, N, ES>> {
  static string
  eval()
  {
    return string("Tensor4<") + StringName<T>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index N,  typename ES>
struct IsEqual< Tensor4<T, N, ES>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index N,  typename ES>
struct IsStaticallySized< Tensor4<T, N, ES>>
{
  static const bool value = true;
};

template <typename T,  typename ES>
struct IsStaticallySized< Tensor4<T, DYNAMIC, ES>>
{
 static const bool value= false;
};

/// Sacado traits specializations for Matrix
template <typename T, Index M, Index N,  typename ES>
struct ScalarType<Matrix<T, M, N, ES>> {
  typedef typename ScalarType<T>::type type;
};

template <typename T, Index M, Index N,  typename ES>
struct ValueType<Matrix<T, M, N, ES>> {
  typedef typename ValueType<T>::type type;
};

template <typename T, Index M, Index N,  typename ES>
struct IsADType<Matrix<T, M, N, ES>> {
  static bool const value = IsADType<T>::value;
};

template <typename T, Index M, Index N,  typename ES>
struct IsScalarType<Matrix<T, M, N, ES>> {
  static bool const value = IsScalarType<T>::value;
};

template <typename T, Index M, Index N,  typename ES>
struct Value<Matrix<T, M, N, ES>> {
  typedef typename ValueType<Matrix<T, M, N, ES>>::type value_type;
  static const Matrix<value_type, M, N, ES>
  eval(Matrix<T, M, N, ES> const & x)
  {
    Matrix<value_type, M, N, ES> v(x.get_num_rows(), x.get_num_cols());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = Value<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index M, Index N,  typename ES>
struct ScalarValue<Matrix<T, M, N, ES>> {
  typedef typename ScalarType<Matrix<T, M, N, ES>>::type scalar_type;
  static const Matrix<scalar_type, M, N, ES>
  eval(Matrix<T, M, N, ES> const & x)
  {
    Matrix<scalar_type, M, N, ES> v(x.get_num_rows(), x.get_num_cols());

    for (Index i = 0; i < x.get_number_components(); ++i) {
      v[i] = ScalarValue<T>::eval(x[i]);
    }

    return v;
  }
};

template <typename T, Index M, Index N,  typename ES>
struct StringName<Matrix<T, M, N, ES>> {
  static string
  eval()
  {
    return string("Matrix<") + StringName<T>::eval() + string(", ") +
        dimension_string<M>::eval() + string(", ") +
        dimension_string<N>::eval() + string(">");
  }
};

template <typename T, Index M, Index N,  typename ES>
struct IsEqual<Matrix<T, M, N, ES>> {
  static bool eval(T const & x, T const & y) { return x == y; }
};

template <typename T, Index M, Index N,  typename ES>
struct IsStaticallySized<Matrix<T, M, N, ES>> {
  static const bool value = true;
};

template <typename T, Index M,  typename ES>
struct IsStaticallySized<Matrix<T, M, DYNAMIC, ES>>
{
  static const bool value = false;
};

template <typename T, Index N,  typename ES>
struct IsStaticallySized<Matrix<T, DYNAMIC, N, ES>>
{
  static const bool value = false;
};

template <typename T,  typename ES>
struct IsStaticallySized<Matrix<T, DYNAMIC, DYNAMIC, ES>>
{
  static const bool value = false;
};

} // namespace Sacado

#endif // Intrepid2_MiniTensor_Definitions_h
