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

namespace Intrepid {

  ///
  /// Indexing type
  ///
  typedef unsigned int Index;

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
  template <typename T> class Vector;
  template <typename T> class Tensor;
  template <typename T> class Tensor3;
  template <typename T> class Tensor4;

  ///
  /// Type promotion
  ///
  template <typename S, typename T>
  struct Promote
  {
    typedef typename Sacado::Promote<S, T> type;
  };

  template <typename S, typename T>
  struct Promote<Intrepid::Tensor<T>, S> {};

  template <typename S, typename T>
  struct Promote<S, Intrepid::Tensor<T> > {};

  /*
  using Sacado::Promote;

  template <typename S, typename T>
  struct Elevate
  {
    typedef typename Sacado::Promote<S, T> type;
  };

  template <typename S, typename T>
  struct Elevate<Intrepid::Tensor<T>, S> {};

  template <typename S, typename T>
  struct Elevate<S, Intrepid::Tensor<T> > {};
  */

} // namespace Intrepid

namespace Sacado {

  ///
  /// Sacado specializations for Vector
  ///
  template <typename T>
  struct ScalarType< Intrepid::Vector<T> > {
    typedef T type;
  };

  template <typename T>
  struct ValueType< Intrepid::Vector<T> > {
    typedef T type;
  };

  template <typename T>
  struct IsADType< Intrepid::Vector<T> > {
    static bool const value = IsADType<T>::value;
  };

  template <typename T>
  struct IsScalarType< Intrepid::Vector<T> > {
    static bool const value = false;
  };

  template <typename T>
  struct Value< Intrepid::Vector<T> > {
    typedef typename ValueType< T >::type value_type;
    static const Intrepid::Vector<value_type> &
    eval(Intrepid::Vector<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Vector<value_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        v(i) = x(i).val();
      }

      return v;
    }
  };

  template <typename T>
  struct ScalarValue< Intrepid::Vector<T> > {
    typedef typename ValueType< T >::type value_type;
    typedef typename ScalarType< T >::type scalar_type;
    static const Intrepid::Vector<scalar_type> &
    eval(Intrepid::Vector<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Vector<scalar_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        v(i) = ScalarValue<value_type>::eval(x(i).val());
      }

      return v;
    }
  };

  template <typename T>
  struct StringName< Intrepid::Vector<T> > {
    static std::string
    eval()
    {
      return std::string("Intrepid::Vector< ") +
          StringName<T>::eval() + std::string(" >");
    }
  };

  template <typename T>
  struct IsEqual< Intrepid::Vector<T> > {
    static bool eval(T const & x, T const & y) { return x == y; }
  };

  template <typename T>
  struct IsStaticallySized< Intrepid::Vector<T> > {
    static const bool value = false;
  };

  ///
  /// Sacado specializations for Tensor
  ///
  template <typename T>
  struct ScalarType< Intrepid::Tensor<T> > {
    typedef T type;
  };

  template <typename T>
  struct ValueType< Intrepid::Tensor<T> > {
    typedef T type;
  };

  template <typename T>
  struct IsADType< Intrepid::Tensor<T> > {
    static bool const value = IsADType<T>::value;
  };

  template <typename T>
  struct IsScalarType< Intrepid::Tensor<T> > {
    static bool const value = false;
  };

  template <typename T>
  struct Value< Intrepid::Tensor<T> > {
    typedef typename ValueType< T >::type value_type;
    static const Intrepid::Tensor<value_type> &
    eval(Intrepid::Tensor<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Tensor<value_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        for (Intrepid::Index j = 0; j < N; ++j) {
          v(i,j) = x(i,j).val();
        }
      }

      return v;
    }
  };

  template <typename T>
  struct ScalarValue< Intrepid::Tensor<T> > {
    typedef typename ValueType< T >::type value_type;
    typedef typename ScalarType< T >::type scalar_type;
    static const Intrepid::Tensor<scalar_type> &
    eval(Intrepid::Tensor<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Tensor<scalar_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        for (Intrepid::Index j = 0; j < N; ++j) {
          v(i,j) = ScalarValue<value_type>::eval(x(i,j).val());
        }
      }

      return v;
    }
  };

  template <typename T>
  struct StringName< Intrepid::Tensor<T> > {
    static std::string
    eval()
    {
      return std::string("Intrepid::Tensor< ") +
          StringName<T>::eval() + std::string(" >");
    }
  };

  template <typename T>
  struct IsEqual< Intrepid::Tensor<T> > {
    static bool eval(T const & x, T const & y) { return x == y; }
  };

  template <typename T>
  struct IsStaticallySized< Intrepid::Tensor<T> > {
    static const bool value = false;
  };

  ///
  /// Sacado specializations for Tensor3
  ///
  template <typename T>
  struct ScalarType< Intrepid::Tensor3<T> > {
    typedef T type;
  };

  template <typename T>
  struct ValueType< Intrepid::Tensor3<T> > {
    typedef T type;
  };

  template <typename T>
  struct IsADType< Intrepid::Tensor3<T> > {
    static bool const value = IsADType<T>::value;
  };

  template <typename T>
  struct IsScalarType< Intrepid::Tensor3<T> > {
    static bool const value = false;
  };

  template <typename T>
  struct Value< Intrepid::Tensor3<T> > {
    typedef typename ValueType< T >::type value_type;
    static const Intrepid::Tensor3<value_type> &
    eval(Intrepid::Tensor3<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Tensor3<value_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        for (Intrepid::Index j = 0; j < N; ++j) {
          for (Intrepid::Index k = 0; k < N; ++k) {
            v(i,j,k) = x(i,j,k).val();
          }
        }
      }

      return v;
    }
  };

  template <typename T>
  struct ScalarValue< Intrepid::Tensor3<T> > {
    typedef typename ValueType< T >::type value_type;
    typedef typename ScalarType< T >::type scalar_type;
    static const Intrepid::Tensor3<scalar_type> &
    eval(Intrepid::Tensor3<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Tensor3<scalar_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        for (Intrepid::Index j = 0; j < N; ++j) {
          for (Intrepid::Index k = 0; k < N; ++k) {
            v(i,j,k) = ScalarValue<value_type>::eval(x(i,j,k).val());
          }
        }
      }

      return v;
    }
  };

  template <typename T>
  struct StringName< Intrepid::Tensor3<T> > {
    static std::string
    eval()
    {
      return std::string("Intrepid::Tensor3< ") +
          StringName<T>::eval() + std::string(" >");
    }
  };

  template <typename T>
  struct IsEqual< Intrepid::Tensor3<T> > {
    static bool eval(T const & x, T const & y) { return x == y; }
  };

  template <typename T>
  struct IsStaticallySized< Intrepid::Tensor3<T> > {
    static const bool value = false;
  };

  ///
  /// Sacado specializations for Tensor4
  ///
  template <typename T>
  struct ScalarType< Intrepid::Tensor4<T> > {
    typedef T type;
  };

  template <typename T>
  struct ValueType< Intrepid::Tensor4<T> > {
    typedef T type;
  };

  template <typename T>
  struct IsADType< Intrepid::Tensor4<T> > {
    static bool const value = IsADType<T>::value;
  };

  template <typename T>
  struct IsScalarType< Intrepid::Tensor4<T> > {
    static bool const value = false;
  };

  template <typename T>
  struct Value< Intrepid::Tensor4<T> > {
    typedef typename ValueType< T >::type value_type;
    static const Intrepid::Tensor4<value_type> &
    eval(Intrepid::Tensor4<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Tensor4<value_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        for (Intrepid::Index j = 0; j < N; ++j) {
          for (Intrepid::Index k = 0; k < N; ++k) {
            for (Intrepid::Index l = 0; l < N; ++l) {
              v(i,j,k,l) = x(i,j,k,l).val();
            }
          }
        }
      }

      return v;
    }
  };

  template <typename T>
  struct ScalarValue< Intrepid::Tensor4<T> > {
    typedef typename ValueType< T >::type value_type;
    typedef typename ScalarType< T >::type scalar_type;
    static const Intrepid::Tensor4<scalar_type> &
    eval(Intrepid::Tensor4<T> const & x)
    {
      Intrepid::Index const
      N = x.get_dimension();

      Intrepid::Tensor4<scalar_type>
      v(N);

      for (Intrepid::Index i = 0; i < N; ++i) {
        for (Intrepid::Index j = 0; j < N; ++j) {
          for (Intrepid::Index k = 0; k < N; ++k) {
            for (Intrepid::Index l = 0; l < N; ++l) {
              v(i,j,k,l) = ScalarValue<value_type>::eval(x(i,j,k,l).val());
            }
          }
        }
      }

      return v;
    }
  };

  template <typename T>
  struct StringName< Intrepid::Tensor4<T> > {
    static std::string
    eval()
    {
      return std::string("Intrepid::Tensor4< ") +
          StringName<T>::eval() + std::string(" >");
    }
  };

  template <typename T>
  struct IsEqual< Intrepid::Tensor4<T> > {
    static bool eval(T const & x, T const & y) { return x == y; }
  };

  template <typename T>
  struct IsStaticallySized< Intrepid::Tensor4<T> > {
    static const bool value = false;
  };

} // namespace Sacado

#endif // Intrepid_MiniTensor_Definitions_h
