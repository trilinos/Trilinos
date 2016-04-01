#ifndef __TACHO_UTIL_HPP__
#define __TACHO_UTIL_HPP__

// standard C includes
#include <stdio.h>
#include <string.h>

// "std" includes
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <memory>

#include <cmath>
//#include <complex>

#include <limits>

/// \file Tacho_Util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  /// \brief Error handling.
  //
#define MSG_NOT_YET_IMPLEMENTED "Not yet implemented"
#define MSG_INVALID_INPUT(what) "Invaid input argument: " #what
#define MSG_NOT_HAVE_PACKAGE(what) "Tacho does not have a package or library: " what
#define MSG_INVALID_TEMPLATE_ARGS "Invaid template arguments"

#define TACHO_TEST_FOR_ABORT(ierr, msg)                                 \
  if ((ierr) != 0) {                                                    \
    fprintf(stderr, ">> Error in file %s, line %d, error %d \n",__FILE__,__LINE__,ierr); \
    fprintf(stderr, "   %s\n", msg);                                    \
    Kokkos::abort(">> Tacho abort\n");                                    \
  }

  /// \brief Control parameter decomposition.
  ///

  // control id
#undef  Ctrl
#define Ctrl(name,algo,variant) name<algo,variant>

  // control leaf
#undef CtrlComponent
#define CtrlComponent(name,algo,variant,component,id)   \
  Ctrl(name,algo,variant)::component[id]

  // control recursion
#undef CtrlDetail
#define CtrlDetail(name,algo,variant,component)                         \
  CtrlComponent(name,algo,variant,component,0),CtrlComponent(name,algo,variant,component,1),name

  template<typename T>
  struct is_complex_type {
    static bool value;
  };

  // default value
  template<typename T> bool is_complex_type<T>::value = false;

  // specialization
  template<> bool is_complex_type<float>::value = false;
  template<> bool is_complex_type<double>::value = false;
  template<> bool is_complex_type<Kokkos::complex<float> >::value = true;
  template<> bool is_complex_type<Kokkos::complex<double> >::value = true;

  template<typename T>
  struct is_scalar_type {
    static bool value;
  };

  // default value
  template<typename T> bool is_scalar_type<T>::value = false;

  template<> bool is_scalar_type<int>::value = true;
  template<> bool is_scalar_type<unsigned int>::value = true;
  template<> bool is_scalar_type<long>::value = true;
  template<> bool is_scalar_type<size_t>::value = true;
  template<> bool is_scalar_type<float>::value = true;
  template<> bool is_scalar_type<double>::value = true;
  template<> bool is_scalar_type<Kokkos::complex<float> >::value = true;
  template<> bool is_scalar_type<Kokkos::complex<double> >::value = true;

  class Util {
  public:
    static const size_t LabelSize = 64;

    template<typename T>
    KOKKOS_INLINE_FUNCTION
    static T min(const T a, const T b) {
      return (a < b ? a : b);
    }

    template<typename T>
    KOKKOS_INLINE_FUNCTION
    static T max(const T a, const T b) {
      return (a > b ? a : b);
    }

    template<typename T>
    KOKKOS_INLINE_FUNCTION
    static T abs(const T a) {
      return (a > 0 ? a : -a);
    }

    template<typename T>
    KOKKOS_INLINE_FUNCTION
    static T real(const T a) {
      return a;
    }

    template<typename T>
    KOKKOS_INLINE_FUNCTION
    static T imag(const T a) {
      return 0;
    }

    template<typename T>
    KOKKOS_INLINE_FUNCTION
    static T conj(const T a) {
      return a;
    }

    template<typename SpT>
    KOKKOS_INLINE_FUNCTION
    static size_t alignDimension(size_t dim, size_t value_type_size) {
      return dim;
    }

    template<typename T1, typename T2, typename T3>
    KOKKOS_FORCEINLINE_FUNCTION
    static void unrollIndex(Kokkos::pair<T1,T1> &idx,
                            const T2 k, const T3 stride) {
      idx.first  = k%stride;
      idx.second = k/stride;
    }

    template<typename T1, typename T2, typename T3, typename T4>
    KOKKOS_FORCEINLINE_FUNCTION
    static void unrollIndex(T1 &i, T2 &j,
                            const T3 k,
                            const T4 stride) {
      i = k%stride;
      j = k/stride;
    }

    template<size_t N, typename Lambda, typename IterT>
    KOKKOS_FORCEINLINE_FUNCTION
    static void unrollLoop(const Lambda &f, const IterT& iter) {
      if (N != 0) unrollLoop<N-1>(f, iter);
      f(iter + N);
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static bool isComplex() {
      return is_complex_type<T>::value;
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static bool isScalar() {
      return is_scalar_type<T>::value;
    }

    template<typename ValueType, typename SpaceType, typename IterType>
    KOKKOS_FORCEINLINE_FUNCTION
    static IterType getLowerBound(const Kokkos::View<ValueType*,SpaceType> &data,
                                  const IterType begin,
                                  const IterType end,
                                  const ValueType val) {
      IterType it = begin, tmp = begin, count = end - begin, step = 0;
      while (count > 0) {
        it = tmp ;
        it += ( step = (count >> 1) );
        if (data[it] < val) {
          tmp = ++it;
          count -= (step+1);
        } else {
          count=step;
        }
      }
      return tmp;
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static void swap(T &a, T &b) {
      T c(a); a = b; b = c;
    }

    template<typename ValueType,
             typename IndexType,
             typename OrdinalType,
             typename SpaceType>
    KOKKOS_INLINE_FUNCTION
    static void sort(Kokkos::View<ValueType*,SpaceType> data,
                     Kokkos::View<IndexType*,SpaceType> idx,
                     const OrdinalType begin,
                     const OrdinalType end) {
      if (begin + 1 < end) {
        const auto piv = data[begin];
        OrdinalType left = (begin + 1), right = end;
        while (left < right) {
          if (data[left] <= piv) {
            ++left;
          } else {
            --right;
            Util::swap(data[left], data[right]);
            Util::swap(idx [left], idx [right]);
          }
        }

        --left;
        Util::swap(data[left], data[begin]);
        Util::swap(idx [left], idx [begin]);

        // recursion
        Util::sort(data, idx, begin, left);
        Util::sort(data, idx, right, end );
      }
    }

    template<typename ValueType,
             typename OrdinalType,
             typename SpaceType>
    KOKKOS_INLINE_FUNCTION
    static void sort(Kokkos::View<ValueType*,SpaceType> data,
                     const OrdinalType begin,
                     const OrdinalType end) {
      if (begin + 1 < end) {
        const auto piv = data[begin];
        OrdinalType left = (begin + 1), right = end;
        while (left < right) {
          if (data[left] <= piv) {
            ++left;
          } else {
            --right;
            Util::swap(data[left], data[right]);
          }
        }

        --left;
        Util::swap(data[left], data[begin]);

        // recursion
        Util::sort(data, begin, left);
        Util::sort(data, right, end );
      }
    }

  };

  /// \class Partition
  /// \brief Matrix partition parameters.
  ///
  class Partition {
  public:
    static constexpr int Top         = 101;
    static constexpr int Bottom      = 102;

    static constexpr int Left        = 201;
    static constexpr int Right       = 202;

    static constexpr int TopLeft     = 401;
    static constexpr int TopRight    = 402;
    static constexpr int BottomLeft  = 403;
    static constexpr int BottomRight = 404;
  };

  /// \class Uplo
  /// \brief Matrix upper/lower parameters.
  ///
  class Uplo {
  public:
    static constexpr int Upper = 501;
    static constexpr int Lower = 502;
  };

  /// \class Side
  /// \brief Matrix left/right parameters.
  class Side {
  public:
    static constexpr int Left  = 601;
    static constexpr int Right = 602;
  };

  /// \class Diag
  /// \brief Matrix unit/non-unit diag parameters.
  class Diag {
  public:
    static constexpr int Unit    = 701;
    static constexpr int NonUnit = 702;
  };

  /// \class Trans
  /// \brief Matrix upper/lower parameters.
  class Trans {
  public:
    static constexpr int Transpose     = 801;
    static constexpr int ConjTranspose = 802;
    static constexpr int NoTranspose   = 803;
  };

  /// \class Variant
  /// \brief Algorithmic variants.
  class Variant {
  public:
    static constexpr int One   = 1;
    static constexpr int Two   = 2;
    static constexpr int Three = 3;
    static constexpr int Four  = 4;
  };

  /// \class AlgoChol
  /// \brief Various Cholesky algorithms for sparse and dense factorization.
  class AlgoChol {
  public:
    // - Flat sparse matrix
    static constexpr int Dummy                  = 1000;
    static constexpr int Unblocked              = 1001;
    static constexpr int ExternalPardiso        = 1002;

    // - Block sparse matrix
    static constexpr int ByBlocks               = 1101;

    // - Flat dense matrix
    static constexpr int ExternalLapack         = 1202;

    // - Hier dense matrix
    static constexpr int DenseByBlocks          = 1211;

    // - Flat sparse with nested dense matrices
    static constexpr int SuperNodes             = 1301;
    static constexpr int SuperNodesByBlocks     = 1302;
  };

  /// \class AlgoBlas
  /// \brief Various matrix BLAS algorithms for sparse and dense operations.
  class AlgoBlas {
  public:
    // - Flat sparse-sparse matrix
    static constexpr int SparseSparseUnblocked          = 2001;
    static constexpr int SparseDenseUnblocked           = 2002;

    // - Flat dense matrix
    static constexpr int ExternalBlas                   = 2011;
    static constexpr int InternalBlas                   = 2012;

    // - Flat sparse with nested dense matrices
    static constexpr int SparseSparseSuperNodes         = 2021;
    static constexpr int SparseDenseSuperNodes          = 2022;
    static constexpr int SparseSparseSuperNodesByBlocks = 2033;
    static constexpr int SparseDenseSuperNodesByBlocks  = 2034;
  };

  class AlgoGemm : public AlgoBlas {
  public:
    static constexpr int DenseByBlocks = 2101;
  };

  class AlgoTrsm : public AlgoBlas {
  public:
    static constexpr int DenseByBlocks = 2201;
  };

  class AlgoHerk : public AlgoBlas {
  public:
    static constexpr int DenseByBlocks = 2301;
  };

  /// \class Coo
  /// \brief Sparse coordinate format; (i, j, val).
  template<typename OrdinalType, typename ValueType>
  class Coo {
  public:
    typedef OrdinalType ordinal_type;
    typedef ValueType   value_type;

  public:
    ordinal_type _i,_j;
    value_type _val;

  public:
    ordinal_type& Row() { return _i;   }
    ordinal_type& Col() { return _j;   }
    value_type&   Val() { return _val; }

    ordinal_type  Row() const { return _i;   }
    ordinal_type  Col() const { return _j;   }
    value_type    Val() const { return _val; }

    Coo() = default;
    Coo(const ordinal_type i,
        const ordinal_type j,
        const value_type val)
      : _i(i), _j(j), _val(val) {}
    Coo(const Coo& b) = default;

    /// \brief Compare "less" index i and j only.
    bool operator<(const Coo &y) const {
      ordinal_type r_val = (this->_i - y._i);
      return (r_val == 0 ? this->_j < y._j : r_val < 0);
    }

    /// \brief Compare "equality" only index i and j.
    bool operator==(const Coo &y) const {
      return (this->_i == y._i) && (this->_j == y._j);
    }

    /// \brief Compare "in-equality" only index i and j.
    bool operator!=(const Coo &y) const {
      return !(*this == y);
    }
  };



}

#endif
