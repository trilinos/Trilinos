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
    printf(">> Error in file %s, line %d, error %d \n   %s\n",__FILE__,__LINE__,ierr,msg); \
    Kokkos::abort(">> Tacho abort\n");                                  \
  }

#define TACHO_TEST_FOR_EXCEPTION(ierr, x, msg)                           \
  if ((ierr) != 0) {                                                    \
    fprintf(stderr, ">> Error in file %s, line %d, error %d \n",__FILE__,__LINE__,ierr); \
    fprintf(stderr, "   %s\n", msg);                                    \
    throw x(msg);                                                       \
  }

#define TACHO_TEST_FOR_WARNING(ierr, msg)                                 \
  if ((ierr) != 0) {                                                    \
    printf(">> Warning in file %s, line %d, error %d \n   %s\n",__FILE__,__LINE__,ierr,msg); \
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

  // default value
  template<typename T>
  struct is_complex_type { enum : bool { value = false }; };

  // specialization
  template< typename T >
  struct is_complex_type< Kokkos::complex<T> >
    { enum : bool { value = true }; };

  // default value
  template<typename T>
  struct is_scalar_type { enum : bool { value = false }; };

  template<>
  struct is_scalar_type<int> { enum : bool { value = true }; };
  template<>
  struct is_scalar_type<unsigned int> { enum : bool { value = true }; };
  template<>
  struct is_scalar_type<long> { enum : bool { value = true }; };
  template<>
  struct is_scalar_type<unsigned long> { enum : bool { value = true }; };
  template<>
  struct is_scalar_type<float> { enum : bool { value = true }; };
  template<>
  struct is_scalar_type<double> { enum : bool { value = true }; };
  template<>
  struct is_scalar_type<Kokkos::complex<float> > { enum : bool { value = true }; };
  template<>
  struct is_scalar_type<Kokkos::complex<double> > { enum : bool { value = true }; };


  class Util {
  public:
    static constexpr size_t LabelSize = 64;

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

    // uses range [first, last)
    template<typename ValueType, typename SpaceType, typename IterType>
    KOKKOS_FORCEINLINE_FUNCTION
    static IterType getLowerBound(const Kokkos::View<ValueType*,SpaceType> &data,
                                  IterType first,
                                  IterType last,
                                  const ValueType val) {
      IterType it, count = last - first, step = 0;
      while (count > 0) {
        it = first;
        it += ( step = (count >> 1) );
        if (data[it] < val) {
          first = ++it;
          count -= step+1;
        } else {
          count=step;
        }
      }
      return first;
    }

    template<typename ValueType, typename SpaceType, typename IterType>
    KOKKOS_FORCEINLINE_FUNCTION
    static IterType getUpperBound(const Kokkos::View<ValueType*,SpaceType> &data,
                                  IterType first,
                                  IterType last,
                                  const ValueType val) {
      IterType it, count = last - first, step;
      while (count > 0) {
        it = first;
        it += ( step = (count >> 1) );
        if (!(val < data[it])) {
          first = ++it;
          count -= step+1;
        } else {
          count = step;
        }
      }
      return first;
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static void swap(T &a, T &b) {
      T c(a); a = b; b = c;
    }

    template<typename dataValueType, class ... dataProperties,
             typename idxValueType, class ... idxProperties,
             typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void sort(/**/  Kokkos::View<dataValueType*,dataProperties...> data,
                     /**/  Kokkos::View<idxValueType*,idxProperties...> idx,
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

    template<typename dataValueType, class ...dataProperties,
             typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void sort(/**/  Kokkos::View<dataValueType*,dataProperties...> data,
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
    static constexpr int Five  = 5;
    static constexpr int Six   = 6;
    static constexpr int Seven = 7;
    static constexpr int Eight = 8;
    static constexpr int Nine  = 9;
  };

  class TaskWindow {
  public:
    static constexpr unsigned int CholByBlocks     = 4096;
    static constexpr unsigned int TriSolveByBlocks = 4096;
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

  /// \class AlgoTriSolve
  /// \brief Various Cholesky algorithms for sparse and dense factorization.
  class AlgoTriSolve {
  public:
    // - Flat sparse matrix
    static constexpr int Dummy                  = 1400;
    static constexpr int Unblocked              = 1401;
    static constexpr int ExternalPardiso        = 1402;

    // - Block sparse matrix
    static constexpr int ByBlocks               = 1501;

    // - Flat sparse with nested dense matrices
    static constexpr int SuperNodes             = 1601;
    static constexpr int SuperNodesByBlocks     = 1602;
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


  struct Stat {
    double flop;

    Stat() : flop(0.0) {}
    Stat& operator+=(const Stat &b) {
      flop += b.flop;
      return *this;
    }
  };
  inline Stat operator+(Stat a, const Stat &b) {
    return a += b;
  }
  
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
