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
  if (ierr != 0) {                                                      \
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
    
  };
  
  /// \class Partition
  /// \brief Matrix partition parameters.
  ///
  class Partition {
  public:
    static const int Top         = 101;
    static const int Bottom      = 102;

    static const int Left        = 201;
    static const int Right       = 202;

    static const int TopLeft     = 401;
    static const int TopRight    = 402;
    static const int BottomLeft  = 403;
    static const int BottomRight = 404;
  };

  /// \class Uplo
  /// \brief Matrix upper/lower parameters.
  ///
  class Uplo {
  public:
    static const int Upper = 501;
    static const int Lower = 502;
  };

  /// \class Side
  /// \brief Matrix left/right parameters.
  class Side {
  public:
    static const int Left  = 601;
    static const int Right = 602;
  };

  /// \class Diag
  /// \brief Matrix unit/non-unit diag parameters.
  class Diag {
  public:
    static const int Unit    = 701;
    static const int NonUnit = 702;
  };

  /// \class Trans
  /// \brief Matrix upper/lower parameters.
  class Trans {
  public:
    static const int Transpose     = 801;
    static const int ConjTranspose = 802;
    static const int NoTranspose   = 803;
  };

  /// \class Variant
  /// \brief Algorithmic variants.
  class Variant {
  public:
    static const int One   = 1;
    static const int Two   = 2;
    static const int Three = 3;
    static const int Four  = 4;
  };

  /// \class AlgoChol
  /// \brief Various Cholesky algorithms for sparse and dense factorization.
  class AlgoChol {
  public:
    // - Flat sparse matrix
    static const int Dummy                  = 1000;
    static const int Unblocked              = 1001;
    static const int ExternalPardiso        = 1002;

    // - Block sparse matrix
    static const int RightLookByBlocks      = 1101;
    static const int ByBlocks               = RightLookByBlocks;

    // - Flat matrix
    static const int DenseTeamParallel      = 1201;
    static const int ExternalLapack         = 1202;

    // - Block matrix
    static const int DenseRightLookByBlocks = 1211;
    static const int DenseByBlocks          = DenseRightLookByBlocks;
  };

  /// \class AlgoBlas
  /// \brief Various matrix BLAS algorithms for sparse and dense operations.
  class AlgoBlas {
  public:
    static const int ForFactorBlocked   = 2001;
    static const int ForTriSolveBlocked = 2011;
    static const int ExternalBlas       = 2021;
    static const int InternalBlas       = 2022;
  };
  class AlgoGemm : public AlgoBlas {
  public:
    static const int DenseByBlocks = 2101;
  };

  class AlgoTrsm : public AlgoBlas {
  public:
    static const int DenseByBlocks = 2201;
  };

  class AlgoHerk : public AlgoBlas {
  public:
    static const int DenseByBlocks = 2301;
  };

}

#endif
