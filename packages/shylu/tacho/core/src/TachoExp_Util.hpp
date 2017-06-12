#ifndef __TACHOEXP_UTIL_HPP__
#define __TACHOEXP_UTIL_HPP__
#define __KK__ 1

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
#include <complex>

#include <limits>

#include "Teuchos_BLAS_types.hpp"

/// \file TachoExp_Util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  namespace Experimental {

    ///
    /// error macros
    ///
    // #define MSG_NOT_YET_IMPLEMENTED(what) "Not yet implemented: " #what
    // #define MSG_INVALID_INPUT(what) "Invaid input argument: " #what
#define MSG_NOT_HAVE_PACKAGE(what) "Tacho does not have a package or library: " #what
#define MSG_INVALID_TEMPLATE_ARGS "Invaid template arguments"

#define TACHO_TEST_FOR_ABORT(ierr, msg)                                 \
    if ((ierr) != 0) {                                                  \
      printf(">> Error in file %s, line %d, error %d \n   %s\n",__FILE__,__LINE__,ierr,msg); \
      Kokkos::abort(">> Tacho abort\n");                                \
    }

#define TACHO_TEST_FOR_EXCEPTION(ierr, x, msg)                          \
    if ((ierr) != 0) {                                                  \
      fprintf(stderr, ">> Error in file %s, line %d, error %d \n",__FILE__,__LINE__,ierr); \
      fprintf(stderr, "   %s\n", msg);                                  \
      throw x(msg);                                                     \
    }

    ///
    /// print execution spaces
    ///
    template<typename SpT>
    void printExecSpaceConfiguration(std::string name, const bool detail = false) {
      TACHO_TEST_FOR_EXCEPTION(!Kokkos::Impl::is_space<SpT>::value,
                               std::logic_error,
                               "SpT is not Kokkos execution space.");
      std::cout << std::setw(16) << name << "::  ";
      if (std::is_same<SpT,Kokkos::Serial>::value)
        std::cout << "Kokkos::Serial " << std::endl;
      else
        SpT::print_configuration(std::cout, detail);
    }

    ///
    /// default ordinal and size type
    ///
    typedef int    ordinal_type;
    typedef size_t size_type;

    ///
    /// label size used to identify object name
    ///
    enum : int { LabelSize = 64,
                 MaxDependenceSize = 10 };

    ///
    /// complex query
    ///
    template<typename T>
    struct is_complex_type { enum : bool { value = false }; };

    template< typename T >
    struct is_complex_type< Kokkos::complex<T> > { enum : bool { value = true }; };

    template< typename T >
    struct is_complex_type< std::complex<T> > { enum : bool { value = true }; };

    ///
    /// Coo : Sparse coordinate format; (i, j, val).
    ///
    template<typename ValueType>
    struct Coo {
      typedef ValueType value_type;

      ordinal_type i,j;
      value_type val;

      Coo() = default;
      Coo(const ordinal_type i,
          const ordinal_type j,
          const value_type val)
        : i(i), j(j), val(val) {}
      Coo(const Coo& b) = default;

      /// \brief Compare "less" index i and j only.
      bool operator<(const Coo &y) const {
        const auto r_val = (this->i - y.i);
        return (r_val == 0 ? this->j < y.j : r_val < 0);
      }

      /// \brief Compare "equality" only index i and j.
      bool operator==(const Coo &y) const {
        return (this->i == y.i) && (this->j == y.j);
      }

      /// \brief Compare "in-equality" only index i and j.
      bool operator!=(const Coo &y) const {
        return !(*this == y);
      }
    };

    ///
    /// util
    ///
    template<typename Ta, typename Tb>
    KOKKOS_FORCEINLINE_FUNCTION
    static Ta min(const Ta a, const Tb b) {
      return (a < b ? a : b);
    }

    template<typename Ta, typename Tb>
    KOKKOS_FORCEINLINE_FUNCTION
    static Ta max(const Ta a, const Tb b) {
      return (a > b ? a : b);
    }

    template<typename Ta, typename Tb>
    KOKKOS_FORCEINLINE_FUNCTION
    static void swap(Ta &a, Tb &b) {
      Ta c(a); a = b; b = c;
    }

    ///
    /// Tag struct
    ///
    struct NullTag { enum : int { tag = 0 }; };
    struct Partition {
      struct Top          { enum : int { tag = 101 }; };
      struct Bottom       { enum : int { tag = 102 }; };

      struct Left         { enum : int { tag = 201 }; };
      struct Right        { enum : int { tag = 202 }; };

      struct TopLeft      { enum : int { tag = 301 }; };
      struct TopRight     { enum : int { tag = 302 }; };
      struct BottomLeft   { enum : int { tag = 303 }; };
      struct BottomRight  { enum : int { tag = 304 }; };
    };
    template<typename T>
    struct is_valid_partition_tag {
      enum : bool { value = (std::is_same<T,Partition::Top>::value        ||
                             std::is_same<T,Partition::Bottom>::value     ||
                             std::is_same<T,Partition::Left>::value       ||
                             std::is_same<T,Partition::Right>::value      ||
                             std::is_same<T,Partition::TopLeft>::value    ||
                             std::is_same<T,Partition::TopRight>::value   ||
                             std::is_same<T,Partition::BottomLeft>::value ||
                             std::is_same<T,Partition::BottomRight>::value)
      };
    };

    struct Uplo {
      struct Upper        { enum : int { tag = 401 }; static constexpr char param = 'U'; static constexpr int teuchos = Teuchos::UPPER_TRI; };
      struct Lower        { enum : int { tag = 402 }; static constexpr char param = 'L'; static constexpr int teuchos = Teuchos::LOWER_TRI; };
    };
    template<typename T>
    struct is_valid_uplo_tag {
      enum : bool { value = (std::is_same<T,Uplo::Upper>::value ||
                             std::is_same<T,Uplo::Lower>::value )
      };
    };

    struct Side {
      struct Left         { enum : int { tag = 501 }; static constexpr char param = 'L'; static constexpr int teuchos = Teuchos::LEFT_SIDE; };
      struct Right        { enum : int { tag = 502 }; static constexpr char param = 'R'; static constexpr int teuchos = Teuchos::RIGHT_SIDE; };
    };
    template<typename T>
    struct is_valid_side_tag {
      enum : bool { value = (std::is_same<T,Side::Left>::value ||
                             std::is_same<T,Side::Right>::value )
      };
    };

    struct Diag {
      struct Unit         { enum : int { tag = 601 }; static constexpr char param = 'U'; static constexpr int teuchos = Teuchos::UNIT_DIAG; };
      struct NonUnit      { enum : int { tag = 602 }; static constexpr char param = 'N'; static constexpr int teuchos = Teuchos::NON_UNIT_DIAG; };
    };
    template<typename T>
    struct is_valid_diag_tag {
      enum : bool { value = (std::is_same<T,Diag::Unit>::value ||
                             std::is_same<T,Diag::NonUnit>::value )
      };
    };

    struct Trans {
      struct Transpose      { enum : int { tag = 701 }; static constexpr char param = 'T'; static constexpr int teuchos = Teuchos::TRANS; };
      struct ConjTranspose  { enum : int { tag = 702 }; static constexpr char param = 'C'; static constexpr int teuchos = Teuchos::CONJ_TRANS; };
      struct NoTranspose    { enum : int { tag = 703 }; static constexpr char param = 'N'; static constexpr int teuchos = Teuchos::NO_TRANS; };
    };
    template<typename T>
    struct is_valid_trans_tag {
      enum : bool { value = (std::is_same<T,Trans::Transpose>::value ||
                             std::is_same<T,Trans::ConjTranspose>::value ||
                             std::is_same<T,Trans::NoTranspose>::value)
      };
    };

    struct Algo {
      struct External { enum : int { tag = 1001 }; };
      struct SuperNodes { enum : int { tag = 1002 }; };
    };

    template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
    using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::Unmanaged |
                                              MemoryTraitsType::RandomAccess |
                                              MemoryTraitsType::Atomic |
                                              flag>;

    template <typename ViewType>
    using UnmanagedViewType = Kokkos::View<typename ViewType::data_type,
                                           typename ViewType::array_layout,
                                           typename ViewType::device_type,
                                           MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged> >;
    template <typename ViewType>
    using ConstViewType = Kokkos::View<typename ViewType::const_data_type,
                                       typename ViewType::array_layout,
                                       typename ViewType::device_type,
                                       typename ViewType::memory_traits>;
    template <typename ViewType>
    using ConstUnmanagedViewType = ConstViewType<UnmanagedViewType<ViewType> >;

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if<!std::is_integral<T>::value,ordinal_type>::type
    get_team_rank(const T& member) {
      return member.team_rank();
    }

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    typename std::enable_if<std::is_integral<T>::value,ordinal_type>::type
    get_team_rank(const T& member) {
      return 0;
    }

  }
}

#endif
