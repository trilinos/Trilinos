#ifndef __TACHOEXP_UTIL_HPP__
#define __TACHOEXP_UTIL_HPP__

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
#include <tuple>

#include <cmath>
#include <complex>

#include <limits>

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#ifdef TACHO_HAVE_MKL
#include "mkl.h"
#endif

/// \file TachoExp_Util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  const char* Version();

  namespace Experimental {

    ///
    /// error macros
    ///
    // #define MSG_NOT_YET_IMPLEMENTED(what) "Not yet implemented: " #what
    // #define MSG_INVALID_INPUT(what) "Invaid input argument: " #what
#define MSG_NOT_HAVE_PACKAGE(what) "Tacho does not have a package or library: " #what
#define MSG_INVALID_TEMPLATE_ARGS "Invaid template arguments"
#define MSG_INVALID_INPUT "Invaid input arguments"
#define MSG_NOT_IMPLEMENTED "Not yet implemented"

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

#if defined( KOKKOS_ENABLE_ASM ) && !defined( _WIN32 ) && !defined( __arm__ ) && !defined( __aarch64__ )
#define KOKKOS_IMPL_PAUSE asm volatile("pause\n":::"memory")
#else
#define KOKKOS_IMPL_PAUSE
#endif

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

#if defined( TACHO_USE_INT_INT )
    typedef int    ordinal_type;
    typedef int    size_type;
#elif defined( TACHO_USE_INT_SIZE_T )
    typedef int    ordinal_type;
    typedef size_t size_type;
#else
    typedef int    ordinal_type;
    typedef size_t size_type;
#endif

    ///
    /// label size used to identify object name
    ///
    enum : int { 
      LabelSize = 64,
      MaxDependenceSize = 4,
      ThresholdSolvePhaseUsingBlas3 = 12 
    };
    
    template<typename T>
    struct TypeTraits;

    template<>
    struct TypeTraits<float> {
      typedef float type;
      typedef float value_type;
      typedef float std_value_type;

      typedef float magnitude_type;
      typedef float scalar_type;
    };

    template<>
    struct TypeTraits<double> {
      typedef double type;
      typedef double value_type;
      typedef double std_value_type;

      typedef double magnitude_type;
      typedef double scalar_type;
    };

    template<>
    struct TypeTraits<std::complex<float> > {
      typedef std::complex<float> type;
      typedef std::complex<float> value_type;
      typedef std::complex<float> std_value_type;

      typedef float magnitude_type;
      typedef float scalar_type;
    };

    template<>
    struct TypeTraits<std::complex<double> > {
      typedef std::complex<double> type;
      typedef std::complex<double> value_type;
      typedef std::complex<double> std_value_type;

      typedef double magnitude_type;
      typedef double scalar_type;
    };
    
    template<>
    struct TypeTraits<Kokkos::complex<float> > {
      typedef Kokkos::complex<float> type;
      typedef Kokkos::complex<float> value_type;
      typedef std::complex<float> std_value_type;

      typedef float magnitude_type;
      typedef float scalar_type;
    };

    template<>
    struct TypeTraits<Kokkos::complex<double> > {
      typedef Kokkos::complex<double> type;
      typedef Kokkos::complex<double> value_type;
      typedef std::complex<double> std_value_type;

      typedef double magnitude_type;
      typedef double scalar_type;
    };

    // template<typename ValueType, typename ExecSpace>
    // struct DenseMatrixView;

    // template<typename ValueType, typename ExecSpace>
    // struct TypeTraits<DenseMatrixView<ValueType,ExecSpace> > {
    //   typedef DenseMatrixView<ValueType,ExecSpace> type;
    //   typedef ValueType value_type;
    //   typedef typename TypeTraits<value_type>::magnitude_type magnitude_type;
    // };

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
      return (a < static_cast<Ta>(b) ? a : static_cast<Ta>(b));
    }

    template<typename Ta, typename Tb>
    KOKKOS_FORCEINLINE_FUNCTION
    static Ta max(const Ta a, const Tb b) {
      return (a > static_cast<Ta>(b) ? a : static_cast<Ta>(b));
    }

    template<typename Ta, typename Tb>
    KOKKOS_FORCEINLINE_FUNCTION
    static void swap(Ta &a, Tb &b) {
      Ta c(a); a = static_cast<Ta>(b); b = static_cast<Tb>(c);
    }

    /// complex conj

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static T conj(const T a);

    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    double conj<double>(const double a) { 
      return a;
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    float conj<float>(const float a) { 
      return a;
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    Kokkos::complex<double> conj<Kokkos::complex<double> >(const Kokkos::complex<double> a) { 
      return Kokkos::complex<double>(a.real(), -a.imag());
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    Kokkos::complex<float> conj<Kokkos::complex<float> >(const Kokkos::complex<float> a) { 
      return Kokkos::complex<float>(a.real(), -a.imag());
    }

    /// complex real

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static typename TypeTraits<T>::scalar_type real(const T a);

    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    double real<double>(const double a) { 
      return a;
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    float real<float>(const float a) { 
      return a;
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    double real<Kokkos::complex<double> >(const Kokkos::complex<double> a) { 
      return a.real();
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    float real<Kokkos::complex<float> >(const Kokkos::complex<float> a) { 
      return a.real();
    }

    /// complex imag

    template<typename T>
    KOKKOS_FORCEINLINE_FUNCTION
    static typename TypeTraits<T>::scalar_type imag(const T a);

    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    double imag<double>(const double a) { 
      return 0.0;
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    float imag<float>(const float a) { 
      return 0.0;
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    double imag<Kokkos::complex<double> >(const Kokkos::complex<double> a) { 
      return a.imag();
    }
    
    template<>
    KOKKOS_FORCEINLINE_FUNCTION
    float imag<Kokkos::complex<float> >(const Kokkos::complex<float> a) { 
      return a.imag();
    }


    KOKKOS_FORCEINLINE_FUNCTION
    static void clear(char *buf, size_type bufsize) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      memset(buf, 0, bufsize);
#else
      for (size_type i=0;i<bufsize;++i) buf[i] = 0;
#endif
    } 
    
    template<size_t BufSize, typename SpaceType = Kokkos::DefaultExecutionSpace>
    struct Flush {
      typedef double value_type;
      
      // flush a large host buffer
      Kokkos::View<value_type*,SpaceType> _buf;
      Flush() : _buf("Flush::buf", BufSize/sizeof(double)) {
        Kokkos::deep_copy(_buf, 1);
      }

      KOKKOS_INLINE_FUNCTION
      void init(value_type &update) {
        update = 0;
      }
      
      KOKKOS_INLINE_FUNCTION
      void join(volatile value_type &update,
                const volatile value_type &input) {
        update += input;
      }
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const int i, value_type &update) const {
        update += _buf[i];
      }
      
      void run() {
        double sum = 0;
        Kokkos::parallel_reduce(Kokkos::RangePolicy<SpaceType>(0,BufSize/sizeof(double)), *this, sum);
        SpaceType::fence();
        FILE *fp = fopen("/dev/null", "w");
        fprintf(fp, "%f\n", sum);
        fclose(fp);
      }

    };

    template<typename ValueType>
    struct Random;

    template<>
    struct Random<double> {
      Random(const unsigned int seed = 0) { srand(seed); }
      double value() { return rand()/((double) RAND_MAX + 1.0); }
    };

    template<>
    struct Random<std::complex<double> > {
      Random(const unsigned int seed = 0) { srand(seed); }
      std::complex<double> value() {
        return std::complex<double>(rand()/((double) RAND_MAX + 1.0),
                                    rand()/((double) RAND_MAX + 1.0));
      }
    };

    template<>
    struct Random<Kokkos::complex<double> > {
      Random(const unsigned int seed = 0) { srand(seed); }
      Kokkos::complex<double> value() {
        return Kokkos::complex<double>(rand()/((double) RAND_MAX + 1.0),
                                       rand()/((double) RAND_MAX + 1.0));
      }
    };


    ///
    /// Tag struct
    ///
    struct NullTag { enum : int { tag = 0 }; };
    struct Partition {
      enum : int { Top = 101,
                   Bottom,
                   Left = 201,
                   Right,
                   TopLeft = 301,
                   TopRight,
                   BottomLeft,
                   BottomRight };
      
      // struct Top          { enum : int { tag = 101 }; };
      // struct Bottom       { enum : int { tag = 102 }; };
      
      // struct Left         { enum : int { tag = 201 }; };
      // struct Right        { enum : int { tag = 202 }; };
      
      // struct TopLeft      { enum : int { tag = 301 }; };
      // struct TopRight     { enum : int { tag = 302 }; };
      // struct BottomLeft   { enum : int { tag = 303 }; };
      // struct BottomRight  { enum : int { tag = 304 }; };
    };
    // template<typename T>
    // struct is_valid_partition_tag {
    //   enum : bool { value = (std::is_same<T,Partition::Top>::value        ||
    //                          std::is_same<T,Partition::Bottom>::value     ||
    //                          std::is_same<T,Partition::Left>::value       ||
    //                          std::is_same<T,Partition::Right>::value      ||
    //                          std::is_same<T,Partition::TopLeft>::value    ||
    //                          std::is_same<T,Partition::TopRight>::value   ||
    //                          std::is_same<T,Partition::BottomLeft>::value ||
    //                          std::is_same<T,Partition::BottomRight>::value)
    //   };
    // };

    struct Uplo {
      enum : int { tag = 400 };
      struct Upper        { 
        enum : int { tag = 401 }; 
        static constexpr char param = 'U'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_UPLO mkl_param = CblasUpper;
#endif
      };
      struct Lower        { 
        enum : int { tag = 402 }; 
        static constexpr char param = 'L'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_UPLO mkl_param = CblasLower;
#endif
      };
    };
    template<typename T>
    struct is_valid_uplo_tag {
      enum : bool { value = (std::is_same<T,Uplo::Upper>::value ||
                             std::is_same<T,Uplo::Lower>::value )
      };
    };
    template<typename T> struct transpose_uplo_tag;
    template<>           struct transpose_uplo_tag<Uplo::Lower> { typedef Uplo::Upper type; };
    template<>           struct transpose_uplo_tag<Uplo::Upper> { typedef Uplo::Lower type; };
    

    struct Side {
      enum : int { tag = 500 };
      struct Left         { 
        enum : int { tag = 501 }; 
        static constexpr char param = 'L'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_SIDE mkl_param = CblasLeft;
#endif
      };
      struct Right        { 
        enum : int { tag = 502 }; 
        static constexpr char param = 'R'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_SIDE mkl_param = CblasRight;
#endif
      };
    };
    template<typename T>
    struct is_valid_side_tag {
      enum : bool { value = (std::is_same<T,Side::Left>::value ||
                             std::is_same<T,Side::Right>::value )
      };
    };
    template<typename T> struct flip_side_tag;
    template<>           struct flip_side_tag<Side::Left>  { typedef Side::Right type; };
    template<>           struct flip_side_tag<Side::Right> { typedef Side::Left type; };

    struct Diag {
      enum : int { tag = 600 };
      struct Unit         { 
        enum : int { tag = 601 }; 
        static constexpr char param = 'U'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_DIAG mkl_param = CblasUnit;
#endif
      };
      struct NonUnit      { 
        enum : int { tag = 602 }; 
        static constexpr char param = 'N'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_DIAG mkl_param = CblasNonUnit;
#endif
      };
    };
    template<typename T>
    struct is_valid_diag_tag {
      enum : bool { value = (std::is_same<T,Diag::Unit>::value ||
                             std::is_same<T,Diag::NonUnit>::value )
      };
    };

    struct Trans {
      enum : int { tag = 700 };
      struct Transpose      { 
        enum : int { tag = 701 }; 
        static constexpr char param = 'T'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_TRANSPOSE mkl_param = CblasTrans;
#endif
      };
      struct ConjTranspose  { 
        enum : int { tag = 702 }; 
        static constexpr char param = 'C'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_TRANSPOSE mkl_param = CblasConjTrans;
#endif
      };
      struct NoTranspose    { 
        enum : int { tag = 703 }; 
        static constexpr char param = 'N'; 
#ifdef TACHO_HAVE_MKL  
        static constexpr CBLAS_TRANSPOSE mkl_param = CblasNoTrans;
#endif
      };
    };
    template<typename T>
    struct is_valid_trans_tag {
      enum : bool { value = (std::is_same<T,Trans::Transpose>::value ||
                             std::is_same<T,Trans::ConjTranspose>::value ||
                             std::is_same<T,Trans::NoTranspose>::value)
      };
    };
    template<typename T> struct      transpose_trans_tag;
    template<typename T> struct conj_transpose_trans_tag;
    
    template<>           struct      transpose_trans_tag<Trans::Transpose>      { typedef Trans::NoTranspose type; };
    template<>           struct      transpose_trans_tag<Trans::ConjTranspose>  { typedef Trans::NoTranspose type; };
    template<>           struct      transpose_trans_tag<Trans::NoTranspose>    { typedef Trans::Transpose type; };

    template<>           struct conj_transpose_trans_tag<Trans::Transpose>      { typedef Trans::NoTranspose type; };
    template<>           struct conj_transpose_trans_tag<Trans::ConjTranspose>  { typedef Trans::NoTranspose type; };
    template<>           struct conj_transpose_trans_tag<Trans::NoTranspose>    { typedef Trans::ConjTranspose type; };

    struct Algo {
      struct External { enum : int { tag = 1001 }; };
      struct Internal { enum : int { tag = 1002 }; };
      struct ByBlocks { enum : int { tag = 1003 }; };

      struct Workflow {
        struct Serial { enum : int { tag = 2001 }; };
      };
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
