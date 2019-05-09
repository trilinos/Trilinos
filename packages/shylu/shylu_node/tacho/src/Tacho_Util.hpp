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
#include <tuple>

#include <cmath>
#include <complex>

#include <limits>

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#if defined (__INTEL_MKL__)
#include "mkl.h"
#endif

/// \file Tacho_Util.hpp
/// \brief Utility functions and constant integer class like an enum class.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  const char* Version();

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


    #if defined( KOKKOS_ENABLE_ASM )
      #if defined( __amd64 )  || defined( __amd64__ ) || \
          defined( __x86_64 ) || defined( __x86_64__ )
          #if !defined( _WIN32 ) /* IS NOT Microsoft Windows */
             #define KOKKOS_IMPL_PAUSE asm volatile( "pause\n":::"memory" );
          #else
             #define KOKKOS_IMPL_PAUSE __asm__ __volatile__( "pause\n":::"memory" );
          #endif
      #elif defined(__PPC64__)
            #define KOKKOS_IMPL_PAUSE  asm volatile( "or 27, 27, 27" ::: "memory" );
      #endif
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
      SpT::print_configuration(std::cout, detail);
    }

    ///
    /// default ordinal and size type
    ///

#if defined( TACHO_USE_INT_INT )
    typedef int ordinal_type;
    typedef int size_type;
#elif defined( TACHO_USE_INT_SIZE_T )
    typedef int ordinal_type;
    typedef size_t size_type;
#else
    typedef int ordinal_type;
    typedef size_t size_type;
#endif

    ///
    /// label size used to identify object name
    ///
    enum : int { 
      MaxDependenceSize = 4,
      ThresholdSolvePhaseUsingBlas3 = 12,
      CudaVectorSize = 4
    };

    ///
    /// later, this would be replaced with Kokkos::ArithTraits
    ///
    template<typename T>
    struct ArithTraits;

    template<>
    struct ArithTraits<float> {
      typedef float val_type;
      typedef float mag_type;

      enum : bool { is_complex = false };
      static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return x > 0 ? x : -x; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x; }
      static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return x; }
    };

    template<>
    struct ArithTraits<double> {
      typedef double val_type;
      typedef double mag_type;

      enum : bool { is_complex = false };
      static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return x > 0 ? x : -x; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x; }
      static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return x; }
    };

    template<>
    struct ArithTraits<std::complex<float> > {
      typedef std::complex<float> val_type;
      typedef float mag_type;

      enum : bool { is_complex = true };
      static inline mag_type abs (const val_type& x) { return std::abs(x); }
      static inline mag_type real(const val_type& x) { return x.real(); }
      static inline mag_type imag(const val_type& x) { return x.imag(); }
      static inline val_type conj(const val_type& x) { return std::conj(x); }
    };

    template<>
    struct ArithTraits<std::complex<double> > {
      typedef std::complex<double> val_type;
      typedef double mag_type;

      enum : bool { is_complex = true };
      static inline mag_type abs (const val_type& x) { return std::abs(x); }
      static inline mag_type real(const val_type& x) { return x.real(); }
      static inline mag_type imag(const val_type& x) { return x.imag(); }
      static inline val_type conj(const val_type& x) { return std::conj(x); }
    };
    
    template<>
    struct ArithTraits<Kokkos::complex<float> > {
      typedef Kokkos::complex<float> val_type;
      typedef float mag_type;

      enum : bool { is_complex = true };
      static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return Kokkos::abs(x); }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x.real(); }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x.imag(); }
      static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return Kokkos::conj(x); }
    };

    template<>
    struct ArithTraits<Kokkos::complex<double> > {
      typedef Kokkos::complex<double> val_type;
      typedef double mag_type;

      enum : bool { is_complex = true };
      static KOKKOS_FORCEINLINE_FUNCTION mag_type abs (const val_type& x) { return Kokkos::abs(x); }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type& x) { return x.real(); }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type& x) { return x.imag(); }
      static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type& x) { return Kokkos::conj(x); }
    };

    ///
    /// Coo : Sparse coordinate format; (i, j, val).
    ///
    template<typename ValueType>
    struct Coo {
      typedef ValueType value_type;

      ordinal_type i,j;
      value_type val;

      Coo() = default;
      Coo(const ordinal_type ii,
          const ordinal_type jj,
          const value_type vval)
        : i(ii), j(jj), val(vval) {}
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

    KOKKOS_FORCEINLINE_FUNCTION
    static void clear(char *buf, size_type bufsize) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      memset(buf, 0, bufsize);
#else
      for (size_type i=0;i<bufsize;++i) buf[i] = 0;
#endif
    } 
    
    template<typename MemberType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void clear(MemberType &member, char *buf, size_type bufsize) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      memset(buf, 0, bufsize);
#else
      const ordinal_type team_index_range = (bufsize/CudaVectorSize) + (bufsize%CudaVectorSize > 0);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,team_index_range),[&](const int &idx) {
          const int ioff = idx * CudaVectorSize;
          const int itmp = bufsize - ioff;
          const int icnt = itmp > CudaVectorSize ? CudaVectorSize : itmp;
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,icnt),[&](const int &ii) {
              const int i = ioff + ii;
              buf[i] = 0;
            });
        });
#endif
    } 

    template<typename T1, typename T2, typename CompareType>
    KOKKOS_INLINE_FUNCTION
    static T1* lower_bound(T1* first, T1* last, const T2& val, 
                           CompareType compare) {
      T1 *it;
      ordinal_type step = 0, count = last - first;
      while (count > 0) {
        it = first; step = count/2; it += step;
        if (compare(*it,val)) {
          first = ++it;
          count -= step + 1;
        } else {
          count = step;
        }
      }
      return first;
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
    };

    struct Uplo {
      enum : int { tag = 400 };
      struct Upper        { 
        enum : int { tag = 401 }; 
        static constexpr char param = 'U'; 
#if defined(__INTEL_MKL__)  
        static constexpr CBLAS_UPLO mkl_param = CblasUpper;
#endif
      };
      struct Lower        { 
        enum : int { tag = 402 }; 
        static constexpr char param = 'L'; 
#if defined(__INTEL_MKL__)  
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
#if defined(__INTEL_MKL__)  
        static constexpr CBLAS_SIDE mkl_param = CblasLeft;
#endif
      };
      struct Right        { 
        enum : int { tag = 502 }; 
        static constexpr char param = 'R'; 
#if defined(__INTEL_MKL__)  
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
#if defined(__INTEL_MKL__)  
        static constexpr CBLAS_DIAG mkl_param = CblasUnit;
#endif
      };
      struct NonUnit      { 
        enum : int { tag = 602 }; 
        static constexpr char param = 'N'; 
#if defined(__INTEL_MKL__)  
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
#if defined(__INTEL_MKL__)  
        static constexpr CBLAS_TRANSPOSE mkl_param = CblasTrans;
#endif
      };
      struct ConjTranspose  { 
        enum : int { tag = 702 }; 
        static constexpr char param = 'C'; 
#if defined(__INTEL_MKL__)  
        static constexpr CBLAS_TRANSPOSE mkl_param = CblasConjTrans;
#endif
      };
      struct NoTranspose    { 
        enum : int { tag = 703 }; 
        static constexpr char param = 'N'; 
#if defined(__INTEL_MKL__)  
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

    ///
    /// helper functions
    ///
    struct Conjugate {
      enum : int { tag = 801 };
      
      KOKKOS_FORCEINLINE_FUNCTION Conjugate() {} 
      KOKKOS_FORCEINLINE_FUNCTION Conjugate(const Conjugate &b) {} 

      KOKKOS_FORCEINLINE_FUNCTION float operator()(const float &v) const { return v; }
      KOKKOS_FORCEINLINE_FUNCTION double operator()(const double &v) const { return v; }
      inline std::complex<float> operator()(const std::complex<float> &v) const { return std::conj(v); }
      inline std::complex<double> operator()(const std::complex<double> &v) const { return std::conj(v); }
      KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<float> operator()(const Kokkos::complex<float> &v) const { return Kokkos::conj(v); }
      KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<double> operator()(const Kokkos::complex<double> &v) const { return Kokkos::conj(v); }
    };
    
    struct NoConjugate {
      enum : int {tag = 802 };

      KOKKOS_FORCEINLINE_FUNCTION NoConjugate() {}
      KOKKOS_FORCEINLINE_FUNCTION NoConjugate(const NoConjugate &b) {}
      
      KOKKOS_FORCEINLINE_FUNCTION float operator()(const float &v) const { return v; }
      KOKKOS_FORCEINLINE_FUNCTION double operator()(const double &v) const { return v; }
      inline std::complex<float> operator()(const std::complex<float> &v) const { return v; }
      inline std::complex<double> operator()(const std::complex<double> &v) const { return v; }
      KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<float> operator()(const Kokkos::complex<float> &v) const { return v; }
      KOKKOS_FORCEINLINE_FUNCTION Kokkos::complex<double> operator()(const Kokkos::complex<double> &v) const { return v; }
    };

    struct Algo {
      struct External { enum : int { tag = 1001 }; };
      struct Internal { enum : int { tag = 1002 }; };
      struct ByBlocks { enum : int { tag = 1003 }; };

      struct Workflow {
        struct Serial      { enum : int { tag = 2001 }; };
        struct SerialPanel { enum : int { tag = 2002 }; };
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

}

#endif
