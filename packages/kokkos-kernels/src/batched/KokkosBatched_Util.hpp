#ifndef __KOKKOSBATCHED_UTIL_HPP__
#define __KOKKOSBATCHED_UTIL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

//#define __KOKKOSBATCHED_PROMOTION__ 1

#include <iomanip>
#include <random>
#include <string>

#include <limits>
#include <cmath>
#include <ctime>

#include <complex>

#include "Kokkos_Core.hpp"
#include "Kokkos_Complex.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosKernels_config.h"

namespace KokkosBatched {
  namespace Experimental {

    // TPL macros
#if defined (KOKKOSKERNELS_ENABLE_TPL_MKL) 
#define __KOKKOSBATCHED_INTEL_MKL__ 1
#include "mkl_version.h"
#if __INTEL_MKL__ >= 2018
#define __KOKKOSBATCHED_INTEL_MKL_BATCHED__ 1    
#define __KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__ 1
#include "mkl.h"
//#include "mkl_types.h"
#endif
#endif

#if defined (KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#define __KOKKOSKERNELS_NVIDIA_CUBLAS__ 1
#endif

#define Int2StringHelper(A) #A
#define Int2String(A) Int2StringHelper(A)
#define StringCat(A,B) A B

    void print_compiler_info();

    template<typename T> struct is_vector : public std::false_type {};

    template<typename Ta, typename Tb>
    struct is_same_mag_type {
      static const bool is_specialized = ( Kokkos::Details::ArithTraits<Ta>::is_specialized &&
                                           Kokkos::Details::ArithTraits<Tb>::is_specialized );
      
      static const bool is_mag_type_same = std::is_same<typename Kokkos::Details::ArithTraits<Ta>::mag_type,
                                                        typename Kokkos::Details::ArithTraits<Tb>::mag_type>::value;
      
      static const bool value = is_specialized && is_mag_type_same;
    };
    
    // to use double, std::complex<double>, Kokkos::complex<double>
    using std::abs;
    using std::min;
    using std::max;

    // view manipulation
    template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
    using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::Unmanaged |
					      MemoryTraitsType::RandomAccess |
					      //  MemoryTraitsType::Atomic |
					      flag>;

    template <typename ViewType>
    using UnmanagedViewType
    = Kokkos::View<typename ViewType::data_type,
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

    template <typename ViewType>
    using ScratchViewType
    = Kokkos::View<typename ViewType::data_type,
		   typename ViewType::array_layout,
		   typename ViewType::execution_space::scratch_memory_space,
		   MemoryTraits<typename ViewType::memory_traits, Kokkos::Unmanaged> >;


    // helper for vector type
    template<typename T>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if<std::is_fundamental<T>::value,size_t>::type
    adjustDimension(const size_t &m) {
      return m;
    }

    template<typename T>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!std::is_fundamental<T>::value,size_t>::type
    adjustDimension(const size_t &m) {
      return (m/T::vector_length + (m%T::vector_length > 0));
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

    template<typename T, typename dummy = T>
    struct Random;

    template<typename T>
    struct Random<T, typename std::enable_if<std::is_same<T,double>::value ||
                                             std::is_same<T,float>::value, T>::type> {
      Random(const unsigned int seed = 0) { srand(seed); }
      T value() { 
        const auto val = (rand()/((T) RAND_MAX) - 0.5)*2.0;
        return val > 0 ? val + 1.0e-3 : val - 1.0e-3;
      }
    };

    template<typename T>
    struct Random<T, typename std::enable_if<std::is_same<T,std::complex<float> >::value ||
                                             std::is_same<T,std::complex<double> >::value ||
                                             std::is_same<T,Kokkos::complex<float> >::value ||
                                             std::is_same<T,Kokkos::complex<double> >::value, T>::type> {
      Random(const unsigned int seed = 0) { srand(seed); }
      T value() {
        const auto rval = (rand()/((double) RAND_MAX) - 0.5)*2.0;        
        const auto ival = (rand()/((double) RAND_MAX) - 0.5)*2.0;        
	return T(rval > 0 ? rval + 1.0e-3 : rval - 1.0e-3,
                 ival > 0 ? ival + 1.0e-3 : ival - 1.0e-3);
      }
    };

    struct Timer {
      std::string _label;
      Kokkos::Impl::Timer _clock;
      Timer (const std::string label)
	: _label(label), _clock() {};

      void reset() { _clock.reset(); }
      double seconds() { return _clock.seconds(); }
      ~Timer() {
	Kokkos::fence();
	const double t = _clock.seconds();
	std::string label = _label; label.resize(24);
	std::cout << "KokkosKernels::Timer:: " << std::setw(26) << label
		  << std::setw(15) << std::scientific << t << " [sec] " << std::endl;
      }
    };

    // Implicit vectorization
    template<typename T>
    struct SIMD {
      static_assert( std::is_same<T,bool>::value                     ||
                     std::is_same<T,int>::value                      ||
                     std::is_same<T,size_t>::value                   ||
                     std::is_same<T,double>::value                   ||
		     std::is_same<T,float>::value                    ||
		     std::is_same<T,Kokkos::complex<float> >::value  ||
		     std::is_same<T,std::complex<float> >::value     ||
		     std::is_same<T,Kokkos::complex<double> >::value ||
		     std::is_same<T,std::complex<double> >::value,
		     "KokkosKernels:: Invalid SIMD<> type." );
      using value_type = T;
    };

    // Intel AVX instruction device (explicit vectorization)
    template<typename T>
    struct AVX {
      static_assert( std::is_same<T,double>::value                   ||
		     std::is_same<T,float>::value                    ||
		     std::is_same<T,Kokkos::complex<double> >::value ||
		     std::is_same<T,std::complex<double> >::value,
		     "KokkosKernels:: Invalid AVX<> type." );
      using value_type = T;
    };

    // Tags for BLAS
    struct Trans {
      struct Transpose {};
      struct NoTranspose {};
      struct ConjTranspose {};
    };

    struct Side {
      struct Left {};
      struct Right {};
    };

    struct Uplo {
      struct Upper {};
      struct Lower {};
    };

    struct Diag {
      struct Unit    { static const bool use_unit_diag = true;  };
      struct NonUnit { static const bool use_unit_diag = false; };
    };

    struct Algo {
      struct Level3 {
	struct Unblocked {
	  static const char* name() { return "Unblocked"; }
	};
	struct Blocked {
	  static const char* name() { return "Blocked"; }
	  // TODO:: for now harwire the blocksizes; this should reflect
	  // regieter blocking (not about team parallelism).
	  // this mb should vary according to
	  // - team policy (smaller) or range policy (bigger)
	  // - space (cuda vs host)
	  // - blocksize input (blk <= 4 mb = 2, otherwise mb = 4), etc.
#if defined(KOKKOS_ENABLE_CUDA)
	  template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION static constexpr
	  typename std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::CudaSpace>::value,int>
	  ::type mb() { return 2; }
#endif
	  template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION static constexpr
	  typename std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::HostSpace>::value,int>
	  ::type mb() { return 4; }
	};
	struct MKL {
	  static const char* name() { return "MKL"; }
	};
	struct CompactMKL {
	  static const char* name() { return "CompactMKL"; }
	};
      };

      using Gemm = Level3;
      using Trsm = Level3;
      using LU   = Level3;
      using InverseLU   = Level3;

      struct Level2 {
	struct Unblocked {};
	struct Blocked {
	  // TODO:: for now harwire the blocksizes; this should reflect
	  // regieter blocking (not about team parallelism).
	  // this mb should vary according to
	  // - team policy (smaller) or range policy (bigger)
	  // - space (cuda vs host)
	  // - blocksize input (blk <= 4 mb = 2, otherwise mb = 4), etc.
#if defined(KOKKOS_ENABLE_CUDA)
	  template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION static constexpr
	  typename std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::CudaSpace>::value,int>
	  ::type mb() { return 1; }
#endif
	  template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION static constexpr
	  typename std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::HostSpace>::value,int>
	  ::type mb() { return 4; }
	};
	struct MKL {};
	struct CompactMKL {};
      };

      using Gemv = Level2;
      using Trsv = Level2;

      //         struct Level1 {
      //           struct Unblocked {};
      //           struct Blocked {
      //             // TODO:: for now harwire the blocksizes; this should reflect
      //             // regieter blocking (not about team parallelism).
      //             // this mb should vary according to
      //             // - team policy (smaller) or range policy (bigger)
      //             // - space (cuda vs host)
      //             // - blocksize input (blk <= 4 mb = 2, otherwise mb = 4), etc.
      // #if defined(KOKKOS_ENABLE_CUDA)
      //             template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION static constexpr
      //             typename std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::CudaSpace>::value,int>
      //             ::type mb() { return 4; }
      // #endif
      //             template<typename ActiveMemorySpaceType> KOKKOS_INLINE_FUNCTION static constexpr
      //             typename std::enable_if<std::is_same<ActiveMemorySpaceType,Kokkos::HostSpace>::value,int>
      //             ::type mb() { return 4; }
      //           };
      //           //struct MKL {};
      //           //struct CompactMKL {};
      //         };

    };

    struct Util {

      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void
      packColMajor(ValueType *__restrict__ A,
		   const int m,
		   const int n,
		   const ValueType *__restrict__ B,
		   const int bs0,
		   const int bs1) {
	for (int j=0;j<n;++j)
	  for (int i=0;i<m;++i)
	    A[i+j*m] = B[i*bs0+j*bs1];
      }

      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void
      packRowMajor(ValueType *__restrict__ A,
		   const int m,
		   const int n,
		   const ValueType *__restrict__ B,
		   const int bs0,
		   const int bs1) {
	for (int i=0;i<m;++i)
	  for (int j=0;j<n;++j)
	    A[i*n+j] = B[i*bs0+j*bs1];
      }

    };

  }
}

#endif
