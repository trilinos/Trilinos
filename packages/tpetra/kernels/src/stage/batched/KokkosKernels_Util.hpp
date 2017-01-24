#ifndef __KOKKOSKERNELS_UTIL_HPP__
#define __KOKKOSKERNELS_UTIL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>
#include <random>
#include <string>
#include <ctime>

#include <complex>
#include "Kokkos_Complex.hpp"

namespace KokkosKernels {

  struct Random {
    Random(const unsigned int seed = 0) { srand(seed); }
    double value() { return rand()/((double) RAND_MAX + 1.0); }
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
      std::cout << "KokkosKernels::Timer:: " << std::setw(12) << _label 
                << std::setw(10) << std::scientific << t << " [sec] " << std::endl;
    }
  };
  
  // Intel, IBM SMT device (implicit vectorization)
  template<typename T, 
           typename SpT = Kokkos::DefaultExecutionSpace>
  struct SIMD { 
    static_assert( std::is_same<T,double>::value                   ||
                   std::is_same<T,float>::value                    ||
                   std::is_same<T,Kokkos::complex<double> >::value || 
                   std::is_same<T,std::complex<double> >::value,
                   "KokkosKernels:: Invalid SIMD<> type." );

    static_assert( std::is_same<SpT,Kokkos::Serial>::value ||
                   std::is_same<SpT,Kokkos::OpenMP>::value,
                   "KokkosKernels:: Invalid SIMD<> exec space." );

    using value_type = T; 
    using exec_space = SpT;
  };
  
  // Intel AVX instruction device (explicit vectorization)
  template<typename T,
           typename SpT = Kokkos::DefaultExecutionSpace>
  struct AVX {
    static_assert( std::is_same<T,double>::value                   ||
                   std::is_same<T,float>::value                    ||
                   std::is_same<T,Kokkos::complex<double> >::value || 
                   std::is_same<T,std::complex<double> >::value,
                   "KokkosKernels:: Invalid AVX<> type." );

    static_assert( std::is_same<SpT,Kokkos::Serial>::value ||
                   std::is_same<SpT,Kokkos::OpenMP>::value,
                   "KokkosKernels:: Invalid AVX<> exec space." );

    using value_type = T;
    using exec_space = SpT;
  };

  // Cuda threading (explicit thread vectorization)
  template<typename T, 
           typename SpT = Kokkos::DefaultExecutionSpace>
  struct SIMT { 
    static_assert( std::is_same<T,double>::value ||
                   std::is_same<T,float>::value,
                   "KokkosKernels:: Invalid SIMT<> type." );

    static_assert( !std::is_same<SpT,Kokkos::Serial>::value &&
                   !std::is_same<SpT,Kokkos::OpenMP>::value,
                   "KokkosKernels:: Invalid AVX<> exec space." );

#if defined(KOKKOS_HAVE_CUDA)
    static_assert( std::is_same<SpT,Kokkos::Cuda>::value,
                   "KokkosKernels:: Invalid SIMT<> exec space." );
#endif

    using value_type = T; 
    using exec_space = SpT;
  };
  
  template<class T, int l>
  struct VectorTag {
    using value_type = typename T::value_type;
    using exec_space = typename T::exec_space;
    using member_type = typename Kokkos::Impl::TeamPolicyInternal<exec_space>::member_type;
    
    static_assert( std::is_same<T,SIMD<value_type,exec_space> >::value || // host compiler vectorization
                   std::is_same<T,AVX<value_type, exec_space> >::value || // host AVX vectorization
                   std::is_same<T,SIMT<value_type,exec_space> >::value,   // cuda thread vectorization
                   "KokkosKernels:: Invalid VectorUnitTag<> type." );
    
    using type = VectorTag;
    enum : int { length = l };
  };
  

  // Tags for BLAS
  struct Trans {
    struct Transpose {};
    struct ConjTranspose {};
    struct NoTranspose {};
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
    struct Unit    { enum : bool { use_unit_diag = true }; };
    struct NonUnit { enum :  bool{ use_unit_diag = false}; };
  };

  struct Algo {
    struct Gemm { 
      struct Triple {};
      struct Blocked {
	enum : int { mb = 4,
		     nb = 4 };
      };
    };

    struct Trsm {
      struct Unblocked {};
      struct Blocked {
	enum : int { mb = 4 };
      };
    };

    struct LU {
      struct Unblocked {};
      struct Blocked {
	enum : int { mb = 4 };
      };
    };
  };

  struct Util {
    template<typename ViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static void 
    set(const ViewType A, const ScalarType alpha) {
      typedef typename ViewType::value_type value_type;

      const int 
        iend = A.dimension(0),
        jend = A.dimension(1),
        kend = iend*jend;
      
      
      const int 
        as0 = A.stride_0(),
        as1 = A.stride_1();

      value_type 
        *__restrict__ pA = A.data();

      if ( (iend == as0 && as1 == 1) ||
           (jend == as1 && as0 == 1) )
        for (int k=0;k<kend;++k)
          pA[k] = alpha;
      else
        if (as0 > as1) 
          for (int i=0;i<iend;++i) 
            for (int j=0;j<jend;++j)
              pA[i*as0+j*as1] = alpha;
        else 
          for (int j=0;j<jend;++j) 
            for (int i=0;i<iend;++i)
              pA[i*as0+j*as1] = alpha;
    }

    template<typename ViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static void 
    scale(ViewType A, const ScalarType alpha) {
      typedef typename ViewType::value_type value_type;

      const int 
        iend = A.dimension(0),
        jend = A.dimension(1),
        kend = iend*jend;      
      const int 
        as0 = A.stride_0(),
        as1 = A.stride_1();

      value_type 
        *__restrict__ pA = A.data();
      if ( (iend == as0 && as1 == 1) ||
           (jend == as1 && as0 == 1) )
        for (int k=0;k<kend;++k)
          pA[k] *= alpha;
      else
        if (as0 > as1) 
          for (int i=0;i<iend;++i) 
            for (int j=0;j<jend;++j)
              pA[i*as0+j*as1] *= alpha;
        else 
          for (int j=0;j<jend;++j) 
            for (int i=0;i<iend;++i)
              pA[i*as0+j*as1] *= alpha;
    }

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

#endif
