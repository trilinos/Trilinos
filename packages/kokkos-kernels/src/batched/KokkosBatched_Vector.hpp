#ifndef __KOKKOSBATCHED_VECTOR_HPP__
#define __KOKKOSBATCHED_VECTOR_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"


// forward declaration
namespace KokkosBatched {
  namespace Experimental {
    template<typename T, int l>
    class Vector;

    template<typename T, int l>
    struct is_vector<Vector<SIMD<T>,l> > : public std::true_type {};

    template<typename ValueType, typename MemorySpace>
    struct DefaultVectorLength {
      enum : int { value = 1 };
    };
    
    template<>
    struct DefaultVectorLength<float,Kokkos::HostSpace> {
#if   defined(__AVX512F__)
      enum : int { value = 16 };
#elif defined(__AVX__) || defined(__AVX2__)
      enum : int { value = 8 };
#else
      enum : int { value = 16 };      
#endif
    };
    template<>
    struct DefaultVectorLength<double,Kokkos::HostSpace> {
#if   defined(__AVX512F__)
      enum : int { value = 8 };
#elif defined(__AVX__) || defined(__AVX2__)
      enum : int { value = 4 };
#else
      enum : int { value = 8 };      
#endif
    };
    template<>    
    struct DefaultVectorLength<Kokkos::complex<float>,Kokkos::HostSpace> {
#if   defined(__AVX512F__)
      enum : int { value = 8 };
#elif defined(__AVX__) || defined(__AVX2__)
      enum : int { value = 4 };
#else
      enum : int { value = 8 };      
#endif
    };
    template<>
    struct DefaultVectorLength<Kokkos::complex<double>,Kokkos::HostSpace> {
#if   defined(__AVX512F__)
      enum : int { value = 4 };
#elif defined(__AVX__) || defined(__AVX2__)
      enum : int { value = 2 };
#else 
      enum : int { value = 4 };      
#endif
    };

#if defined(KOKKOS_ENABLE_CUDA)
    template<>
    struct DefaultVectorLength<float,Kokkos::CudaSpace> {
      enum : int { value = 16 };
    };
    template<>
    struct DefaultVectorLength<double,Kokkos::CudaSpace> {
      enum : int { value = 16 };
    };
    template<>
    struct DefaultVectorLength<Kokkos::complex<float>,Kokkos::CudaSpace> {
      enum : int { value = 16 };
    };
    template<>
    struct DefaultVectorLength<Kokkos::complex<double>,Kokkos::CudaSpace> {
      enum : int { value = 16 };
    };
#endif

    template<typename T>
    struct MagnitudeScalarType;

    template<> struct MagnitudeScalarType<float> { typedef float type; };
    template<> struct MagnitudeScalarType<double> { typedef double type; };
    template<> struct MagnitudeScalarType<Kokkos::complex<float> > { typedef float type; };
    template<> struct MagnitudeScalarType<Kokkos::complex<double> > { typedef double type; };

    template<int l> struct MagnitudeScalarType<Vector<SIMD<float>,l> > { typedef float type; };
    template<int l> struct MagnitudeScalarType<Vector<SIMD<double>,l> > { typedef double type; };
    template<int l> struct MagnitudeScalarType<Vector<SIMD<Kokkos::complex<float> >,l> > { typedef float type; };
    template<int l> struct MagnitudeScalarType<Vector<SIMD<Kokkos::complex<double> >,l> > { typedef double type; };

    template<typename T>
    struct ValueScalarType;

    template<> struct ValueScalarType<float> { typedef float type; };
    template<> struct ValueScalarType<double> { typedef double type; };
    template<> struct ValueScalarType<Kokkos::complex<float> > { typedef Kokkos::complex<float> type; };
    template<> struct ValueScalarType<Kokkos::complex<double> > { typedef Kokkos::complex<double> type; };

    template<int l> struct ValueScalarType<Vector<SIMD<float>,l> > { typedef float type; };
    template<int l> struct ValueScalarType<Vector<SIMD<double>,l> > { typedef double type; };
    template<int l> struct ValueScalarType<Vector<SIMD<Kokkos::complex<float> >,l> > { typedef Kokkos::complex<float> type; };
    template<int l> struct ValueScalarType<Vector<SIMD<Kokkos::complex<double> >,l> > { typedef Kokkos::complex<double> type; };

    template<typename T> typename Kokkos::Details::ArithTraits<T>::mag_type RealPart(const T &val) { return val; }
    template<typename T> typename Kokkos::Details::ArithTraits<Kokkos::complex<T> >::mag_type RealPart(const Kokkos::complex<T> &val) { return val.real(); }
    template<typename T, int l> typename Kokkos::Details::ArithTraits<Vector<SIMD<Kokkos::complex<T> >,l> >::mag_type
    RealPart(const Vector<SIMD<Kokkos::complex<T> >,l> &val) {
      typename Kokkos::Details::ArithTraits<Vector<SIMD<Kokkos::complex<T> >,l> >::mag_type r_val;
      for (int i=0;i<l;++i) { r_val[i] = val[i].real(); }
      return r_val;
    }
    template<typename T> typename Kokkos::Details::ArithTraits<T>::mag_type ImagPart(const T &val) { return 0; }
    template<typename T> typename Kokkos::Details::ArithTraits<Kokkos::complex<T> >::mag_type ImagPart(const Kokkos::complex<T> &val) { return val.imag(); }
    template<typename T, int l> typename Kokkos::Details::ArithTraits<Vector<SIMD<Kokkos::complex<T> >,l> >::mag_type
    ImagPart(const Vector<SIMD<Kokkos::complex<T> >,l> &val) {
      typename Kokkos::Details::ArithTraits<Vector<SIMD<Kokkos::complex<T> >,l> >::mag_type r_val;
      for (int i=0;i<l;++i) { r_val[i] = val[i].imag(); }
      return r_val;
    }

  }
}

#include "KokkosBatched_Vector_SIMD.hpp"

// arith traits overload for vector types
namespace Kokkos {
  namespace Details {

    using namespace KokkosBatched::Experimental;

    template<typename T, int l>
    class ArithTraits<Vector<SIMD<T>,l> > { 
    public:
      typedef typename ArithTraits<T>::val_type val_scalar_type;
      typedef typename ArithTraits<T>::mag_type mag_scalar_type;

      typedef Vector<SIMD<val_scalar_type>,l> val_type;
      typedef Vector<SIMD<mag_scalar_type>,l> mag_type;
      
      static const bool is_specialized = ArithTraits<T>::is_specialized;
      static const bool is_signed = ArithTraits<T>::is_signed;
      static const bool is_integer = ArithTraits<T>::is_integer;
      static const bool is_exact = ArithTraits<T>::is_exact;
      static const bool is_complex = ArithTraits<T>::is_complex;
    };

  }
}

#endif
