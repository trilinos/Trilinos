#ifndef __KOKKOSBATCHED_VECTOR_SIMD_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Complex.hpp"

namespace KokkosBatched {
  namespace Experimental {

    template<typename T, typename SpT, int l>
    class Vector<VectorTag<SIMD<T,SpT>,l> > {
    public:
      using tag_type = VectorTag<SIMD<T,SpT>,l>;

      using type = Vector<tag_type>;
      using value_type = typename tag_type::value_type;
      using mag_type = typename Kokkos::Details::ArithTraits<value_type>::mag_type;

      // TODO:: do I need this ? 
      //   - for cuda, this is necessary (if this is used on cuda) and we mostly use layout left directly
      //   - for host, the lambda does not allow further compiler optimization 
      using member_type = typename tag_type::member_type;

      static const int vector_length = tag_type::length;

      typedef value_type data_type[vector_length];

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "SIMD"; }

    private:
      mutable data_type _data;

    public:
      KOKKOS_INLINE_FUNCTION Vector() {
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = 0;
          });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = 0;
#endif
      }
      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const ArgValueType val) {
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = val;
          });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = val;
#endif
      }
      KOKKOS_INLINE_FUNCTION Vector(const type &b) {
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = b._data[i];
          });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = b._data[i];
#endif
      }

      KOKKOS_INLINE_FUNCTION
      type& loadAligned(value_type const *p) {
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = p[i];
          });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
        for (int i=0;i<vector_length;++i)
          _data[i] = p[i];
#endif
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadUnaligned(value_type const *p) {
        return loadAligned(p);
      }
      
      // AVX has aligned version and unaligned version;
      // aligned load store are recommended if memory is aligned
      // in this version, there is no difference between aligned and unaligned

      KOKKOS_INLINE_FUNCTION
      void storeAligned(value_type *p) const {
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            p[i] = _data[i];
          });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
        for (int i=0;i<vector_length;++i)
          p[i] = _data[i];
#endif
      }

      KOKKOS_INLINE_FUNCTION
      void storeUnaligned(value_type *p) const {
        storeAligned(p);
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int i) const {
        return _data[i];
      }

    };

    /// ---------------------------------------------------------------------------------------------------    

    /// simd, simd

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> >
    operator + (Vector<VectorTag<SIMD<T,SpT>,l> > const & a,  Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(l),
         [&](const int &i) {
          r_val[i] = a[i] + b[i];
        });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] + b[i];
#endif
      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > &
    operator += (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a + b;
      return a;
    }

    /// simd, real

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> >
    operator + (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const T b) {
      return a + Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }
    
    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator + (const T a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) + b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > & 
    operator += (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const T b) {
      a = a + b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> >
    operator ++ (Vector<VectorTag<SIMD<T,SpT>,l> > & a, int) {
      Vector<VectorTag<SIMD<T,SpT>,l> > a0 = a;
      a = a + typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a0;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > & 
    operator ++ (Vector<VectorTag<SIMD<T,SpT>,l> > & a) {
      a = a + typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a;
    }

    /// simd complex, real

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >
    operator + (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const T b) {
      return a + Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }
    
    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator + (const T a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) + b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & 
    operator += (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const T b) {
      a = a + b;
      return a;
    }

    /// simd complex, complex 

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >
    operator + (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const Kokkos::complex<T> b) {
      return a + Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }
    
    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator + (const Kokkos::complex<T> a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) + b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & 
    operator += (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const Kokkos::complex<T> b) {
      a = a + b;
      return a;
    }

    /// ---------------------------------------------------------------------------------------------------

    /// simd, simd

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator - (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(l),
         [&](const int &i) {
          r_val[i] = a[i] - b[i];
        });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] - b[i];
#endif
      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator - (Vector<VectorTag<SIMD<T,SpT>,l> > const & a) {
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(l),
         [&](const int &i) {
          a[i] = -a[i];
        });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
      for (int i=0;i<l;++i)
        a[i] = -a[i];
#endif
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > &
    operator -= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a - b;
      return a;
    }

    /// simd, real

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator - (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const T b) {
      return a - Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator - (const T a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) - b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > &
    operator -= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const T b) {
      a = a - b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator -- (Vector<VectorTag<SIMD<T,SpT>,l> > & a, int) {
      Vector<VectorTag<SIMD<T,SpT>,l> > a0 = a;
      a = a - typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a0;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > & 
    operator -- (Vector<VectorTag<SIMD<T,SpT>,l> > & a) {
      a = a - typename Kokkos::Details::ArithTraits<T>::mag_type(1);
      return a;
    }

    /// simd complex, real
    
    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator - (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const T b) {
      return a - Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator - (const T a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) - b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > &
    operator -= (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const T b) {
      a = a - b;
      return a;
    }

    /// simd complex, complex
    
    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator - (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const Kokkos::complex<T> b) {
      return a - Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator - (const Kokkos::complex<T> a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) - b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > &
    operator -= (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const Kokkos::complex<T> b) {
      a = a - b;
      return a;
    }

    /// ---------------------------------------------------------------------------------------------------    

    /// simd, simd

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator * (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(l),
         [&](const int &i) {
          r_val[i] = a[i] * b[i];
        });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] * b[i];
#endif
      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > &
    operator *= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a * b;
      return a;
    }


    /// simd, real

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator * (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const T b) {
      return a * Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator * (const T a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) * b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > &
    operator *= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const T b) {
      a = a * b;
      return a;
    }

    /// simd complex, real

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator * (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const T b) {
      return a * Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator * (const T a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) * b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > &
    operator *= (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const T b) {
      a = a * b;
      return a;
    }

    /// simd complex, complex

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator * (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const Kokkos::complex<T> b) {
      return a * Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator * (const Kokkos::complex<T> a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) * b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > &
    operator *= (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const Kokkos::complex<T> b) {
      a = a * b;
      return a;
    }

    /// ---------------------------------------------------------------------------------------------------    

    /// simd, simd

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator / (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
#if                                                     \
  defined (KOKKOS_HAVE_CUDA) &&                         \
  defined (KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(l),
         [&](const int &i) {
          r_val[i] = a[i] / b[i];
        });
#else
#if defined( KOKKOS_ENABLE_PRAGMA_IVDEP )
#pragma ivdep
#endif
#if defined( KOKKOS_ENABLE_PRAGMA_VECTOR )
#pragma vector always
#endif
#if !defined( KOKKOS_DEBUG ) & defined( KOKKOS_ENABLE_PRAGMA_SIMD )
#pragma simd
#endif
      for (int i=0;i<l;++i)
        r_val[i] = a[i] / b[i];
#endif
      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > &
    operator /= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a / b;
      return a;
    }

    /// simd, real

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator / (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const T b) {
      return a / Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<T,SpT>,l> > 
    operator / (const T a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) / b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<T,SpT>,l> > &
    operator /= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const T b) {
      a = a / b;
      return a;
    }

    /// simd complex, real

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator / (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const T b) {
      return a / Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator / (const T a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) / b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > &
    operator /= (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const T b) {
      a = a / b;
      return a;
    }

    /// simd complex, complex

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator / (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & a, const Kokkos::complex<T> b) {
      return a / Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static 
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > 
    operator / (const Kokkos::complex<T> a, Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> >(a) / b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static
    Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > &
    operator /= (Vector<VectorTag<SIMD<Kokkos::complex<T>,SpT>,l> > & a, const Kokkos::complex<T> b) {
      a = a / b;
      return a;
    }

  }
}

#endif
