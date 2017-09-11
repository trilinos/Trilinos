#ifndef __KOKKOSBATCHED_VECTOR_SIMD_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {
  namespace Experimental {

    template<typename T, typename SpT, int l>
    class Vector<VectorTag<SIMD<T,SpT>,l> > {
    public:
      using tag_type = VectorTag<SIMD<T,SpT>,l>;

      using type = Vector<tag_type>;
      using value_type = typename tag_type::value_type;
      using member_type = typename tag_type::member_type;

      enum : int { vector_length = tag_type::length };

      typedef value_type data_type[vector_length];

      KOKKOS_INLINE_FUNCTION
      static const char* label() { return "SIMD"; }

    private:
      mutable data_type _data;

    public:
      KOKKOS_INLINE_FUNCTION Vector() {
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = 0;
          });
      }
      template<typename ArgValueType>
      KOKKOS_INLINE_FUNCTION Vector(const ArgValueType val) {
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = val;
          });
      }
      KOKKOS_INLINE_FUNCTION Vector(const type &b) {
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = b._data[i];
          });
      }

      KOKKOS_INLINE_FUNCTION
      type& loadAligned(value_type const *p) {
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = p[i];
          });
        return *this;
      }

      KOKKOS_INLINE_FUNCTION
      type& loadUnaligned(value_type const *p) {
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            _data[i] = p[i];
          });
        return *this;
      }

      // AVX has aligned version and unaligned version;
      // aligned load store are recommended if memory is aligned
      // in this version, there is no difference between aligned and unaligned

      KOKKOS_INLINE_FUNCTION
      void storeAligned(value_type *p) const {
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            p[i] = _data[i];
          });
      }

      KOKKOS_INLINE_FUNCTION
      void storeUnaligned(value_type *p) const {
        Kokkos::parallel_for
          (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,member_type>(vector_length),
           [&](const int &i) {
            p[i] = _data[i];
          });
      }

      KOKKOS_INLINE_FUNCTION
      value_type& operator[](const int i) const {
        return _data[i];
      }

    };

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator + (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(VectorTag<SIMD<T,SpT>,l>::length),
         [&](const int &i) {
          r_val[i] = a[i] + b[i];
        });
      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator + (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      return a + Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator + (const typename VectorTag<SIMD<T,SpT>,l>::value_type a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) + b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator += (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a + b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator += (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      a = a + b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator ++ (Vector<VectorTag<SIMD<T,SpT>,l> > & a, int) {
      Vector<VectorTag<SIMD<T,SpT>,l> > a0 = a;
      a = a + 1.0;
      return a0;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator ++ (Vector<VectorTag<SIMD<T,SpT>,l> > & a) {
      a = a + 1.0;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator - (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(VectorTag<SIMD<T,SpT>,l>::length),
         [&](const int &i) {
          r_val[i] = a[i] - b[i];
        });
      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator - (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      return a - Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator - (const typename VectorTag<SIMD<T,SpT>,l>::value_type a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) - b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator - (Vector<VectorTag<SIMD<T,SpT>,l> > const & a) {
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(VectorTag<SIMD<T,SpT>,l>::length),
         [&](const int &i) {
          a[i] = -a[i];
        });
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator -= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a - b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator -= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      a = a - b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator -- (Vector<VectorTag<SIMD<T,SpT>,l> > & a, int) {
      Vector<VectorTag<SIMD<T,SpT>,l> > a0 = a;
      a = a - 1.0;
      return a0;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator -- (Vector<VectorTag<SIMD<T,SpT>,l> > & a) {
      a = a - 1.0;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator * (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(VectorTag<SIMD<T,SpT>,l>::length),
         [&](const int &i) {
          r_val[i] = a[i] * b[i];
        });

      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator * (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      return a * Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator * (const typename VectorTag<SIMD<T,SpT>,l>::value_type a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) * b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator *= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a * b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator *= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      a = a * b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator / (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      Vector<VectorTag<SIMD<T,SpT>,l> > r_val;
      Kokkos::parallel_for
        (Kokkos::Impl::ThreadVectorRangeBoundariesStruct<int,typename VectorTag<SIMD<T,SpT>,l>::member_type>(VectorTag<SIMD<T,SpT>,l>::length),
         [&](const int &i) {
          r_val[i] = a[i] / b[i];
        });
      return r_val;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator / (Vector<VectorTag<SIMD<T,SpT>,l> > const & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      return a / Vector<VectorTag<SIMD<T,SpT>,l> >(b);
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > operator / (const typename VectorTag<SIMD<T,SpT>,l>::value_type a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      return Vector<VectorTag<SIMD<T,SpT>,l> >(a) / b;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator /= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, Vector<VectorTag<SIMD<T,SpT>,l> > const & b) {
      a = a / b;
      return a;
    }

    template<typename T, typename SpT, int l>
    KOKKOS_INLINE_FUNCTION
    static Vector<VectorTag<SIMD<T,SpT>,l> > & operator /= (Vector<VectorTag<SIMD<T,SpT>,l> > & a, const typename VectorTag<SIMD<T,SpT>,l>::value_type b) {
      a = a / b;
      return a;
    }

  }
}

#endif
