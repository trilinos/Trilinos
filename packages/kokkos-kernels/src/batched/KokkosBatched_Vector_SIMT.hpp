#ifndef __KOKKOSBATCHED_VECTOR_SIMT_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMT_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {
  namespace Experimental {

    /// I am not sure if this is ever useful

    template<typename T, typename SpT, int l>
    class Vector<VectorTag<SIMT<T,SpT>,l> > {
    public:
      using tag_type = VectorTag<SIMT<T,SpT>,l>;

      using type = Vector<tag_type>;
      using value_type = typename tag_type::value_type;
      using member_type = typename tag_type::member_type;

      enum : int { vector_length = tag_type::length };

      typedef value_type data_type[vector_length];

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
      KOKKOS_INLINE_FUNCTION Vector(const value_type val) {
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

  }
}
#endif
