#ifndef __KOKKOSBATCHED_VECTOR_SIMD_VIEW_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_VIEW_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch"

namespace KokkosBatched {
  namespace Experimental {
    
    template<int dim>
    struct PackDim {
      enum : int { value = dim };
    }; 
    
    // temporary solution until kokkos support SIMD layout or I do support it 
    template<typename ViewType, typename PackDim>
    struct SimdViewAccess {
    private:
      ViewType _a;
      
    public:
      typedef typename ViewType::reference_type  reference_simd_type ;
      typedef typename ViewType::pointer_type    pointer_simd_type ;
      typedef typename ViewType::value_type      value_simd_type;

      typedef typename value_simd_type::value_type value_type;
      typedef value_type& reference_type;
      typedef value_type* pointer_type; 
      
      enum : int { rank = ViewType::rank };
      enum : int { pack_dim = PackDim::value };
      enum : int { vector_length = value_simd_type::vector_length };
      
      SimdViewAccess() : _a() {} 
      SimdViewAccess(const ViewType &a) : _a(a) {} 

      SimdViewAccess & operator=(const ViewType &b) {
        _a = b;
        return *this;
      }
      SimdViewAccess & operator=(const SimdViewAccess &b) {
        if (this != &b) {
          _a = b._a;
        }
        return *this;
      }

      template< typename iType >
      KOKKOS_INLINE_FUNCTION constexpr
      typename std::enable_if< std::is_integral<iType>::value , size_t >::type
      extent( const iType & r ) const
      { return _a.extent(r)*(r == PackDim::value ? vector_length : 1); }
    
      template< typename iType >
      KOKKOS_INLINE_FUNCTION constexpr
      typename std::enable_if< std::is_integral<iType>::value , int >::type
      extent_int( const iType & r ) const
      { return static_cast<int>(_a.extent(r)*(r == PackDim::value ? vector_length : 1)); }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE

      template< typename iType >
      KOKKOS_INLINE_FUNCTION constexpr
      typename std::enable_if< std::is_integral<iType>::value , size_t >::type
      dimension( const iType & r ) const { return extent( r ); }

      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const { return _a.extent(0)*(0 == PackDim::value ? vector_length : 1); }
      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const { return _a.extent(1)*(1 == PackDim::value ? vector_length : 1); }
      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const { return _a.extent(2)*(2 == PackDim::value ? vector_length : 1); }
      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const { return _a.extent(3)*(3 == PackDim::value ? vector_length : 1); }
      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const { return _a.extent(4)*(4 == PackDim::value ? vector_length : 1); }
      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const { return _a.extent(5)*(5 == PackDim::value ? vector_length : 1); }
      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const { return _a.extent(6)*(6 == PackDim::value ? vector_length : 1); }
      KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const { return _a.extent(7)*(7 == PackDim::value ? vector_length : 1); }

#endif

      KOKKOS_INLINE_FUNCTION constexpr size_t size() const { 
        return (_a.size() * vector_length);
      }

      KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return _a.span()*vector_length; }
      KOKKOS_INLINE_FUNCTION constexpr bool   span_span_is_contiguous() const { return _a.span_span_is_contiguous(); }
      KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const { return _a.data(); }

      /// rank 0
      /// this does not make sense as this is flat view to simd view
      
      /// rank 1
      template< typename I0 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,Args...>::value && 1 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , Args ... args ) const {
        return _a(i0/vector_length)[i0%vector_length];
      }
      
      /// rank 2 
      template< typename I0 , typename I1 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,I1,Args...>::value && 2 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , const I1 & i1 , Args ... args ) const {
        switch (PackDim::value) {
        case 0: return _a(i0/vector_length,i1)[i0%vector_length];
        case 1: break;
        default:break;
        }
        return _a(i0,i1/vector_length)[i1%vector_length];        
      }

      /// rank 3
      template< typename I0 , typename I1 , typename I2 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,I1,I2,Args...>::value && 3 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , Args ... args ) const {
        switch (PackDim::value) {
        case 0: return _a(i0/vector_length,i1,i2)[i0%vector_length];
        case 1: return _a(i0,i1/vector_length,i2)[i1%vector_length];
        case 2: break;
        default:break;
        }
        return _a(i0,i1,i2/vector_length)[i2%vector_length];
      }
      
      /// rank 4
      template< typename I0 , typename I1 , typename I2 , typename I3 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,I1,I2,I3,Args...>::value && 4 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 , Args ... args ) const {
        switch (PackDim::value) {
        case 0: return _a(i0/vector_length,i1,i2,i3)[i0%vector_length];
        case 1: return _a(i0,i1/vector_length,i2,i3)[i1%vector_length];
        case 2: return _a(i0,i1,i2/vector_length,i3)[i2%vector_length];
        case 3: break;
        default:break;
        }
        return _a(i0,i1,i2,i3/vector_length)[i3%vector_length];
      }

      /// rank 5
      template< typename I0 , typename I1 , typename I2 , typename I3 , typename I4 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,Args...>::value && 5 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 , const I4 & i4 , Args ... args ) const {
        switch (PackDim::value) {
        case 0: return _a(i0/vector_length,i1,i2,i3,i4)[i0%vector_length];
        case 1: return _a(i0,i1/vector_length,i2,i3,i4)[i1%vector_length];
        case 2: return _a(i0,i1,i2/vector_length,i3,i4)[i2%vector_length];
        case 3: return _a(i0,i1,i2,i3/vector_length,i4)[i3%vector_length];
        case 4: break;
        default:break;
        }
        return _a(i0,i1,i2,i3,i4/vector_length)[i4%vector_length];
      }

      /// rank 6
      template< typename I0 , typename I1 , typename I2 , typename I3 , typename I4 , typename I5 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,Args...>::value && 6 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 , const I4 & i4 , const I5 & i5 , Args ... args ) const {
        switch (PackDim::value) {
        case 0: return _a(i0/vector_length,i1,i2,i3,i4,i5)[i0%vector_length];
        case 1: return _a(i0,i1/vector_length,i2,i3,i4,i5)[i1%vector_length];
        case 2: return _a(i0,i1,i2/vector_length,i3,i4,i5)[i2%vector_length];
        case 3: return _a(i0,i1,i2,i3/vector_length,i4,i5)[i3%vector_length];
        case 4: return _a(i0,i1,i2,i3,i4/vector_length,i5)[i4%vector_length];
        case 5: break;
        default:break;
        }
        return _a(i0,i1,i2,i3,i4,i5/vector_length)[i5%vector_length];
      }
    
      /// rank 7
      template< typename I0 , typename I1 , typename I2 , typename I3 , typename I4 , typename I5 , typename I6 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6,Args...>::value && 7 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 , const I4 & i4 , const I5 & i5 , const I6 & i6 , Args ... args ) const {
        switch (PackDim::value) {
        case 0: return _a(i0/vector_length,i1,i2,i3,i4,i5,i6)[i0%vector_length];
        case 1: return _a(i0,i1/vector_length,i2,i3,i4,i5,i6)[i1%vector_length];
        case 2: return _a(i0,i1,i2/vector_length,i3,i4,i5,i6)[i2%vector_length];
        case 3: return _a(i0,i1,i2,i3/vector_length,i4,i5,i6)[i3%vector_length];
        case 4: return _a(i0,i1,i2,i3,i4/vector_length,i5,i6)[i4%vector_length];
        case 5: return _a(i0,i1,i2,i3,i4,i5/vector_length,i6)[i5%vector_length];
        case 6: break;
        default:break;
        }
        return _a(i0,i1,i2,i3,i4,i5,i6/vector_length)[i6%vector_length];
      }

      /// rank 8
      template< typename I0 , typename I1 , typename I2 , typename I3 , typename I4 , typename I5 , typename I6 , typename I7 , class ... Args >
      KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6,I7,Args...>::value && 8 == ViewType::rank, reference_type >::type
      operator()( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 , const I4 & i4 , const I5 & i5 , const I6 & i6 , const I7 & i7 , Args ... args ) const {
        switch (PackDim::value) {
        case 0: return _a(i0/vector_length,i1,i2,i3,i4,i5,i6,i7)[i0%vector_length];
        case 1: return _a(i0,i1/vector_length,i2,i3,i4,i5,i6,i7)[i1%vector_length];
        case 2: return _a(i0,i1,i2/vector_length,i3,i4,i5,i6,i7)[i2%vector_length];
        case 3: return _a(i0,i1,i2,i3/vector_length,i4,i5,i6,i7)[i3%vector_length];
        case 4: return _a(i0,i1,i2,i3,i4/vector_length,i5,i6,i7)[i4%vector_length];
        case 5: return _a(i0,i1,i2,i3,i4,i5/vector_length,i6,i7)[i5%vector_length];
        case 6: return _a(i0,i1,i2,i3,i4,i5,i6/vector_length,i7)[i6%vector_length];
        case 7: break;
        default:break;
        }
        return _a(i0,i1,i2,i3,i4,i5,i6,i7/vector_length)[i7%vector_length];
      }
    };
  }
}

#pragma GCC diagnostic pop
  
#endif
