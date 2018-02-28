// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_UTILS_HPP
#define KOKKOS_VIEW_UTILS_HPP

namespace Kokkos {

namespace Impl {

KOKKOS_INLINE_FUNCTION
void raise_error(const char *msg)
{
#if defined(__CUDACC__) && defined(__CUDA_ARCH__)
  Kokkos::abort(msg);
#else
  throw std::runtime_error(msg);
#endif
}

template< class T , class Device > struct RebindStokhosStorageDevice ;

template< class T , class Device >
struct RebindStokhosStorageDevice< T * , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , class Device >
struct RebindStokhosStorageDevice< T [] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , unsigned N , class Device >
struct RebindStokhosStorageDevice< T[N] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type type[N] ;
};

// Get Sacado size from a list of dimensions
template <unsigned Rank> struct GetSacadoSize {};
template <> struct GetSacadoSize<0> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 = 0 ,
                      const size_t n2 = 0 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n0;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[0];
  }
};
template <> struct GetSacadoSize<1> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 = 0 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n1;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[1];
  }
};
template <> struct GetSacadoSize<2> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n2;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[2];
  }
};
template <> struct GetSacadoSize<3> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n3;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[3];
  }
};
template <> struct GetSacadoSize<4> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n4;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[4];
  }
};
template <> struct GetSacadoSize<5> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n5;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[5];
  }
};
template <> struct GetSacadoSize<6> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 ,
                      const size_t n7 = 0 ) {
    return n6;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[6];
  }
};
template <> struct GetSacadoSize<7> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 ,
                      const size_t n7 ) {
    return n7;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[7];
  }
};

} // namespace Impl

// Typename of flat array where sacado dimension is folded into neighbor
template <typename view_type, typename Enabled = void>
struct FlatArrayType {
  typedef view_type type;
};

// Typename of the intrinsic scalar type in a view
template <typename view_type, typename Enabled = void>
struct IntrinsicScalarType {
  typedef typename view_type::array_type::non_const_value_type type;
};

template <typename ViewType>
ViewType
make_view(const std::string& label,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(label, N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(const ViewAllocateWithoutInitializing& init,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(init, N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(typename ViewType::pointer_type ptr,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(ptr, N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(const std::string& label,
          const Impl::WithoutInitializing_t& init,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(view_alloc(label,init),
                  N0, N1, N2, N3, N4, N5, N6, N7);
}

namespace Impl {

// Specialization for deep_copy( view, view::value_type ) for Cuda

template <class OutputView, typename Enabled = void>
struct StokhosViewFill
{
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::execution_space execution_space ;

  const OutputView output ;
  const_value_type input ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_t i0 ) const
  {
    const size_t n1 = output.dimension_1();
    const size_t n2 = output.dimension_2();
    const size_t n3 = output.dimension_3();
    const size_t n4 = output.dimension_4();
    const size_t n5 = output.dimension_5();
    const size_t n6 = output.dimension_6();
    const size_t n7 = output.dimension_7();

    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_t i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_t i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_t i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_t i7 = 0 ; i7 < n7 ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input ;
    }}}}}}}
  }

  StokhosViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      const size_t n0 = output.dimension_0();
      Kokkos::RangePolicy<execution_space> policy( 0, n0 );
      Kokkos::parallel_for( policy, *this );
      execution_space::fence();
    }
};
} // namespace Impl

} // namespace Kokkos

#endif // KOKKOS_VIEW_UTILS_HPP
