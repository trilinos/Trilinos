/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_SPARSE_MV_HPP_
#define KOKKOS_SPARSE_MV_HPP_

#ifdef KOKKOS_HAVE_CXX11
#include <type_traits>
#endif // KOKKOS_HAVE_CXX11

namespace KokkosSparse {

template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void
spmv(const char mode[],
     const AlphaType& alpha_in,
     const AMatrix& A,
     const XVector& x,
     const BetaType& beta_in,
     const YVector& y,
     const RANK_TWO) {

  typedef typename Impl::GetCoeffView<AlphaType, typename XVector::device_type>::view_type alpha_view_type;
  typedef typename Impl::GetCoeffView<BetaType,  typename XVector::device_type>::view_type beta_view_type;

  //alpha_view_type alpha = Impl::GetCoeffView<AlphaType, typename XVector::device_type>::get_view(alpha_in,x.dimension_1());
  //beta_view_type  beta =  Impl::GetCoeffView<AlphaType, typename XVector::device_type>::get_view(beta_in, x.dimension_1());

#ifdef KOKKOS_HAVE_CXX11
  // Make sure that both x and y have the same rank.
  static_assert (XVector::rank == YVector::rank, "KokkosBlas::spmv: Vector ranks do not match.");
  // Make sure that y is non-const.
  static_assert (Kokkos::Impl::is_same<typename YVector::value_type,
                                       typename YVector::non_const_value_type>::value,
                 "KokkosBlas::spmv: Output Vector must be non-const.");

#else
  // We prefer to use C++11 static_assert, because it doesn't give
  // "unused typedef" warnings, like the constructs below do.
  //
  // Make sure that both x and y have the same rank.
  typedef typename
    Kokkos::Impl::StaticAssert<XVector::rank == YVector::rank>::type Blas1_spmv_vector_ranks_do_not_match;
#endif // KOKKOS_HAVE_CXX11

  // Check compatibility of dimensions at run time.
  if((mode[0]==NoTranspose[0])||(mode[0]==Conjugate[0])) {
    if ((x.dimension_1 () != y.dimension_1 ()) ||
        (A.numCols() > x.dimension_0()) ||
        (A.numRows() > y.dimension_0()) ){
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match: "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.dimension_0 () << " x " << x.dimension_1()
         << ", y: " << y.dimension_0 () << " x " << y.dimension_1()
         ;

      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  } else {
    if ((x.dimension_1 () != y.dimension_1 ()) ||
        (A.numCols() > y.dimension_0()) ||
        (A.numRows() > x.dimension_0()) ){
      std::ostringstream os;
      os << "KokkosBlas::spmv: Dimensions do not match (transpose): "
         << ", A: " << A.numRows () << " x " << A.numCols()
         << ", x: " << x.dimension_0 () << " x " << x.dimension_1()
         << ", y: " << y.dimension_0 () << " x " << y.dimension_1()
         ;

      Kokkos::Impl::throw_runtime_exception (os.str ());
    }
  }

  typedef KokkosSparse::CrsMatrix<typename AMatrix::const_value_type,
    typename AMatrix::const_ordinal_type,
    typename AMatrix::device_type,
    typename AMatrix::memory_traits,
    typename AMatrix::const_size_type> AMatrix_Internal;
  AMatrix_Internal A_i = A;

  // Call single vector version if appropriate
  if( x.dimension_1() == 1) {
    typedef Kokkos::View<typename XVector::const_value_type*,
      typename Kokkos::Impl::if_c<Kokkos::Impl::is_same<typename YVector::array_layout,Kokkos::LayoutLeft>::value,
                                  Kokkos::LayoutLeft,Kokkos::LayoutStride>::type,
      typename XVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,
      typename XVector::specialize> XVector_SubInternal;
    typedef Kokkos::View<typename YVector::non_const_value_type*,
      typename Kokkos::Impl::if_c<Kokkos::Impl::is_same<typename YVector::array_layout,Kokkos::LayoutLeft>::value,
                                  Kokkos::LayoutLeft,Kokkos::LayoutStride>::type,
      typename YVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename YVector::specialize> YVector_SubInternal;

    XVector_SubInternal x_i = Kokkos::subview(x,Kokkos::ALL(),0);
    YVector_SubInternal y_i = Kokkos::subview(y,Kokkos::ALL(),0);

    alpha_view_type alpha = Impl::GetCoeffView<AlphaType, typename XVector::device_type>::get_view(alpha_in,x.dimension_1());
    beta_view_type  beta =  Impl::GetCoeffView<BetaType, typename XVector::device_type>::get_view(beta_in, x.dimension_1());

    typename alpha_view_type::non_const_type::HostMirror h_alpha =
        Kokkos::create_mirror_view(alpha);
    Kokkos::deep_copy(h_alpha,alpha);
    typename beta_view_type::non_const_type::HostMirror h_beta =
        Kokkos::create_mirror_view(beta);
    Kokkos::deep_copy(h_beta,beta);
    spmv(mode,h_alpha(0),A,x_i,h_beta(0),y_i);
    return;
  } else {

    typedef Kokkos::View<typename XVector::const_value_type**,
      typename XVector::array_layout,
      typename XVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>,
      typename XVector::specialize> XVector_Internal;
    typedef Kokkos::View<typename YVector::non_const_value_type**,
      typename YVector::array_layout,
      typename YVector::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename YVector::specialize> YVector_Internal;

    XVector_Internal x_i = x;
    YVector_Internal y_i = y;

    typedef Kokkos::View<typename alpha_view_type::const_value_type*,
      typename alpha_view_type::array_layout,
      typename alpha_view_type::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename alpha_view_type::specialize> alpha_view_type_Internal;
    typedef Kokkos::View<typename beta_view_type::const_value_type*,
      typename beta_view_type::array_layout,
      typename beta_view_type::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>,
      typename beta_view_type::specialize> beta_view_type_Internal;

    //alpha_view_type_Internal alpha_c = alpha;
    //beta_view_type_Internal  beta_c  = beta;
    return Impl::SPMV_MV<typename alpha_view_type_Internal::value_type*,
                         typename alpha_view_type_Internal::array_layout,
                         typename alpha_view_type_Internal::device_type,
                         typename alpha_view_type_Internal::memory_traits,
                         typename alpha_view_type_Internal::specialize,
                         typename AMatrix_Internal::value_type,
                         typename AMatrix_Internal::ordinal_type,
                         typename AMatrix_Internal::device_type,
                         typename AMatrix_Internal::memory_traits,
                         typename AMatrix_Internal::size_type,
                         typename XVector_Internal::value_type**,
                         typename XVector_Internal::array_layout,
                         typename XVector_Internal::device_type,
                         typename XVector_Internal::memory_traits,
                         typename XVector_Internal::specialize,
                         typename beta_view_type_Internal::value_type*,
                         typename beta_view_type_Internal::array_layout,
                         typename beta_view_type_Internal::device_type,
                         typename beta_view_type_Internal::memory_traits,
                         typename beta_view_type_Internal::specialize,
                         typename YVector_Internal::value_type**,
                         typename YVector_Internal::array_layout,
                         typename YVector_Internal::device_type,
                         typename YVector_Internal::memory_traits,
                         typename YVector_Internal::specialize>::spmv_mv(mode,alpha_in,A,x,beta_in,y);
  }
}

template <class AlphaType, class AMatrix, class XVector, class BetaType, class YVector>
void
spmv(const char mode[],
     const AlphaType& alpha,
     const AMatrix& A,
     const XVector& x,
     const BetaType& beta,
     const YVector& y) {
  typedef typename Kokkos::Impl::if_c<XVector::rank==2,RANK_TWO,RANK_ONE>::type RANK_SPECIALISE;
  spmv(mode,alpha,A,x,beta,y,RANK_SPECIALISE());
}

}


#endif

