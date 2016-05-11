#ifndef __TACHO_GEMM_NOTRANS_NOTRANS_INTERNAL_BLAS_HPP__
#define __TACHO_GEMM_NOTRANS_NOTRANS_INTERNAL_BLAS_HPP__

/// \file Tacho_Gemm_NoTrans_NoTrans_InternalBlas.hpp
/// \brief BLAS matrix-matrix multiplication 
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {
  /// BLAS Gemm
  /// =========
  /// Properties:
  /// - Compile with Device (o),
  /// - Callable in KokkosFunctors (o)
  /// - For now, this is for HostSpace only.
  template<>
  template<typename PolicyType,
           typename MemberType,
           typename ScalarType,
           typename DenseExecViewTypeA,
           typename DenseExecViewTypeB,
           typename DenseExecViewTypeC>
  KOKKOS_INLINE_FUNCTION
  int
  Gemm<Trans::NoTranspose,Trans::NoTranspose,
       AlgoGemm::InternalBlas,Variant::One>
  ::invoke(PolicyType &policy,
           const MemberType &member,
           const ScalarType alpha,
           DenseExecViewTypeA &A,
           DenseExecViewTypeB &B,
           const ScalarType beta,
           DenseExecViewTypeC &C) {
    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                Kokkos::Cuda
    //                >::value,
    //                "Cuda space is not available for calling external BLAS" );

    // static_assert( Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeA::space_type,
    //                typename DenseMatrixTypeB::space_type
    //                >::value && 
    //                Kokkos::Impl::is_same<
    //                typename DenseMatrixTypeB::space_type,
    //                typename DenseMatrixTypeC::space_type
    //                >::value,
    //                "Space type of input matrices does not match" );
    
    //typedef typename DenseExecViewTypeA::space_type   space_type;
    typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
    typedef typename DenseExecViewTypeA::value_type   value_type;

    if (member.team_rank() == 0) {
      const ordinal_type m = C.NumRows();
      const ordinal_type n = C.NumCols();
      const ordinal_type k = B.NumRows();

      // for now simple implementation
      if (m == 0 || n == 0 || ((alpha == 0 || k == 0) && (beta == 1))) return 0;
      
      if (alpha == 0) {
        if (beta == 0) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, n),
                               [&](const ordinal_type j) {
                                 for (ordinal_type i=0;i<m;++i)
                                   C.Value(i, j) = 0.0;
                               });
        } else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, n),
                               [&](const ordinal_type j) {
                                 for (ordinal_type i=0;i<m;++i)
                                   C.Value(i, j) = beta*C.Value(i, j);
                               });
        }
      } else {

        // scale beta
        if      (beta == 0.0) 
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, n),
                               [&](const ordinal_type j) {
                                 for (ordinal_type i=0;i<m;++i)
                                   C.Value(i, j) = 0.0;
                               });
        else if (beta != 1.0) 
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, n),
                               [&](const ordinal_type j) {
                                 for (ordinal_type i=0;i<m;++i)
                                   C.Value(i, j) = beta*C.Value(i, j);
                               });
        
        // gemm
        for (ordinal_type l=0;l<k;++l) {      
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, n),
                               [&](const ordinal_type j) {
                                 const value_type tmp = B.Value(l, j);
                                 //#pragma unroll
                                 for (ordinal_type i=0;i<m;++i)
                                   C.Value(i, j) += A.Value(i, l)*tmp;
                               });
          member.team_barrier();
        }
      } 
    }

    return 0;
  }
}

#endif
