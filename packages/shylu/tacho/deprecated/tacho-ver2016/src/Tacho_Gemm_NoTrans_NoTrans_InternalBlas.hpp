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
       AlgoGemm::InternalBlas,Variant::Two>
  ::invoke(PolicyType &policy,
           MemberType &member,
           const ScalarType alpha,
           DenseExecViewTypeA &A,
           DenseExecViewTypeB &B,
           const ScalarType beta,
           DenseExecViewTypeC &C) {
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
      
      // gemm (triple loop)
      {
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
           MemberType &member,
           const ScalarType alpha,
           DenseExecViewTypeA &A,
           DenseExecViewTypeB &B,
           const ScalarType beta,
           DenseExecViewTypeC &C) {
    typedef typename DenseExecViewTypeA::ordinal_type ordinal_type;
    typedef typename DenseExecViewTypeA::value_type   value_type;

    const ordinal_type m = C.NumRows();
    const ordinal_type n = C.NumCols();
    const ordinal_type k = B.NumRows();
    
    // for now simple implementation
    if (m == 0 || n == 0 || ((alpha == 0 || k == 0) && (beta == 1))) return 0;

    // C = beta C + alpha AB
    
    if (member.team_rank() == 0) {
      if (alpha == 0) {
        if (beta == 0) {
          for (ordinal_type j=0;j<n;++j)
            for (ordinal_type i=0;i<m;++i)
              C.Value(i, j) = 0.0;
        } else {
          for (ordinal_type j=0;j<n;++j)
            for (ordinal_type i=0;i<m;++i)
              C.Value(i, j) = beta*C.Value(i, j);
        }
      } else {
        // scale beta
        if      (beta == 0.0) 
          for (ordinal_type j=0;j<n;++j)
            for (ordinal_type i=0;i<m;++i)
              C.Value(i, j) = 0.0;
        else if (beta != 1.0) 
          for (ordinal_type j=0;j<n;++j)
            for (ordinal_type i=0;i<m;++i)
              C.Value(i, j) = beta*C.Value(i, j);
        
        // gemm blocked 
        {
          constexpr ordinal_type mc = 128, nr = 128, kc = 32, nnr = 16;
          {
            // block update
            const ordinal_type mm = m/mc, nn = n/nr, kk = k/kc;
            for (ordinal_type l=0;l<kk;++l)      
              for (ordinal_type i=0;i<mm;++i) 
                for (ordinal_type j=0;j<nn;++j) {
                  const ordinal_type loff = l*kc, moff = i*mc, noff = j*nr;
               
                  // GEBP : C_ij += A_il B_lj; 
                  {
                    constexpr ordinal_type np = (nr/nnr);
                    for (ordinal_type p=0;p<np;++p) {
                      const ordinal_type poff = p*nnr;
                      for (ordinal_type ll=0;ll<kc;++ll)      
                        for (ordinal_type ii=0;ii<mc;++ii) 
                          for (ordinal_type jj=0;jj<nnr;++jj) 
                            C.Value(ii+moff, jj+noff+poff) 
                              += A.Value(ii+moff, ll+loff)*B.Value(ll+loff, jj+noff+poff);
                    }
                  }
                }
          }
          {
            // remainder
            const ordinal_type lbegin = (k - k%kc), ibegin = (m - m%mc), jbegin = (n - n%nr);
            for (ordinal_type l=lbegin;l<k;++l)       
              for (ordinal_type i=ibegin;i<m;++i)
                for (ordinal_type j=jbegin;j<n;++j) 
                  C.Value(i, j) += A.Value(i, l)*B.Value(l, j);
          }
        }
      }        
    } 
    
    return 0;
  }


}

#endif
