#ifndef __TACHO_SCALE_2X2_BLOCK_INVERSE_DIAGONALS_INTERNAL_HPP
#define __TACHO_SCALE_2X2_BLOCK_INVERSE_DIAGONALS_INTERNAL_HPP


/// \file  Tacho_Scale2x2_BlockInverseDiagonals_Internal.hpp
/// \brief Inverse scale
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  /// row exchange
  template<>
  struct Scale2x2_BlockInverseDiagonals<Side::Left,Algo::Internal> {
    template<typename MemberType,
             typename ViewTypeD,
             typename ViewTypeA>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(MemberType &member,
           const ViewTypeD &D,
           const ViewTypeA &A) {
      typedef typename ViewTypeA::non_const_value_type value_type;        
      
      if (A.extent(0) == D.extent(0)) {
        if (A.span() > 0) {
          const ordinal_type m = A.extent(0), n = A.extent(1);
          const value_type zero(0), one(1);
#if defined(__CUDA_ARCH__)
          Kokkos::parallel_for
            (Kokkos::TeamThreadRange(member, n),
             [&](const ordinal_type &j) {
              Kokkos::parallel_for
                (Kokkos::ThreadVectorRange(member, m),
                 [&](const ordinal_type &i) {
                  const value_type prev_offdiag = i > 0 ? D(i-1,1) : zero;
                  const value_type offidag = D(i,1);
                  if (offidag == zero) {
                    /// 1x1 block
                    const value_type inv_diag = one/D(i,0);
                    A(i,j) *= inv_diag;
                  } else if (prev_offdiag == zero) {
                    /// 2x2 block
                    const value_type a = D(i,0), b = D(i,1), c = D(i+1,0), d = D(i+1,1);
                    const value_type det = a*d-b*c;
                    const value_type ia = d/det, ib = -b/det, ic = -c/det, id = a/det;
                    const value_type x0 = A(i,j), x1 = A(i+1,j);
                    A(i,  j) = ia*x0 + ib*x1;
                    A(i+1,j) = ic*x0 + id*x1;
                  } else {
                    /// this is taken care in other threads
                  }
                });
            });
#else
          for (ordinal_type i=0;i<m;) {
            const ordinal_type offidag = D(i,1);
            if (offidag == zero) {
              /// 1x1 block
              const value_type diag = D(i,0);              
              if (diag == zero) {
                printf("Error: Scale2x2_BlockInverseDiagonals<Side::Left,Algo::Internal> Encounters zero diagonal\n");
              }
              const value_type inv_diag = one/diag;
              for (ordinal_type j=0;j<n;++j) {
                A(i,j) *= inv_diag;
              }
              ++i;
            } else {
              /// 2x2 block
              const value_type a = D(i,0), b = D(i,1), c = D(i+1,0), d = D(i+1,1);
              const value_type det = a*d-b*c;
              if (det == zero) {
                printf("Error: Scale2x2_BlockInverseDiagonals<Side::Left,Algo::Internal> Encounters zero determinant\n");
              }
              const value_type ia = d/det, ib = -b/det, ic = -c/det, id = a/det;
              for (ordinal_type j=0;j<n;++j) {
                const value_type x0 = A(i,j), x1 = A(i+1,j);
                A(i,  j) = ia*x0 + ib*x1;
                A(i+1,j) = ic*x0 + id*x1;
              }   
              i+=2;           
            }
          }
#endif
        }
      } else {
        printf("Error: Scale2x2_BlockInverseDiagonals<Side::Left,Algo::Internal> A is not square\n");
      }
      return 0;
    }
  };


}
#endif
