#pragma once
#ifndef __ICHOL_LEFT_BLOCKED_HPP__
#define __ICHOL_LEFT_BLOCKED_HPP__

/// \file ichol_left_blocked.hpp
/// \brief Blocked incomplete Chloesky factorization.
/// \author Kyungjoo Kim (kyukim@sandia.gov)
///
/// Unlike the dense matrix algebra, this blocked version of sparse 
/// factorization does not lead to an efficient computation. 
/// This algorithm is only for the testing and debugging purpose.

#include "partition.hpp"

#include "scale.hpp"
#include "dot.hpp"

#include "gemv.hpp"
#include "gemm.hpp"

#include "trsv.hpp"
#include "trsm.hpp"

namespace Example { 

  using namespace std;
  
  // use Lower Triangular part only
  template<>
  template<typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int 
  IChol<Uplo::Lower,Algo::LeftBlocked>::invoke(const CrsMatViewType A) {
    typedef typename CrsMatViewType::value_type   value_type;
    typedef typename CrsMatViewType::ordinal_type ordinal_type;

    ordinal_type mb = blocksize;

    CrsMatViewType ATL, ATR,      A00, A01, A02,
      /**/         ABL, ABR,      A10, A11, A12,
      /**/                        A20, A21, A22;
    
    Part_2x2(A,  ATL, ATR,
             /**/ABL, ABR, 
             0, 0, Partition::TopLeft);
    
    while (ATL.NumRows() < A.NumRows()) {
      Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                      /*******/ /**/  A10, A11, A12,
                      ABL, ABR, /**/  A20, A21, A22,  
                      mb, mb, Partition::BottomRight);
      // -----------------------------------------------------
      CrsMatViewType AB0, AB1;
      
      Merge_2x1(A11,
                A21, AB1);
      
      Merge_2x1(A10,
                A20, AB0);

      Gemm<Trans::NoTranspose,Trans::Transpose>::invoke(-1.0, AB0, A10, 1.0, AB1);

      int r_val = IChol<Uplo::Lower,Algo::LeftUnblocked>::invoke(A11);
      if (r_val) 
        return r_val;

      Trsm<Side::Right,Uplo::Lower,Trans::Transpose>::invoke(Diag::NonUnit, 1.0, A11, A21);
      
      // -----------------------------------------------------
      Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                       A10, A11, A12, /**/ /******/
                       A20, A21, A22, /**/ ABL, ABR,
                       Partition::TopLeft);
    }

    return 0;
  }

}

#endif
