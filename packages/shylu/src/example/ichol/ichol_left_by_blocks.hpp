#pragma once
#ifndef __ICHOL_LEFT_BY_BLOCKS_HPP__
#define __ICHOL_LEFT_BY_BLOCKS_HPP__

#include "partition.hpp"

#include "scale.hpp"
#include "dot.hpp"

#include "gemv.hpp"
#include "gemm.hpp"

#include "trsv.hpp"
#include "trsm.hpp"

#include "ichol_left_unblocked.hpp"

#include "gen_tasks_ichol_gemm.hpp"
#include "gen_tasks_ichol_scalar.hpp"
#include "gen_tasks_ichol_trsm.hpp"

namespace Example { 

  using namespace std;
  
  // use Lower Triangular part only
  template<typename CrsTaskViewType>
  inline int 
  ichol_left_by_blocks_lower(const CrsTaskViewType A) {
    // if succeed, return 0 
    int r_val = 0;

    CrsTaskViewType ATL, ATR,      A00, A01, A02,
      /**/          ABL, ABR,      A10, A11, A12,
      /**/                         A20, A21, A22;
    
    Part_2x2(A,  ATL, ATR,
             /**/ABL, ABR, 
             0, 0, Partition::TopLeft);
    
    while (ATL.NumRows() < A.NumRows()) {
      Part_2x2_to_3x3(ATL, ATR, /**/  A00, A01, A02,
                      /*******/ /**/  A10, A11, A12,
                      ABL, ABR, /**/  A20, A21, A22,  
                      1, 1, Partition::BottomRight);
      // -----------------------------------------------------
      CrsTaskViewType AB0, AB1;
      
      Merge_2x1(A11,
                A21, AB1);
      
      Merge_2x1(A10,
                A20, AB0);

      // sparse gemm
      gen_tasks_ichol_gemm(AB0, A10, AB1);

      // cholesky on diagonal block
      gen_tasks_ichol_scalar(A11);

      // trsm
      gen_tasks_ichol_trsm(A11, A21);

      // -----------------------------------------------------
      Merge_3x3_to_2x2(A00, A01, A02, /**/ ATL, ATR,
                       A10, A11, A12, /**/ /******/
                       A20, A21, A22, /**/ ABL, ABR,
                       Partition::TopLeft);
    }

    return r_val;
  }

}

#endif
