#pragma once
#ifndef __ICHOL_LEFT_UNBLOCKED_HPP__
#define __ICHOL_LEFT_UNBLOCKED_HPP__

namespace Example { 

  using namespace std;
  
  // use Lower Triangular part only
  template<typename CrsMatrixView>
  int factorizeLeftUnblocked(const CrsMatrixView &A) {
    CrsMatrixView 
      ATL, ATR,
      ABL, ABR;
    
    CrsMatrixView
      A00,  a01,     A02,
      a10t, alpha11, a12t,
      A20,  a21,     A22;
    
    Part_2x2(A,  ATL, ATR,
             /**/ABL, ABR, 0, Partition::TopLeft);
    
    while (ATL.NumRows() < A.NumRows()) {
      Part_2x2_to_3x3(ATL, ATR, /**/  A00,  a01,     A02,
                      /*******/ /**/  a10t, alpha11, a12t,
                      ABL, ABR, /**/  A20,  a21,     A22,  
                      1, 1, Partition::BottomRight);
      // -----------------------------------------------------
      CrsMatrixView AB0, ab1;
      
      Merge_2x1(alpha11,
                a12,     ab1);
      
      Merge_2x1(a10t,
                A20,     AB0);
      
      // sparse gemv
      for (ordinal_type k=0;k<ab1.NumRows();++k) {
        auto &alpha = ab1(k,0);
        alpha = alpha - SpDot(AB0.ExtractRow(k), a10t.ExtractRow(0));
      }
      
      for (ordinal_type k=0;k<ab1.NumRows();++k) {
        auto &alpha = ab1(k,0);
        alpha = sqrt(alpha);
        
        // sparse inverse scal
        for (ordinal_type m=k+1;m<ab1.NumRows();++m) {
          ab1(m,0) = ab1(m,0)/alpha;
        } 
      }
      
      // -----------------------------------------------------
      Partition_3x3_to_2x2(A00,  a01,     A02,  /**/ ATL, ATR,
                           a10t, alpha11, a12t, /**/ /******/
                           A20,  a21,     A22,  /**/ ABL, ABR  
                           1, 1, Partition::TopLeft);
    }

    return 0;
  }

}

#endif
