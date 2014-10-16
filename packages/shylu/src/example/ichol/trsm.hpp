#pragma once
#ifndef __TRSM_HPP__
#define __TRSM_HPP__

namespace Example { 

  using namespace std;
  
  template<typename CrsMatViewType>
  inline int
  trsm_r_l_t(const int diag,
             const typename CrsMatViewType::value_type alpha,
             const CrsMatViewType A,
             const CrsMatViewType B) {
    
    CrsMatViewType BT,   B0,
      /**/         BB,   b1t,
      /**/               B2;

    scale(alpha, B);
    
    Part_2x1(B,   BT,
             /**/ BB,
             0, Partition::Top);

    while (BT.NumRows() < B.NumRows()) {
      Part_2x1_to_3x1(BT,  B0,
                      /**/ b1t,
                      BB,  B2,
                      1, Partition::Bottom);
      // -----------------------------------------------------

      trsv_l_n_t(diag, A, b1t);

      // -----------------------------------------------------      
      Merge_3x1_to_2x1(B0,  BT,
                       b1t, 
                       B2,  BB,
                       Partition::Top);
    }
    
    return 0;
  }

}

#endif
