#pragma once
#ifndef __TRSV_HPP__
#define __TRSV_HPP__

/// \file trsv.hpp
/// \brief Sparse triangular matrix solve on given sparse patterns and a single rhs.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;
  
  template<typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int
  trsv_l_n_t(const int diag,
             const CrsMatViewType A,
             const CrsMatViewType x) {
    typedef typename CrsMatViewType::ordinal_type ordinal_type;
    typedef typename CrsMatViewType::value_type   value_type;
    
    // case that A.lower.no_transpose.non_unit, b.transpose (row vector)

    int r_val = 0;
    value_type zero = 0.0, one = 1.0;

    CrsMatViewType ATL, ATR,      A00,  a01,     A02,
      /**/         ABL, ABR,      a10t, alpha11, a12t,
      /**/                        A20,  a21,     A22;

    CrsMatViewType xL, xR,        x0,   xi1,     x2;
    
    Part_2x2(A,   ATL, ATR,
             /**/ ABL, ABR,
             0, 0, Partition::TopLeft);

    Part_1x2(x,   xL, xR,
             0, Partition::Left);

    while (ATL.NumRows() < A.NumRows()) {
      Part_2x2_to_3x3(ATL, ATR, /**/  A00,  a01,     A02,
                      /*******/ /**/  a10t, alpha11, a12t,
                      ABL, ABR, /**/  A20,  a21,     A22,
                      1, 1, Partition::BottomRight);    
      
      Part_1x2_to_1x3(xL, xR,   /**/  x0, xi1, x2,
                      1, Partition::Right);
      
      // -----------------------------------------------------

      auto xi = xi1.extractRow(0);
      ordinal_type id_xi = xi.Index(0);
      if (id_xi >= 0) {
        // unit diag is default
        value_type alpha_val = one;

        // Diag::NonUnit is used, retrive alpha
        if (diag == Diag::NonUnit) {
          auto alpha = alpha11.extractRow(0);
          ordinal_type id_alpha = alpha.Index(0);
          alpha_val = (id_alpha < 0 ? zero : alpha.Value(id_alpha));
          
          // if encounter null diag, return the -(row + 1)
          if (abs(alpha_val) == 0.0) {
            r_val = -(ATL.NumRows() + 1);
            break;
          }
        }

        // update xi: xi = (xi/alpha);
        value_type &xi_val = xi.Value(id_xi);
        xi_val /= alpha_val;

        // update x2: x2 = x2 - a21*xi;
        auto xx = x2.extractRow(0);
        for (ordinal_type j=0;j<xx.NumNonZeros();++j) {
          auto aa = a21.extractRow(xx.Col(j));
          xx.Value(j) -= (aa.ValueAtColumn(0)*xi_val);
        }
      }
      
      // -----------------------------------------------------
      
      Merge_3x3_to_2x2(A00,  a01,     A02,  /**/ ATL, ATR,
                       a10t, alpha11, a12t, /**/ /******/
                       A20,  a21,     A22,  /**/ ABL, ABR,
                       Partition::TopLeft);         
      
      Merge_1x3_to_1x2(x0, xi1, x2, /**/  xL, xR,   
                       Partition::Left);
    }    

    return r_val;
  }

}

#endif
