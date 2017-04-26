#ifndef __KOKKOSKERNELS_TRSM_SERIAL_IMPL_HPP__
#define __KOKKOSKERNELS_TRSM_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_Util.hpp"

#include "KokkosKernels_Gemm_Decl.hpp"
#include "KokkosKernels_Gemm_Serial_Impl.hpp"

namespace KokkosKernels {

  namespace Serial {



    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftLowerUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, const int n,
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m x n) = B (m x n)

      if (!(m && n)) return 0;

      for (int p=0;p<m;++p) {
        const ValueType
          *__restrict__ a21 = A + (p+1)*_as0 + p*_as1;
        ValueType 
          *__restrict__ b1t = B + p*_bs0,
          *__restrict__ B2  = b1t + _bs0;
        
        const int
          iend = m-p-1,
          jend = n;
        
        for (int i=0;i<iend;++i)
          for (int j=0;j<jend;++j) 
            B2[i*_bs0+j*_bs1] -= a21[i*_as0] * b1t[j*_bs1];
      }
      return 0;
    }

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftLowerUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           const int n,
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m x n) = B (m x n)

      if (!n) return 0;

      const ValueType 
        a_10 = A[1*_as0+0*_as1], 
        a_20 = A[2*_as0+0*_as1], a_21 = A[2*_as0+1*_as1], 
        a_30 = A[3*_as0+0*_as1], a_31 = A[3*_as0+1*_as1], a_32 = A[3*_as0+2*_as1];

      auto trsv = [&](const int p,
                      ValueType &b_0p, 
                      ValueType &b_1p, 
                      ValueType &b_2p, 
                      ValueType &b_3p) {
        // load
        b_0p = B[0*_bs0+p*_bs1];
        b_1p = B[1*_bs0+p*_bs1];
        b_2p = B[2*_bs0+p*_bs1];
        b_3p = B[3*_bs0+p*_bs1];
        
        // 0 iteration
        b_1p -= a_10 * b_0p; 
        b_2p -= a_20 * b_0p; 
        b_3p -= a_30 * b_0p; 

        // 1 iteration
        b_2p -= a_21 * b_1p;
        b_3p -= a_31 * b_1p;

        // 2 iteration
        b_3p -= a_32 * b_2p; 

        // store
        B[1*_bs0+p*_bs1] = b_1p;
        B[2*_bs0+p*_bs1] = b_2p;
        B[3*_bs0+p*_bs1] = b_3p;
      };

      ValueType b_p[4];
      for (int p=0;p<n;++p) 
        trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);

      // unroll right hand side twice
      // ValueType b_p[2][4];
      // const int nn = (n/2)*2;
      // for (int p=0;p<nn;p+=2) {
      //   trsv(p,   b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);
      //   trsv(p+1, b_p[1][0], b_p[1][1], b_p[1][2], b_p[1][3]);
      // }

      // if (n%2) 
      //   trsv(nn,  b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);

      return 0;
    }

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftLowerNonUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, const int n,
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m x n) = B (m x n)

      if (!(m && n)) return 0;

      for (int p=0;p<m;++p) {
        const ValueType
          inv_alpha11 = 1.0/A[p*_as0+p*_as1], //alpha11 = A[p*_as0 + p*_as1],
          *__restrict__ a21 = A + (p+1)*_as0 + p*_as1;
        ValueType 
          *__restrict__ b1t = B + p*_bs0,
          *__restrict__ B2  = b1t + _bs0;
        
        const int
          iend = m-p-1,
          jend = n;

        // inverse scale
        for (int j=0;j<jend;++j)
          b1t[j*_bs1] *= inv_alpha11;

        // division
        // for (int j=0;j<jend;++j)
        //   b1t[j*_bs1] /= alpha11;

        for (int i=0;i<iend;++i)
          for (int j=0;j<jend;++j) 
            B2[i*_bs0+j*_bs1] -= a21[i*_as0] * b1t[j*_bs1];
      }
      return 0;
    }

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftLowerNonUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           const int n,
           /**/  ValueType *__restrict__ B) {
      // L(m x m) X(m x n) = B (m x n)

      if (!n) return 0;

      const ValueType 
        a_10 = A[1*_as0+0*_as1], 
        a_20 = A[2*_as0+0*_as1], a_21 = A[2*_as0+1*_as1], 
        a_30 = A[3*_as0+0*_as1], a_31 = A[3*_as0+1*_as1], a_32 = A[3*_as0+2*_as1];

      // const ValueType 
      //   a_00 = A[0*_as0+0*_as1], 
      //   a_11 = A[1*_as0+1*_as1], 
      //   a_22 = A[2*_as0+2*_as1], 
      //   a_33 = A[3*_as0+3*_as1];

      const ValueType           
        inv_a_00 = 1.0/A[0*_as0+0*_as1],
        inv_a_11 = 1.0/A[1*_as0+1*_as1],
        inv_a_22 = 1.0/A[2*_as0+2*_as1],
        inv_a_33 = 1.0/A[3*_as0+3*_as1];
      
      auto trsv = [&](const int p,
                      ValueType &b_0p, 
                      ValueType &b_1p, 
                      ValueType &b_2p, 
                      ValueType &b_3p) {
        // load
        b_0p = B[0*_bs0+p*_bs1]; 
        b_1p = B[1*_bs0+p*_bs1]; 
        b_2p = B[2*_bs0+p*_bs1]; 
        b_3p = B[3*_bs0+p*_bs1]; 

        // 0 iteration
        b_0p *= inv_a_00; /* b_0p /= a_00;*/   
        b_1p -= a_10 * b_0p;                  
        b_2p -= a_20 * b_0p;                  
        b_3p -= a_30 * b_0p;                    

        // 1 iteration                         
        b_1p *= inv_a_11; /* b_1p /= a_11; */  
        b_2p -= a_21 * b_1p;                   
        b_3p -= a_31 * b_1p;                   
                                                
        // 2 iteration                         
        b_2p *= inv_a_22; /* b_2p /= a_22; */  
        b_3p -= a_32 * b_2p;                   
                                                
        // 3 iteration                         
        b_3p *= inv_a_33; /* b_3p /= a_33; */     

        // store
        B[0*_bs0+p*_bs1] = b_0p; 
        B[1*_bs0+p*_bs1] = b_1p; 
        B[2*_bs0+p*_bs1] = b_2p; 
        B[3*_bs0+p*_bs1] = b_3p; 
      };

      ValueType b_p[4];
      for (int p=0;p<n;++p) 
        trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);

      // ValueType b_p[2][4];
      // const int nn = (n/2)*2;
      // for (int p=0;p<nn;p+=2) {
      //   trsv(p,   b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);
      //   trsv(p+1, b_p[1][0], b_p[1][1], b_p[1][2], b_p[1][3]);
      // }

      // if (n%2) 
      //   trsv(nn,  b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);

      return 0;
    }

    
    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmRightUpperUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, const int n,
           /**/  ValueType *__restrict__ B) {
      // X(m x n) U(n x n) = B (m x n)

      if (!(m && n)) return 0;

      for (int p=0;p<n;++p) {
        const ValueType
          alpha11 = A[p*_as0 + p*_as1],
          *__restrict__ a12t = A + p*_as0 + (p+1)*_as1;
        ValueType
          *__restrict__ b1 = B + p*_bs1,
          *__restrict__ B2 = b1 + _bs1;
        
        const int
          iend = m,
          jend = n-p-1;
        
        for (int i=0;i<iend;++i) 
          for (int j=0;j<jend;++j)
            B2[i*_bs0+j*_bs1] -= a12t[j*_as1] * b1[i*_bs0];
      }
      return 0;
    }


    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmRightUpperUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           const int m,
           /**/  ValueType *__restrict__ B) {
      // X(m x n) U(n x n) = B (m x n)

      if (!m) return 0;

      const ValueType 
        a_01 = A[0*_as0+1*_as1], a_02 = A[0*_as0+2*_as1], a_03 = A[0*_as0+3*_as1], 
        /**/                     a_12 = A[1*_as0+2*_as1], a_13 = A[1*_as0+3*_as1], 
        /**/                                              a_23 = A[2*_as0+3*_as1];      
      
      auto trsv = [&](const int p,
                      ValueType &b_p0, 
                      ValueType &b_p1, 
                      ValueType &b_p2, 
                      ValueType &b_p3) {
        // load
        b_p0 = B[p*_bs0+0*_bs1];
        b_p1 = B[p*_bs0+1*_bs1];
        b_p2 = B[p*_bs0+2*_bs1];
        b_p3 = B[p*_bs0+3*_bs1];

        // 0 iteration
        b_p1 -= b_p0 * a_01; 
        b_p2 -= b_p0 * a_02; 
        b_p3 -= b_p0 * a_03;
 
        // 1 iteration
        b_p2 -= b_p1 * a_12;
        b_p3 -= b_p1 * a_13;

        // 2 iteration
        b_p3 -= b_p2 * a_23; 

        // store
        B[p*_bs0+1*_bs1] = b_p1;
        B[p*_bs0+2*_bs1] = b_p2;
        B[p*_bs0+3*_bs1] = b_p3;
      };

      ValueType b_p[4];
      for (int p=0;p<m;++p) 
        trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);
      
      // ValueType b_p[2][4];
      // const int mm = (m/2)*2;
      // for (int p=0;p<mm;p+=2) {
      //   trsv(p,   b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);
      //   trsv(p+1, b_p[1][0], b_p[1][1], b_p[1][2], b_p[1][3]);
      // }

      // if (m%2) 
      //   trsv(mm,  b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);

      return 0;
    }

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmRightUpperNonUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, const int n,
           /**/  ValueType *__restrict__ B) {
      // X(m x n) U(n x n) = B (m x n)

      if (!(m && n)) return 0;

      for (int p=0;p<n;++p) {
        const ValueType
          inv_alpha11 = 1.0/A[p*_as0+p*_as1], // alpha11 = A[p*_as0 + p*_as1],
          *__restrict__ a12t = A + p*_as0 + (p+1)*_as1;
        ValueType
          *__restrict__ b1 = B + p*_bs1,
          *__restrict__ B2 = b1 + _bs1;
        
        const int
          iend = m,
          jend = n-p-1;

        for (int i=0;i<iend;++i)
          b1[i*_bs0] *= inv_alpha11;

        // for (int i=0;i<iend;++i)
        //   b1[i*_bs0] /= alpha11;
        
        for (int i=0;i<iend;++i) 
          for (int j=0;j<jend;++j)
            B2[i*_bs0+j*_bs1] -= a12t[j*_as1] * b1[i*_bs0];
      }
      return 0;
    }


    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmRightUpperNonUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           const int m,
           /**/  ValueType *__restrict__ B) {
      // X(m x n) U(n x n) = B (m x n)

      if (!m) return 0;

      const ValueType 
        a_01 = A[0*_as0+1*_as1], a_02 = A[0*_as0+2*_as1], a_03 = A[0*_as0+3*_as1], 
        /**/                     a_12 = A[1*_as0+2*_as1], a_13 = A[1*_as0+3*_as1], 
        /**/                                              a_23 = A[2*_as0+3*_as1];

      // const ValueType 
      //   a_00 = A[0*_as0+0*_as1], 
      //   a_11 = A[1*_as0+1*_as1], 
      //   a_22 = A[2*_as0+2*_as1], 
      //   a_33 = A[3*_as0+3*_as1];

      const ValueType 
        inv_a_00 = 1.0/A[0*_as0+0*_as1], 
        inv_a_11 = 1.0/A[1*_as0+1*_as1], 
        inv_a_22 = 1.0/A[2*_as0+2*_as1], 
        inv_a_33 = 1.0/A[3*_as0+3*_as1];      
      
      auto trsv = [&](const int p,
                      ValueType &b_p0, 
                      ValueType &b_p1, 
                      ValueType &b_p2, 
                      ValueType &b_p3) {
        // load
        b_p0 = B[p*_bs0+0*_bs1];
        b_p1 = B[p*_bs0+1*_bs1];
        b_p2 = B[p*_bs0+2*_bs1];
        b_p3 = B[p*_bs0+3*_bs1];

        // 0 iteration
        b_p0 *= inv_a_00;  // b_p0 /= a_00; 
        b_p1 -= b_p0 * a_01; 
        b_p2 -= b_p0 * a_02; 
        b_p3 -= b_p0 * a_03;
 
        // 1 iteration
        b_p1 *= inv_a_11; // b_p1 /= a_11;
        b_p2 -= b_p1 * a_12;
        b_p3 -= b_p1 * a_13;

        // 2 iteration
        b_p2 *= inv_a_22; // b_p2 /= a_22;
        b_p3 -= b_p2 * a_23; 

        // 3 iteration
        b_p3 *= inv_a_33; // b_p3 /= a_33;

        // store
        B[p*_bs0+0*_bs1] = b_p0;
        B[p*_bs0+1*_bs1] = b_p1;
        B[p*_bs0+2*_bs1] = b_p2;
        B[p*_bs0+3*_bs1] = b_p3;
      };

      ValueType b_p[4];
      for (int p=0;p<m;++p) 
        trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);

      // ValueType b_p[2][4];
      // const int mm = (m/2)*2;
      // for (int p=0;p<mm;p+=2) {
      //   trsv(p,   b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);
      //   trsv(p+1, b_p[1][0], b_p[1][1], b_p[1][2], b_p[1][3]);
      // }

      // if (m%2) 
      //   trsv(mm,  b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);

      return 0;
    }


    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftUpperUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, const int n,
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m x n) = B (m x n)

      if (!(m && n)) return 0;

      ValueType 
        *__restrict__ B0  = B;

      for (int p=(m-1);p>=0;--p) {
        const int
          iend = p,
          jend = n;

        const ValueType
          *__restrict__ a01 = A + p*_as1;

        ValueType 
          *__restrict__ b1t = B + p*_bs0;
        
        for (int i=0;i<iend;++i)
          for (int j=0;j<jend;++j) 
            B0[i*_bs0+j*_bs1] -= a01[i*_as0] * b1t[j*_bs1];
      }
      return 0;
    }

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftUpperUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           const int n,
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m x n) = B (m x n)

      if (!n) return 0;

      const ValueType 
        a_01 = A[0*_as0+1*_as1], a_02 = A[0*_as0+2*_as1], a_03 = A[0*_as0+3*_as1], 
        /**/                     a_12 = A[1*_as0+2*_as1], a_13 = A[1*_as0+3*_as1], 
        /**/                                              a_23 = A[2*_as0+3*_as1];

      auto trsv = [&](const int p,
                      ValueType &b_0p, 
                      ValueType &b_1p, 
                      ValueType &b_2p, 
                      ValueType &b_3p) {
        // load
        b_0p = B[0*_bs0+p*_bs1];
        b_1p = B[1*_bs0+p*_bs1];
        b_2p = B[2*_bs0+p*_bs1];
        b_3p = B[3*_bs0+p*_bs1];
        
        // 0 iteration
        b_0p -= a_03 * b_3p;
        b_1p -= a_13 * b_3p;
        b_2p -= a_23 * b_3p;

        // 1 iteration
        b_0p -= a_02 * b_2p; 
        b_1p -= a_12 * b_2p;

        // 2 iteration
        b_0p -= a_01 * b_1p;

        // store
        B[0*_bs0+p*_bs1] = b_0p;
        B[1*_bs0+p*_bs1] = b_1p;
        B[2*_bs0+p*_bs1] = b_2p;
      };
      
      ValueType b_p[4];
      for (int p=0;p<n;++p) 
        trsv(p, b_p[3], b_p[2], b_p[1], b_p[3]);

      // unroll right hand side twice
      // ValueType b_p[2][4];
      // const int nn = (n/2)*2;
      // for (int p=0;p<nn;p+=2) {
      //   trsv(p,   b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);
      //   trsv(p+1, b_p[1][0], b_p[1][1], b_p[1][2], b_p[1][3]);
      // }

      // if (n%2) 
      //   trsv(nn,  b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);

      return 0;
    }

    template<int bmn>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftUpperNonUnitDiag<bmn>::
    invoke(const ValueType *__restrict__ A,
           const int m, const int n,
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m x n) = B (m x n)

      if (!(m && n)) return 0;

      ValueType 
        *__restrict__ B0  = B;

      for (int p=(m-1);p>=0;--p) {
        const int
          iend = p,
          jend = n;

        const ValueType
          inv_alpha11 = 1.0/A[p*_as0+p*_as1], //alpha11 = A[p*_as0 + p*_as1],
          *__restrict__ a01 = A + p*_as1;

        ValueType 
          *__restrict__ b1t = B + p*_bs0;
        
        // inverse scale
        for (int j=0;j<jend;++j)
          b1t[j*_bs1] *= inv_alpha11;

        // division
        // for (int j=0;j<jend;++j)
        //   b1t[j*_bs1] /= alpha11;
        
        for (int i=0;i<iend;++i)
          for (int j=0;j<jend;++j) 
            B0[i*_bs0+j*_bs1] -= a01[i*_as0] * b1t[j*_bs1];
      }
      return 0;
    }

    template<>
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int
    InnerTrsmLeftUpperNonUnitDiag<4>::
    invoke(const ValueType *__restrict__ A,
           const int n,
           /**/  ValueType *__restrict__ B) {
      // U(m x m) X(m x n) = B (m x n)

      if (!n) return 0;

      const ValueType 
        a_01 = A[0*_as0+1*_as1], a_02 = A[0*_as0+2*_as1], a_03 = A[0*_as0+3*_as1], 
        /**/                     a_12 = A[1*_as0+2*_as1], a_13 = A[1*_as0+3*_as1], 
        /**/                                              a_23 = A[2*_as0+3*_as1];


      // const ValueType 
      //   a_00 = A[0*_as0+0*_as1], 
      //   a_11 = A[1*_as0+1*_as1], 
      //   a_22 = A[2*_as0+2*_as1], 
      //   a_33 = A[3*_as0+3*_as1];

      const ValueType           
        inv_a_00 = 1.0/A[0*_as0+0*_as1],
        inv_a_11 = 1.0/A[1*_as0+1*_as1],
        inv_a_22 = 1.0/A[2*_as0+2*_as1],
        inv_a_33 = 1.0/A[3*_as0+3*_as1];

      auto trsv = [&](const int p,
                      ValueType &b_0p, 
                      ValueType &b_1p, 
                      ValueType &b_2p, 
                      ValueType &b_3p) {
        // load
        b_0p = B[0*_bs0+p*_bs1];
        b_1p = B[1*_bs0+p*_bs1];
        b_2p = B[2*_bs0+p*_bs1];
        b_3p = B[3*_bs0+p*_bs1];
        
        // 0 iteration
        b_3p *= inv_a_33; /* b_3p /= a_33;*/   
        b_2p -= a_23 * b_3p;
        b_1p -= a_13 * b_3p;
        b_0p -= a_03 * b_3p;

        // 1 iteration
        b_2p *= inv_a_22; /* b_2p /= a_22; */  
        b_1p -= a_12 * b_2p;
        b_0p -= a_02 * b_2p; 

        // 2 iteration
        b_1p *= inv_a_11; /* b_1p /= a_11; */  
        b_0p -= a_01 * b_1p;

        // 3 iteration                         
        b_0p *= inv_a_00; /* b_0p /= a_00; */     

        // store
        B[0*_bs0+p*_bs1] = b_0p;
        B[1*_bs0+p*_bs1] = b_1p;
        B[2*_bs0+p*_bs1] = b_2p;
        B[3*_bs0+p*_bs1] = b_3p; 
      };
      
      ValueType b_p[4];
      for (int p=0;p<n;++p) 
        trsv(p, b_p[3], b_p[2], b_p[1], b_p[0]);

      // unroll right hand side twice
      // ValueType b_p[2][4];
      // const int nn = (n/2)*2;
      // for (int p=0;p<nn;p+=2) {
      //   trsv(p,   b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);
      //   trsv(p+1, b_p[1][0], b_p[1][1], b_p[1][2], b_p[1][3]);
      // }

      // if (n%2) 
      //   trsv(nn,  b_p[0][0], b_p[0][1], b_p[0][2], b_p[0][3]);

      return 0;
    }

    template<typename ArgDiag>
    struct Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsm::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        typedef typename BViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = B.dimension(0),
          n = B.dimension(1),
          vl = vector_type::vector_length;

        // no error check
        cblas_dtrsm_compact(CblasRowMajor, 
                            CblasLeft, CblasLower, CblasNoTrans, 
                            ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                            m, n, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (double*)B.data(), B.stride_0(), 
                            (MKL_INT)vl, (MKL_INT)1);
      }
    };
    
    template<typename ArgDiag>
    struct Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsm::Unblocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {

        typedef typename BViewType::value_type value_type;
        
        if (alpha == 0) {
          Util::set(B, 0);
        } else {
          if (alpha != 1)
            Util::scale(B, alpha);
          
          // B (m x n), A(m x m)
          const int
            m = B.dimension(0),
            n = B.dimension(1);

          if (!(m && n)) return 0;

          const int
            as0 = A.stride_0(),
            //as1 = A.stride_1(),
            bs0 = B.stride_0(),
            bs1 = B.stride_1();
          
          for (int p=0;p<m;++p) {
            const int
              iend = m-p-1,
              jend = n;

            const value_type
              *__restrict__ a21 = iend ? &A(p+1, p) : NULL;
            
            value_type
              *__restrict__ b1t =        &B(p,   0),
              *__restrict__ B2  = iend ? &B(p+1, 0) : NULL;
                        
            if (!ArgDiag::use_unit_diag) {
              const value_type alpha11 = A(p, p);
              for (int j=0;j<jend;++j)
                b1t[j*bs1] /= alpha11;
            }
            
            for (int i=0;i<iend;++i)
              for (int j=0;j<jend;++j)
                B2[i*bs0+j*bs1] -= a21[i*as0] * b1t[j*bs1];
          }
        }
        
        return 0;
      }
    };


    template<typename ArgDiag>
    struct Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,ArgDiag,Algo::Trsm::Blocked> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {

        typedef typename BViewType::value_type value_type;
        
        if (alpha == 0) {
          Util::set(B, 0);
        } else {
          if (alpha != 1)
            Util::scale(B, alpha);
          
          // Lower(A) (m x m), B(m x n)
          const int
            m = B.dimension(0),
            n = B.dimension(1);

          if (!(m && n)) return 0;

          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            bs0 = B.stride_0(),
            bs1 = B.stride_1();

	  enum : int {
            mb = Algo::Trsm::Blocked::mb };
          
          InnerTrsmLeftLowerUnitDiag<mb> trsm_llu(as0, as1,
                                                      bs0, bs1);
          
          InnerTrsmLeftLowerNonUnitDiag<mb> trsm_lln(as0, as1,
                                                         bs0, bs1);

          InnerRankUpdate<mb,mb> gemm_nn(as0, as1,
                                         bs0, bs1,
                                         bs0, bs1);
          
          const int 
            mm = (m/mb)*mb,
            nn = (n/mb)*mb;

          for (int p=0;p<mm;p+=mb) {
            // trsm update
            if (ArgDiag::use_unit_diag) 
              trsm_llu.invoke(&A(p,p), nn, &B(p,0));
            else 
              trsm_lln.invoke(&A(p,p), nn, &B(p,0));

            // gemm update
            for (int i=p+mb;i<mm;i+=mb)
              for (int j=0;j<nn;j+=mb)
                gemm_nn.serial_invoke(-1, &A(i,p), &B(p,j), mb, &B(i,j));
          }

          // remainder
          const int 
            mp = (m%mb),
            np = (n%mb);

          if (mp) {
            gemm_nn.serial_invoke(-1, &A(mm,0), &B(0,0), mp, nn, mm, &B(mm,0)); 
            if (ArgDiag::use_unit_diag) 
              trsm_llu.invoke(&A(mm,mm), mp, nn, &B(mm,0));
            else 
              trsm_lln.invoke(&A(mm,mm), mp, nn, &B(mm,0));
            
          }
          if (np) {
            if (ArgDiag::use_unit_diag) 
              trsm_llu.invoke(&A(0,0), m, np, &B(0,nn));
            else 
              trsm_lln.invoke(&A(0,0), m, np, &B(0,nn));
          }
        }
        return 0;
      }
    };

    template<typename ArgDiag>
    struct Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        typedef typename BViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = B.dimension(0),
          n = B.dimension(1),
          vl = vector_type::vector_length;

        // no error check
        cblas_dtrsm_compact(CblasRowMajor, 
                            CblasRight, CblasUpper, CblasNoTrans, 
                            ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                            m, n, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (double*)B.data(), B.stride_0(), 
                            (MKL_INT)vl, (MKL_INT)1);
      }
    };

    template<typename ArgDiag>
    struct Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Unblocked> {
      
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {

        typedef typename BViewType::value_type value_type;        

        if (alpha == 0) {
          Util::set(B, 0);
        } else {
          if (alpha != 1)
            Util::scale(B, alpha);
          
          // B (m x n), A(m x m)
          const int
            m = B.dimension(0),
            n = B.dimension(1);

          if (!(m && n)) return 0;

          const int
            //as0 = A.stride_0(),
            as1 = A.stride_1(),
            bs0 = B.stride_0(),
            bs1 = B.stride_1();
          
          for (int p=0;p<n;++p) {
            const int
              iend = m,
              jend = n-p-1;

            const value_type
              *__restrict__ a12t = jend ? &A(p, p+1) : NULL;

            value_type
              *__restrict__ b1   =        &B(0, p  ),
              *__restrict__ B2   = jend ? &B(0, p+1) : NULL;
                        
            if (!ArgDiag::use_unit_diag) {
              const value_type alpha11 = A(p, p);
              for (int i=0;i<iend;++i)
                b1[i*bs0] /= alpha11;
            }
            
            for (int i=0;i<iend;++i)
              for (int j=0;j<jend;++j)
                B2[i*bs0+j*bs1] -= a12t[j*as1] * b1[i*bs0];
          }
        }
        
        return 0;
      }
    };

    template<typename ArgDiag>
    struct Trsm<Side::Right,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Blocked> {
      
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {

        typedef typename BViewType::value_type value_type;        

        if (alpha == 0) {
          Util::set(B, 0);
        } else {
          if (alpha != 1)
            Util::scale(B, alpha);
          
          // B (m x n), A(m x m)
          const int
            m = B.dimension(0),
            n = B.dimension(1);

          if (!(m && n)) return 0;

          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            bs0 = B.stride_0(),
            bs1 = B.stride_1();

	  enum : int {
            mb = Algo::Trsm::Blocked::mb };
          
          InnerTrsmRightUpperUnitDiag<mb> trsm_ruu(as0, as1,
                                                       bs0, bs1);

          InnerTrsmRightUpperNonUnitDiag<mb> trsm_run(as0, as1,
                                                          bs0, bs1);

          InnerRankUpdate<mb,mb> gemm_nn(as0, as1,
                                         bs0, bs1,
                                         bs0, bs1);
          
          const int 
            mm = (m/mb)*mb,
            nn = (n/mb)*mb;

          for (int p=0;p<nn;p+=mb) {
            // trsm update
            if (ArgDiag::use_unit_diag) 
              trsm_ruu.invoke(&A(p,p), mm, &B(0,p));
            else
              trsm_run.invoke(&A(p,p), mm, &B(0,p));
            
            // gemm update
            for (int j=p+mb;j<nn;j+=mb)
              for (int i=0;i<mm;i+=mb)
                gemm_nn.serial_invoke(-1, &B(i,p), &A(p,j), mb, &B(i,j));
          }

          // remainder
          const int 
            mp = (m%mb),
            np = (n%mb);
          
          if (np) {
            gemm_nn.serial_invoke(-1, &B(0,0), &A(0,nn), mm, np, nn, &B(0,nn)); 
            if (ArgDiag::use_unit_diag) 
              trsm_ruu.invoke(&A(nn,nn), mm, np, &B(0,nn));
            else 
              trsm_run.invoke(&A(nn,nn), mm, np, &B(0,nn));
          }
          if (mp) {
            if (ArgDiag::use_unit_diag) 
              trsm_ruu.invoke(&A(0,0), mp, n, &B(mm,0));
            else 
              trsm_run.invoke(&A(0,0), mp, n, &B(mm,0));
          }

        }
        
        return 0;
      }
      
    };

    template<typename ArgDiag>
    struct Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::CompactMKL> {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {
        typedef typename BViewType::value_type vector_type;
        typedef typename vector_type::value_type value_type;
        
        const int
          m = B.dimension(0),
          n = B.dimension(1),
          vl = vector_type::vector_length;

        // no error check
        cblas_dtrsm_compact(CblasRowMajor, 
                            CblasLeft, CblasUpper, CblasNoTrans, 
                            ArgDiag::use_unit_diag ? CblasUnit : CblasNonUnit,
                            m, n, 
                            alpha, 
                            (const double*)A.data(), A.stride_0(), 
                            (double*)B.data(), B.stride_0(), 
                            (MKL_INT)vl, (MKL_INT)1);
      }
    };

    template<typename ArgDiag>
    struct Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Unblocked> {
      
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {

        typedef typename BViewType::value_type value_type;        

        if (alpha == 0) {
          Util::set(B, 0);
        } else {
          if (alpha != 1)
            Util::scale(B, alpha);
          
          // B (m x n), A(m x m)
          const int
            m = B.dimension(0),
            n = B.dimension(1);

          if (!(m && n)) return 0;

          const int
            as0 = A.stride_0(),
            //as1 = A.stride_1(),
            bs0 = B.stride_0(),
            bs1 = B.stride_1();

          value_type
            *__restrict__ B0 = &B(0, 0);
          
          for (int p=(m-1);p>=0;--p) {
            const int
              iend = p,
              jend = n;

            const value_type
              *__restrict__ a01 = &A(0, p);

            value_type
              *__restrict__ b1t = &B(p, 0);
            
            if (!ArgDiag::use_unit_diag) {
              const value_type alpha11 = A(p, p);
              for (int j=0;j<n;++j)
                b1t[j*bs1] /= alpha11;
            }
            
            for (int i=0;i<iend;++i)
              for (int j=0;j<jend;++j)
                B0[i*bs0+j*bs1] -= a01[i*as0] * b1t[j*bs1];
          }
        }        
        return 0;
      }
    };

    template<typename ArgDiag>
    struct Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,ArgDiag,Algo::Trsm::Blocked> {
      
      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType &A,
             const BViewType &B) {

        typedef typename BViewType::value_type value_type;        

        if (alpha == 0) {
          Util::set(B, 0);
        } else {
          if (alpha != 1)
            Util::scale(B, alpha);
          
          // B (m x n), A(m x m)
          const int
            m = B.dimension(0),
            n = B.dimension(1);

          if (!(m && n)) return 0;

          const int
            as0 = A.stride_0(),
            as1 = A.stride_1(),
            bs0 = B.stride_0(),
            bs1 = B.stride_1();

	  enum : int {
            mb = Algo::Trsm::Blocked::mb };
          
          InnerTrsmLeftUpperUnitDiag<mb> trsm_luu(as0, as1,
                                                      bs0, bs1);

          InnerTrsmLeftUpperNonUnitDiag<mb> trsm_lun(as0, as1,
                                                         bs0, bs1);

          InnerRankUpdate<mb,mb> gemm_nn(as0, as1,
                                         bs0, bs1,
                                         bs0, bs1);
          
          const int 
            mm = (m/mb)*mb,
            nn = (n/mb)*mb;

          for (int pp=0;pp<mm;pp+=mb) {
            const int p = m - pp - mb;

            // trsm update
            if (ArgDiag::use_unit_diag) 
              trsm_luu.invoke(&A(p,p), nn, &B(p,0));
            else
              trsm_lun.invoke(&A(p,p), nn, &B(p,0));
            
            // gemm update
            for (int i=(p-mb);i>=0;i-=mb)
              for (int j=0;j<nn;j+=mb)
                gemm_nn.serial_invoke(-1, &A(i,p), &B(p,j), mb, &B(i,j));
          }

          // remainder
          const int 
            mp = (m%mb),
            np = (n%mb);

          if (mp) {
            if (mm && nn) 
              gemm_nn.serial_invoke(-1, &A(0,mp), &B(mp,0), mp, nn, mm, &B(0,0)); 
            if (ArgDiag::use_unit_diag) 
              trsm_luu.invoke(&A(0,0), mp, nn, &B(0,0));
            else 
              trsm_lun.invoke(&A(0,0), mp, nn, &B(0,0));
            
          }
          if (np) {
            if (ArgDiag::use_unit_diag) 
              trsm_luu.invoke(&A(0,0), m, np, &B(0,nn));
            else 
              trsm_lun.invoke(&A(0,0), m, np, &B(0,nn));
          }
        }
        
        return 0;
      }
      
    };

  }


}

#endif
