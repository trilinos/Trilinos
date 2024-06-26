// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_INTERNAL_COMPADRE_HPP__
#define __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_INTERNAL_COMPADRE_HPP__


/// \author Paul Kuberry (pakuber@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Gemv_TeamVector_Internal.hpp"
#include "KokkosBatched_Trsv_TeamVector_Internal.hpp"

#include "KokkosBatched_Gemm_TeamVector_Internal.hpp"
#include "KokkosBatched_Trsm_TeamVector_Internal.hpp"

namespace KokkosBatched {

    /// TeamVector Internal
    /// =================== 
    //
    // Attention!: This is a fork of TeamVectorSolveUTV_Internal from
    // KokkosKernels that takes advantage of all non-square matrices
    // having a RHS that is diagonal. 
    //
    // This prevents having to actually form U'*B which would be very 
    // large if m >> n for the A matrix (result is m x m)
    //
    struct TeamVectorSolveUTV_Internal_Compadre {

    template<typename MemberType,
             typename ValueType,
         typename IntType>
    KOKKOS_INLINE_FUNCTION
    static int
    invoke(const MemberType &member, 
           const int matrix_rank,
           const int m, const int n, const int nrhs,
           const ValueType * U, const int us0, const int us1,
           const ValueType * T, const int ts0, const int ts1,
           const ValueType * V, const int vs0, const int vs1,
           const IntType   * p, const int ps0,
           /* */ ValueType * B, const int bs0, const int bs1,
           /* */ ValueType * X, const int xs0, const int xs1,
           /* */ ValueType * w, ValueType * wq, 
           const bool implicit_RHS) {
    
        typedef ValueType value_type;
    
        const value_type one(1), zero(0);
    
        value_type * W = w; /// m x nrhs
        value_type * WQ = wq; /// 3m
        const int ws0 = nrhs, ws1=1; // only works with w layout right
    
        bool do_print = false;
        if (do_print) {
            Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                using Kokkos::printf;
#endif
                printf("size is: %d %d %d %d\n", matrix_rank, m, n, nrhs);
                printf("U, us1, us0: %d %d\n", us1, us0);
                printf("T, ts0, ts1: %d %d\n", ts0, ts1);
                printf("B, bs0, bs1: %d %d\n", bs0, bs1);
                printf("W, ws0, ws1: %d %d\n", ws0, ws1);
                printf("B=zeros(%d,%d);\n", m, nrhs);
                for (int i=0; i<m; ++i) {
                    for (int j=0; j<nrhs; ++j) {
                        printf("B(%d,%d)= %f;\n", i+1,j+1,B[i*bs0+j*bs1]);
                    }
                }
            });
        }
    
        if (matrix_rank < n) {
            // U is m x matrix_rank
            // T is matrix_rank x matrix_rank
            // V is matrix_rank x n
            // W = U^T B
            if (!implicit_RHS) { // LU case
                TeamVectorGemmInternal<Algo::Gemm::Unblocked>
                  ::invoke(member,
                       matrix_rank, nrhs, m,
                       one,
                       U, us1, us0,
                       B, bs0, bs1,
                       zero,
                       W, ws0, ws1);
            } else {
                Kokkos::parallel_for
                (Kokkos::ThreadVectorRange(member,nrhs),
                 [&](const int &j) {
                  Kokkos::parallel_for
                    (Kokkos::TeamThreadRange(member,matrix_rank),
                      [&](const int &i) {                                   
                        W[i*ws0+j*ws1] = U[j*us0+i*us1] * B[j];
                  });
    
               });
    
            }
            member.team_barrier();
    
            if (do_print) {
                Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                    using Kokkos::printf;
#endif
                    printf("W=zeros(%d,%d);\n", m, nrhs);
                    for (int i=0; i<m; ++i) {
                        for (int j=0; j<nrhs; ++j) {
                            printf("W(%d,%d)= %f;\n", i+1,j+1,W[i*ws0+j]);
                        }
                    }
                      printf("B=zeros(%d,%d);\n", m, nrhs);
                    for (int i=0; i<m; ++i) {
                        for (int j=0; j<nrhs; ++j) {
                            printf("B(%d,%d)= %f;\n", i+1,j+1,B[i*bs0+j]);
                        }
                    }
                });
            }
    
            /// W = T^{-1} W
            TeamVectorTrsmInternalLeftLower<Algo::Trsm::Unblocked>
              ::invoke(member,
                   false,
                   matrix_rank, nrhs,
                   one,
                   T, ts0, ts1,
                   W, ws0, ws1);
            member.team_barrier();
    
            if (do_print) {
                Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                    using Kokkos::printf;
#endif
                    printf("W=zeros(%d,%d);\n", m, nrhs);
                    for (int i=0; i<m; ++i) {
                        for (int j=0; j<nrhs; ++j) {
                            printf("W(%d,%d)= %f;\n", i+1,j+1,W[i*ws0+j]);
                        }
                    }
                });
            }
            
            /// X = V^T W
            TeamVectorGemmInternal<Algo::Gemm::Unblocked>
              ::invoke(member,
                   n, nrhs, matrix_rank, 
                   one,
                   V, vs1, vs0,
                   W, ws0, ws1,
                   zero,
                   X, xs0, xs1);
            member.team_barrier();
    
            if (do_print) {
                Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                    using Kokkos::printf;
#endif
                    printf("X=zeros(%d,%d);\n", n, nrhs);
                    for (int i=0; i<n; ++i) {
                        for (int j=0; j<nrhs; ++j) {
                            printf("X(%d,%d)= %f;\n", i+1,j+1,X[i*xs0+j*xs1]);
                        }
                    }
                });
            }
        } else {
            if (!implicit_RHS) { // LU case
                /// W = U^T B
                TeamVectorGemmInternal<Algo::Gemm::Unblocked>
                  ::invoke(member,
                       matrix_rank, nrhs, m, 
                       one,
                       U, us1, us0,
                       B, bs0, bs1,
                       zero,
                       W, ws0, ws1);
                member.team_barrier();
                // copy W to X
                Kokkos::parallel_for
                (Kokkos::ThreadVectorRange(member,nrhs),
                 [&](const int &j) {
                  Kokkos::parallel_for
                    (Kokkos::TeamThreadRange(member,matrix_rank),
                      [&](const int &i) {                                   
                        X[i*xs0+j*xs1] = W[i*ws0+j*ws1];
                  });
               });
            } else {
                Kokkos::parallel_for
                (Kokkos::TeamVectorRange(member,nrhs),
                 [&](const int &i) {
                    WQ[i] = B[i];    
                });
                member.team_barrier();

                Kokkos::parallel_for
                (Kokkos::ThreadVectorRange(member,nrhs),
                 [&](const int &j) {
                  Kokkos::parallel_for
                    (Kokkos::TeamThreadRange(member,matrix_rank),
                      [&](const int &i) {                                   
                        X[i*xs0+j*xs1] = U[j*us0+i*us1] * WQ[j];
                  });
               });
            }
            member.team_barrier();
    
            if (do_print) {
                Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                    using Kokkos::printf;
#endif
                    printf("m=zeros(%d,%d);\n", matrix_rank, nrhs);
                    for (int i=0; i<matrix_rank; ++i) {
                        for (int j=0; j<nrhs; ++j) {
                            printf("m(%d,%d)= %f;\n", i+1,j+1,X[i*ws0+j]);
                        }
                    }
                });
            }
    
            if (do_print) {
                Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                    using Kokkos::printf;
#endif
                    printf("T=zeros(%d,%d);\n", m, matrix_rank);
                    for (int i=0; i<m; ++i) {
                        for (int j=0; j<matrix_rank; ++j) {
                            printf("T(%d,%d)= %f;\n", i+1,j+1,T[i*ts0+j]);
                        }
                    }
                });
            }
    
            /// X = T^{-1} X
            TeamVectorTrsmInternalLeftUpper<Algo::Trsm::Unblocked>
              ::invoke(member,
                   false,
                   matrix_rank, nrhs,
                   one,
                   T, ts0, ts1,
                   X, xs0, xs1);        
            member.team_barrier();
    
            if (do_print) {
                Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                    using Kokkos::printf;
#endif
                    printf("x=zeros(%d,%d);\n", n, nrhs);
                    for (int i=0; i<n; ++i) {
                        for (int j=0; j<nrhs; ++j) {
                            printf("x(%d,%d)= %f;\n", i+1,j+1,X[i*xs0+j*xs1]);
                        }
                    }
                });
            }
        }
          
        // X = P^T X
        TeamVectorApplyPivotMatrixBackwardInternal
            ::invoke(member,
                 nrhs, matrix_rank,
                 p, ps0,
                 X, xs0, xs1);
        if (do_print) {
            Kokkos::single(Kokkos::PerTeam(member), [&] () {
#if KOKKOS_VERSION >= 40200
                using Kokkos::printf;
#endif
                printf("X=zeros(%d,%d);\n", n, nrhs);
                for (int i=0; i<n; ++i) {
                    for (int j=0; j<nrhs; ++j) {
                        printf("X(%d,%d)= %f;\n", i+1,j+1,X[i*xs0+j*xs1]);
                    }
                }
            });
        }
    
        return 0;
    }
    };

} // end namespace KokkosBatched


#endif
