// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef TEST_LINEARALGEBRA
#define TEST_LINEARALGEBRA

#include "Compadre_LinearAlgebra_Definitions.hpp"
#include "Compadre_LinearAlgebra_Declarations.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace Compadre;

class LinearAlgebraTest: public ::testing::Test {
public:

    ParallelManager pm;
    Kokkos::View<double*, host_execution_space> A;
    Kokkos::View<double*, device_execution_space> A_d;
    Kokkos::View<double*, host_execution_space> B;
    Kokkos::View<double*, device_execution_space> B_d;

    LinearAlgebraTest( ) {
        // initialization
        pm = ParallelManager();
    }

    // silence warning, letting compiler know we mean
    // to override
    using ::testing::Test::SetUp; 
    virtual void SetUp(const int lda, const int nda, const int ldb, const int ndb, 
            const int M, const int N, const int NRHS, const int num_matrices, 
            const int rank, const bool A_is_LR, const bool B_is_LR) {

        // multiple full rank matrices
        Kokkos::resize(A, lda*nda*num_matrices);
        Kokkos::resize(A_d, lda*nda*num_matrices);
        Kokkos::deep_copy(A, 0.0);

        Kokkos::resize(B, ldb*ndb*num_matrices);
        Kokkos::resize(B_d, ldb*ndb*num_matrices);
        Kokkos::deep_copy(B, 0.0);

        if (M==N) { // square case
            std::vector<int> a_i, a_j, b_i, b_j;
            std::vector<double> a_k, b_k;
            if (rank==M) { // full rank
                // fill in A matrices
                // A = [[-4, 1, 0], [1, -4, 1], [0, 1 -4]]
                a_i ={0,0,1,1,1,2,2};
                a_j ={0,1,0,1,2,1,2};
                a_k ={-4,1,1,-4,1,1,-4};

            } else if (rank==N-1) { // rank deficient by one
                // fill in A matrices
                // A= [[-4, -8, 0], [1, 2, 1], [0, 0 -4]]
                a_i ={0,0,1,1,1,2};
                a_j ={0,1,0,1,2,2};
                a_k ={-4,-8,1,2,1,-4};
            }
            if (NRHS==M) {
                // B = [[0, 0, 0], [0, 1, 0], [0, 0, 2]]
                b_i ={0,1,2};
                b_j ={0,1,2};
                b_k ={0,1,2};
            } else if (NRHS==M+1) {
                // B = [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 3]]
                b_i ={0,1,2,2};
                b_j ={0,1,2,3};
                b_k ={0,1,2,3};
            }
            for (int mat_num=0; mat_num<num_matrices; ++mat_num) {
                if (A_is_LR) {
                    host_scratch_matrix_right_type this_A(A.data() + mat_num*lda*nda, lda, nda);
                    for (size_t ind=0; ind<a_i.size(); ++ind) {
                        this_A(a_i[ind], a_j[ind]) = a_k[ind];
                    }
                } else {
                    host_scratch_matrix_left_type this_A(A.data() + mat_num*lda*nda, lda, nda);
                    for (size_t ind=0; ind<a_i.size(); ++ind) {
                        this_A(a_i[ind], a_j[ind]) = a_k[ind];
                    }
                }
                if (B_is_LR) {
                    host_scratch_matrix_right_type this_B(B.data() + mat_num*ldb*ndb, ldb, ndb);
                    for (size_t ind=0; ind<b_i.size(); ++ind) {
                        this_B(b_i[ind], b_j[ind]) = b_k[ind];
                    }
                } else {
                    host_scratch_matrix_left_type this_B(B.data() + mat_num*ldb*ndb, ldb, ndb);
                    for (size_t ind=0; ind<b_i.size(); ++ind) {
                        this_B(b_i[ind], b_j[ind]) = b_k[ind];
                    }
                }
            }
        }
        Kokkos::deep_copy(A_d, A);
        Kokkos::deep_copy(B_d, B);
        Kokkos::fence();
    }

    void TearDown( ) {
        // after test completes
    }
};


//
// Square FullRank Tests
// Tested for A=LayoutLeft (LLA) LayoutRightA(LRA)
// Tested for B=LayoutLeft (LLB) LayoutRightA(LRB)
// Tested for X=LayoutLeft (LLX) LayoutRightA(LRX)
//
// along with various sizes of lda, nda, ldb, and ndb, and NRHS
// 

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LRA_LLB_LRX) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    //for (int i=0; i<3; ++i) {
    //    for (int j=0; j<3; ++j) {
    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
    //    }
    //}
    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LRA_LLB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=4;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    // LU expects layout left B, so ndb and ldb reverse ordered
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LRA_LLB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=7, nda=12;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LRA_LLB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=7, ndb=13;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LRA_LLB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=12, ndb=8;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LRA_LRB_LRX) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    //for (int i=0; i<3; ++i) {
    //    for (int j=0; j<3; ++j) {
    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
    //    }
    //}
    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LRA_LRB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=4;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    // QR expects layout left B, so ndb and ldb reverse ordered
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LRA_LRB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=7, nda=12;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LRA_LRB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=7, ndb=13;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LRA_LRB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=12, ndb=8;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LLA_LLB_LRX) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    //for (int i=0; i<3; ++i) {
    //    for (int j=0; j<3; ++j) {
    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
    //    }
    //}
    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LLA_LLB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=4;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    // LU expects layout left B, so ndb and ldb reverse ordered
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LLA_LLB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=7, nda=12;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LLA_LLB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=7, ndb=13;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LLA_LLB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=12, ndb=8;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LLA_LRB_LRX) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    //for (int i=0; i<3; ++i) {
    //    for (int j=0; j<3; ++j) {
    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
    //    }
    //}
    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LLA_LRB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=4;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    // QR expects layout left B, so ndb and ldb reverse ordered
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LLA_LRB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=7, nda=12;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LLA_LRB_LRX) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=7, ndb=13;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
}

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LLA_LRB_LRX) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=12, ndb=8;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_right>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
}

//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LRA_LLB_LLX) {
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    //for (int i=0; i<3; ++i) {
//LLX!    //    for (int j=0; j<3; ++j) {
//LLX!    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
//LLX!    //    }
//LLX!    //}
//LLX!    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LRA_LLB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=4;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    // LU expects layout left B, so ndb and ldb reverse ordered
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LRA_LLB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=7, nda=12;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LRA_LLB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=7, ndb=13;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LRA_LLB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=12, ndb=8;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LRA_LRB_LLX) {
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    //for (int i=0; i<3; ++i) {
//LLX!    //    for (int j=0; j<3; ++j) {
//LLX!    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
//LLX!    //    }
//LLX!    //}
//LLX!    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LRA_LRB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=4;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    // QR expects layout left B, so ndb and ldb reverse ordered
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LRA_LRB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=7, nda=12;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LRA_LRB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=7, ndb=13;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LRA_LRB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=12, ndb=8;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LLA_LLB_LLX) {
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    //for (int i=0; i<3; ++i) {
//LLX!    //    for (int j=0; j<3; ++j) {
//LLX!    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
//LLX!    //    }
//LLX!    //}
//LLX!    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LLA_LLB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=4;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    // LU expects layout left B, so ndb and ldb reverse ordered
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LLA_LLB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=7, nda=12;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LLA_LLB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=7, ndb=13;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LLA_LLB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=12, ndb=8;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, false/*B is LL*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LLA_LRB_LLX) {
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    //for (int i=0; i<3; ++i) {
//LLX!    //    for (int j=0; j<3; ++j) {
//LLX!    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
//LLX!    //    }
//LLX!    //}
//LLX!    //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_LLA_LRB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=3, ndb=4;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    // QR expects layout left B, so ndb and ldb reverse ordered
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_LLA_LRB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=7, nda=12;
//LLX!    int ldb=3, ndb=3;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_LLA_LRB_LLX) {
//LLX!    // lda and nda larger than M and N
//LLX!    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=7, ndb=13;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!}
//LLX!
//LLX!TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_LLA_LRB_LLX) {
//LLX!    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//LLX!    int lda=3, nda=3;
//LLX!    int ldb=12, ndb=8;
//LLX!    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is LR*/, true/*B is LR*/);
//LLX!    GMLS_LinearAlgebra::batchQRPivotingSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_scratch_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//LLX!    EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//LLX!    EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//LLX!    EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//LLX!    EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//LLX!    EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
//LLX!}

// TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA) {
//     int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//     int lda=3, nda=3;
//     int ldb=3, ndb=3; // relative to layout right
//     SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank);
//     // LU expects layout left B, so ndb and ldb reverse ordered
//     GMLS_LinearAlgebra::batchQRPivotingSolve(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//     // solution: X = [
//     //                   0  -0.071428571428571  -0.035714285714286
//     //                   0  -0.285714285714286  -0.142857142857143
//     //                   0  -0.071428571428571  -0.535714285714286
//     //               ]
//     //for (int i=0; i<3; ++i) {
//     //    for (int j=0; j<3; ++j) {
//     //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
//     //    }
//     //}
//     //for (int i=0; i<ldb*ndb; ++i) printf("Bi(%d): %.16f\n", i, B(i));
//     host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//     EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//     EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//     EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//     EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//     EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
// }
// 
// TEST_F (LinearAlgebraTest, Square_FullRank_batchQRFactorize_Same_LDA_NDA_Larger_NRHS) {
//     int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//     int lda=3, nda=3;
//     int ldb=4, ndb=3; // relative to layout right
//     SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank);
//     // LU expects layout left B, so ndb and ldb reverse ordered
//     GMLS_LinearAlgebra::batchQRFactorize(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//     // solution: X = [
//     //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//     //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//     //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//     //               ]
//     host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//     EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//     EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//     EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//     EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//     EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//     EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//     EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//     EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
// }
// 
// TEST_F (LinearAlgebraTest, Square_FullRank_batchQRFactorize_Larger_LDA_NDA) {
//     // lda and nda larger than M and N
//     int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//     int lda=7, nda=12;
//     int ldb=3, ndb=3; // relative to layout right
//     SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank);
//     // LU expects layout left B, so ndb and ldb reverse ordered
//     GMLS_LinearAlgebra::batchQRFactorize(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//     // solution: X = [
//     //                   0  -0.071428571428571  -0.035714285714286
//     //                   0  -0.285714285714286  -0.142857142857143
//     //                   0  -0.071428571428571  -0.535714285714286
//     //               ]
//     host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//     EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//     EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//     EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//     EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//     EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
// }
// 
// TEST_F (LinearAlgebraTest, Square_FullRank_batchQRFactorize_Larger_LDB_NDB) {
//     // lda and nda larger than M and N
//     int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
//     int lda=3, nda=3;
//     int ldb=7, ndb=13; // relative to layout right
//     SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank);
//     // LU expects layout left B, so ndb and ldb reverse ordered
//     GMLS_LinearAlgebra::batchQRFactorize(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//     // solution: X = [
//     //                   0  -0.071428571428571  -0.035714285714286
//     //                   0  -0.285714285714286  -0.142857142857143
//     //                   0  -0.071428571428571  -0.535714285714286
//     //               ]
//     host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//     EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//     EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//     EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//     EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//     EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
// }
// 
// TEST_F (LinearAlgebraTest, Square_FullRank_batchQRFactorize_Larger_LDB_NDB_Larger_NRHS) {
//     int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
//     int lda=3, nda=3;
//     int ldb=12, ndb=8; // relative to layout right
//     SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank);
//     // LU expects layout left B, so ndb and ldb reverse ordered
//     GMLS_LinearAlgebra::batchQRFactorize(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//     // solution: X = [
//     //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//     //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//     //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//     //               ]
//     host_scratch_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
//     EXPECT_NEAR(               0.0, B2(0,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(1,0), 1e-14);
//     EXPECT_NEAR(               0.0, B2(2,0), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(0,1), 1e-14);
//     EXPECT_NEAR(-0.285714285714286, B2(1,1), 1e-14);
//     EXPECT_NEAR(-0.071428571428571, B2(2,1), 1e-14);
//     EXPECT_NEAR(-0.035714285714286, B2(0,2), 1e-14);
//     EXPECT_NEAR(-0.142857142857143, B2(1,2), 1e-14);
//     EXPECT_NEAR(-0.535714285714286, B2(2,2), 1e-14);
//     EXPECT_NEAR(-0.053571428571429, B2(0,3), 1e-14);
//     EXPECT_NEAR(-0.214285714285714, B2(1,3), 1e-14);
//     EXPECT_NEAR(-0.803571428571429, B2(2,3), 1e-14);
// }

//TEST_F (LinearAlgebraTest, SquareRankDeficient) {
//}
//TEST_F (LinearAlgebraTest, RectangularFullRank) {
//}
//TEST_F (LinearAlgebraTest, RectangularRankDeficient) {
//}

#endif
