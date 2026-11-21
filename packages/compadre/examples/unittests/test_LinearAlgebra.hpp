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
            const int rank, const bool A_is_tranpose, const bool B_is_transpose) {

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
                if (A_is_tranpose) {
                    host_unmanaged_matrix_right_type this_A(A.data() + mat_num*lda*nda, lda, nda);
                    for (size_t ind=0; ind<a_i.size(); ++ind) {
                        this_A(a_j[ind], a_i[ind]) = a_k[ind];
                    }
                } else {
                    host_unmanaged_matrix_right_type this_A(A.data() + mat_num*lda*nda, lda, nda);
                    for (size_t ind=0; ind<a_i.size(); ++ind) {
                        this_A(a_i[ind], a_j[ind]) = a_k[ind];
                    }
                }
                if (B_is_transpose) {
                    host_unmanaged_matrix_right_type this_B(B.data() + mat_num*ldb*ndb, ldb, ndb);
                    for (size_t ind=0; ind<b_i.size(); ++ind) {
                        this_B(b_j[ind], b_i[ind]) = b_k[ind];
                    }
                } else {
                    host_unmanaged_matrix_right_type this_B(B.data() + mat_num*ldb*ndb, ldb, ndb);
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
// Tested for A=LayoutRightA(LRA)
// Tested for B=LayoutRightA(LRB)
// Tested for X=LayoutRightA(LRX)
//
// with A and B potentially tranposed (AT, BT) or not (NAT, NBT)
//
// along with various sizes of lda, nda, ldb, and ndb, and NRHS
// 

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
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
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchSolveQR_Same_LDA_NDA_Larger_NRHS_BT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=4, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, true/*B is T*/);
    // RHS B was set up tranposed 4 x 3
    //for (int i=0; i<3; ++i) {
    //    for (int j=0; j<3; ++j) {
    //        printf("A(%d,%d): %.16f\n", i, j, A(i*N+j));
    //    }
    //}
    //for (int i=0; i<4; ++i) {
    //    for (int j=0; j<3; ++j) {
    //        printf("B(%d,%d): %.16f\n", i, j, B(i*N+j));
    //    }
    //}
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_BT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=7, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_BT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=13, ndb=7;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false,  B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_BT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=8, ndb=12;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_NBT) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
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
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_NBT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=4;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, false/*B is T*/);
    // QR expects layout left B, so ndb and ldb reverse ordered
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_NBT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=7, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_NBT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=7, ndb=13;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_NBT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=12, ndb=8;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, false/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_LLA_LLB_LRX_AT_BT) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
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
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_AT_BT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=4, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_AT_BT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=7, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_AT_BT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=13, ndb=7;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_AT_BT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=8, ndb=12;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, true/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    // ndb and ldb are switched because solution stored in B is not transposed (even though original B-RHS data was transposed)
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ndb, ldb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_AT_NBT) {
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
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
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Same_LDA_NDA_Larger_NRHS_AT_NBT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=3, ndb=4;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, false/*B is T*/);
    // QR expects layout left B, so ndb and ldb reverse ordered
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDA_NDA_AT_NBT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=7;
    int ldb=3, ndb=3;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_AT_NBT) {
    // lda and nda larger than M and N
    int M=3, N=3, NRHS=3, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=7, ndb=13;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                   0  -0.071428571428571  -0.035714285714286
    //                   0  -0.285714285714286  -0.142857142857143
    //                   0  -0.071428571428571  -0.535714285714286
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

TEST_F (LinearAlgebraTest, Square_FullRank_batchQRPivotingSolve_Larger_LDB_NDB_Larger_NRHS_AT_NBT) {
    int M=3, N=3, NRHS=4, num_matrices=2, rank=3;
    int lda=3, nda=3;
    int ldb=12, ndb=8;
    SetUp(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, rank, true/*A is T*/, false/*B is T*/);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // solution: X = [
    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
    //               ]
    host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
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
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
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
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_right,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
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
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_left,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
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
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                   0  -0.071428571428571  -0.035714285714286
//LLX!    //                   0  -0.285714285714286  -0.142857142857143
//LLX!    //                   0  -0.071428571428571  -0.535714285714286
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//LLX!    GMLS_LinearAlgebra::batchSolve<layout_left,layout_right,layout_left>(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
//LLX!    // solution: X = [
//LLX!    //                  0  -0.071428571428571  -0.035714285714286  -0.053571428571429
//LLX!    //                  0  -0.285714285714286  -0.142857142857143  -0.214285714285714
//LLX!    //                  0  -0.071428571428571  -0.535714285714286  -0.803571428571429
//LLX!    //               ]
//LLX!    host_unmanaged_matrix_left_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//     GMLS_LinearAlgebra::batchSolve(pm, A_d.data(), lda, nda, B_d.data(), ldb, ndb, M, N, NRHS, num_matrices);
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
//     host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//     host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//     host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//     host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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
//     host_unmanaged_matrix_right_type B2(B.data() + 1*ldb*ndb, ldb, ndb);
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

//
// Tests with arbitrary dense A and B (non-symmetric A, rectangular A, SPD A for Cholesky)
//
// a is M x N in row-major order. If A_is_transpose, A is stored as its N x M transpose.
// b is M x NRHS in row-major order (or, if implicit_RHS, a vector of length M stored
// contiguously from entry (0,0)). If B_is_transpose, B is stored as its NRHS x M transpose.
//
class LinearAlgebraDenseTest: public LinearAlgebraTest {
public:
    void SetUpDense(const int lda, const int nda, const int ldb, const int ndb,
            const int M, const int N, const int NRHS, const int num_matrices,
            const std::vector<double>& a, const bool A_is_transpose,
            const std::vector<double>& b, const bool B_is_transpose, const bool implicit_RHS) {

        Kokkos::resize(A, lda*nda*num_matrices);
        Kokkos::resize(A_d, lda*nda*num_matrices);
        Kokkos::deep_copy(A, 0.0);
        Kokkos::resize(B, ldb*ndb*num_matrices);
        Kokkos::resize(B_d, ldb*ndb*num_matrices);
        Kokkos::deep_copy(B, 0.0);

        for (int mat_num=0; mat_num<num_matrices; ++mat_num) {
            host_unmanaged_matrix_right_type this_A(A.data() + mat_num*lda*nda, lda, nda);
            for (int i=0; i<M; ++i) {
                for (int j=0; j<N; ++j) {
                    if (A_is_transpose) this_A(j,i) = a[i*N+j];
                    else                this_A(i,j) = a[i*N+j];
                }
            }
            if (implicit_RHS) {
                for (int i=0; i<M; ++i) B(mat_num*ldb*ndb + i) = b[i];
            } else {
                host_unmanaged_matrix_right_type this_B(B.data() + mat_num*ldb*ndb, ldb, ndb);
                for (int i=0; i<M; ++i) {
                    for (int j=0; j<NRHS; ++j) {
                        if (B_is_transpose) this_B(j,i) = b[i*NRHS+j];
                        else                this_B(i,j) = b[i*NRHS+j];
                    }
                }
            }
        }
        Kokkos::deep_copy(A_d, A);
        Kokkos::deep_copy(B_d, B);
        Kokkos::fence();
    }
};

TEST_F (LinearAlgebraDenseTest, Square_FullRank_NonSymmetricA_batchQRPivotingSolve_Larger_LDA_NDA_AT_NBT) {
    // non-symmetric A, so a transpose that does not actually transpose changes the answer
    int M=3, N=3, NRHS=3, num_matrices=2;
    int lda=7, nda=5;
    int ldb=3, ndb=3;
    // A = [[2, 1, 0], [0, 3, 1], [1, 0, 4]]
    std::vector<double> a = {2, 1, 0, 0, 3, 1, 1, 0, 4};
    // B = [[1, 2, 0], [0, 1, 0], [3, 0, 1]]
    std::vector<double> b = {1, 2, 0, 0, 1, 0, 3, 0, 1};
    SetUpDense(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, a, true, b, false, false);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // X = [[0.6, 0.8, 0.04], [-0.2, 0.4, -0.08], [0.6, -0.2, 0.24]]
    host_unmanaged_matrix_right_type X(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR( 0.600000000000000, X(0,0), 1e-13);
    EXPECT_NEAR( 0.800000000000000, X(0,1), 1e-13);
    EXPECT_NEAR( 0.040000000000000, X(0,2), 1e-13);
    EXPECT_NEAR(-0.200000000000000, X(1,0), 1e-13);
    EXPECT_NEAR( 0.400000000000000, X(1,1), 1e-13);
    EXPECT_NEAR(-0.080000000000000, X(1,2), 1e-13);
    EXPECT_NEAR( 0.600000000000000, X(2,0), 1e-13);
    EXPECT_NEAR(-0.200000000000000, X(2,1), 1e-13);
    EXPECT_NEAR( 0.240000000000000, X(2,2), 1e-13);
}

TEST_F (LinearAlgebraDenseTest, Square_FullRank_NonSymmetricA_batchQRPivotingSolve_Larger_LDA_NDA_Larger_NRHS_AT_BT) {
    // non-symmetric A and non-square B, both transposed, with unequal padding
    int M=3, N=3, NRHS=4, num_matrices=2;
    int lda=7, nda=5;
    int ldb=8, ndb=6;
    // A = [[2, 1, 0], [0, 3, 1], [1, 0, 4]]
    std::vector<double> a = {2, 1, 0, 0, 3, 1, 1, 0, 4};
    // B = [[1, 2, 0, 1], [0, 1, 0, 2], [3, 0, 1, 0]]
    std::vector<double> b = {1, 2, 0, 1, 0, 1, 0, 2, 3, 0, 1, 0};
    SetUpDense(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, a, true, b, true, false);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // X = [[0.6, 0.8, 0.04, 0.16], [-0.2, 0.4, -0.08, 0.68], [0.6, -0.2, 0.24, -0.04]]
    host_unmanaged_matrix_right_type X(B.data() + 1*ldb*ndb, ndb, ldb);
    EXPECT_NEAR( 0.600000000000000, X(0,0), 1e-13);
    EXPECT_NEAR( 0.800000000000000, X(0,1), 1e-13);
    EXPECT_NEAR( 0.040000000000000, X(0,2), 1e-13);
    EXPECT_NEAR( 0.160000000000000, X(0,3), 1e-13);
    EXPECT_NEAR(-0.200000000000000, X(1,0), 1e-13);
    EXPECT_NEAR( 0.400000000000000, X(1,1), 1e-13);
    EXPECT_NEAR(-0.080000000000000, X(1,2), 1e-13);
    EXPECT_NEAR( 0.680000000000000, X(1,3), 1e-13);
    EXPECT_NEAR( 0.600000000000000, X(2,0), 1e-13);
    EXPECT_NEAR(-0.200000000000000, X(2,1), 1e-13);
    EXPECT_NEAR( 0.240000000000000, X(2,2), 1e-13);
    EXPECT_NEAR(-0.040000000000000, X(2,3), 1e-13);
}

TEST_F (LinearAlgebraDenseTest, Rectangular_FullRank_batchQRPivotingSolve_ImplicitRHS_NAT) {
    // M > N least squares, as used by GMLS: solves A X = sqrt(W) in the least-squares sense
    int M=5, N=3, NRHS=5, num_matrices=2;
    int lda=5, nda=3;
    int ldb=7, ndb=7;
    // A = [[1, 0, 2], [0, 1, 1], [1, 1, 0], [2, 0, 1], [0, 3, 1]]
    std::vector<double> a = {1, 0, 2, 0, 1, 1, 1, 1, 0, 2, 0, 1, 0, 3, 1};
    // B is implicit: diagonal entries of sqrt(W) stored as a vector
    std::vector<double> b = {1, 2, 3, 1, 2};
    SetUpDense(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, a, false, b, false, true);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, true);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // X = [[-0.0883721, -0.288372, 0.976744, 0.381395, -0.12093], [-0.144186, 0.055814, 0.488372, -0.00930233, 0.539535], [0.418605, 0.418605, -0.837209, -0.0697674, 0.0465116]]
    host_unmanaged_matrix_right_type X(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(-0.088372093023256, X(0,0), 1e-13);
    EXPECT_NEAR(-0.288372093023256, X(0,1), 1e-13);
    EXPECT_NEAR( 0.976744186046511, X(0,2), 1e-13);
    EXPECT_NEAR( 0.381395348837209, X(0,3), 1e-13);
    EXPECT_NEAR(-0.120930232558140, X(0,4), 1e-13);
    EXPECT_NEAR(-0.144186046511628, X(1,0), 1e-13);
    EXPECT_NEAR( 0.055813953488372, X(1,1), 1e-13);
    EXPECT_NEAR( 0.488372093023256, X(1,2), 1e-13);
    EXPECT_NEAR(-0.009302325581395, X(1,3), 1e-13);
    EXPECT_NEAR( 0.539534883720930, X(1,4), 1e-13);
    EXPECT_NEAR( 0.418604651162791, X(2,0), 1e-13);
    EXPECT_NEAR( 0.418604651162790, X(2,1), 1e-13);
    EXPECT_NEAR(-0.837209302325581, X(2,2), 1e-13);
    EXPECT_NEAR(-0.069767441860465, X(2,3), 1e-13);
    EXPECT_NEAR( 0.046511627906976, X(2,4), 1e-13);
}

TEST_F (LinearAlgebraDenseTest, Rectangular_FullRank_batchQRPivotingSolve_ImplicitRHS_Larger_LDA_NDA_AT) {
    // M != N with A transposed: stored as N x M, transposed to M x N
    int M=5, N=3, NRHS=5, num_matrices=2;
    int lda=4, nda=6;
    int ldb=7, ndb=7;
    // A = [[1, 0, 2], [0, 1, 1], [1, 1, 0], [2, 0, 1], [0, 3, 1]]
    std::vector<double> a = {1, 0, 2, 0, 1, 1, 1, 1, 0, 2, 0, 1, 0, 3, 1};
    // B is implicit: diagonal entries of sqrt(W) stored as a vector
    std::vector<double> b = {1, 2, 3, 1, 2};
    SetUpDense(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, a, true, b, false, true);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::QR, A_d.data(), lda, nda, true, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, true);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // X = [[-0.0883721, -0.288372, 0.976744, 0.381395, -0.12093], [-0.144186, 0.055814, 0.488372, -0.00930233, 0.539535], [0.418605, 0.418605, -0.837209, -0.0697674, 0.0465116]]
    host_unmanaged_matrix_right_type X(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR(-0.088372093023256, X(0,0), 1e-13);
    EXPECT_NEAR(-0.288372093023256, X(0,1), 1e-13);
    EXPECT_NEAR( 0.976744186046511, X(0,2), 1e-13);
    EXPECT_NEAR( 0.381395348837209, X(0,3), 1e-13);
    EXPECT_NEAR(-0.120930232558140, X(0,4), 1e-13);
    EXPECT_NEAR(-0.144186046511628, X(1,0), 1e-13);
    EXPECT_NEAR( 0.055813953488372, X(1,1), 1e-13);
    EXPECT_NEAR( 0.488372093023256, X(1,2), 1e-13);
    EXPECT_NEAR(-0.009302325581395, X(1,3), 1e-13);
    EXPECT_NEAR( 0.539534883720930, X(1,4), 1e-13);
    EXPECT_NEAR( 0.418604651162791, X(2,0), 1e-13);
    EXPECT_NEAR( 0.418604651162790, X(2,1), 1e-13);
    EXPECT_NEAR(-0.837209302325581, X(2,2), 1e-13);
    EXPECT_NEAR(-0.069767441860465, X(2,3), 1e-13);
    EXPECT_NEAR( 0.046511627906976, X(2,4), 1e-13);
}

TEST_F (LinearAlgebraDenseTest, Square_SPD_batchCholeskySolve_Same_LDA_NDA_NBT) {
    // DenseSolverType::CHOL with full rank takes the Cholesky path
    int M=3, N=3, NRHS=3, num_matrices=2;
    int lda=3, nda=3;
    int ldb=3, ndb=3;
    // A = [[4, 1, 0], [1, 3, 1], [0, 1, 2]]
    std::vector<double> a = {4, 1, 0, 1, 3, 1, 0, 1, 2};
    // B = [[1, 2, 0], [0, 1, 0], [3, 0, 1]]
    std::vector<double> b = {1, 2, 0, 0, 1, 0, 3, 0, 1};
    SetUpDense(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, a, false, b, false, false);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::CHOL, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // X = [[0.444444, 0.444444, 0.0555556], [-0.777778, 0.222222, -0.222222], [1.88889, -0.111111, 0.611111]]
    host_unmanaged_matrix_right_type X(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR( 0.444444444444444, X(0,0), 1e-13);
    EXPECT_NEAR( 0.444444444444444, X(0,1), 1e-13);
    EXPECT_NEAR( 0.055555555555556, X(0,2), 1e-13);
    EXPECT_NEAR(-0.777777777777778, X(1,0), 1e-13);
    EXPECT_NEAR( 0.222222222222222, X(1,1), 1e-13);
    EXPECT_NEAR(-0.222222222222222, X(1,2), 1e-13);
    EXPECT_NEAR( 1.888888888888889, X(2,0), 1e-13);
    EXPECT_NEAR(-0.111111111111111, X(2,1), 1e-13);
    EXPECT_NEAR( 0.611111111111111, X(2,2), 1e-13);
}

TEST_F (LinearAlgebraDenseTest, Square_SPD_batchCholeskySolve_Larger_LDA_NDA_Larger_NRHS_BT) {
    // B transposed, as in every CHOL call in GMLS.cpp (P -> P^T)
    int M=3, N=3, NRHS=4, num_matrices=2;
    int lda=7, nda=5;
    int ldb=8, ndb=6;
    // A = [[4, 1, 0], [1, 3, 1], [0, 1, 2]]
    std::vector<double> a = {4, 1, 0, 1, 3, 1, 0, 1, 2};
    // B = [[1, 2, 0, 1], [0, 1, 0, 2], [3, 0, 1, 0]]
    std::vector<double> b = {1, 2, 0, 1, 0, 1, 0, 2, 3, 0, 1, 0};
    SetUpDense(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, a, false, b, true, false);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::CHOL, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, true, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // X = [[0.444444, 0.444444, 0.0555556, 0.0555556], [-0.777778, 0.222222, -0.222222, 0.777778], [1.88889, -0.111111, 0.611111, -0.388889]]
    host_unmanaged_matrix_right_type X(B.data() + 1*ldb*ndb, ndb, ldb);
    EXPECT_NEAR( 0.444444444444444, X(0,0), 1e-13);
    EXPECT_NEAR( 0.444444444444444, X(0,1), 1e-13);
    EXPECT_NEAR( 0.055555555555556, X(0,2), 1e-13);
    EXPECT_NEAR( 0.055555555555556, X(0,3), 1e-13);
    EXPECT_NEAR(-0.777777777777778, X(1,0), 1e-13);
    EXPECT_NEAR( 0.222222222222222, X(1,1), 1e-13);
    EXPECT_NEAR(-0.222222222222222, X(1,2), 1e-13);
    EXPECT_NEAR( 0.777777777777778, X(1,3), 1e-13);
    EXPECT_NEAR( 1.888888888888889, X(2,0), 1e-13);
    EXPECT_NEAR(-0.111111111111111, X(2,1), 1e-13);
    EXPECT_NEAR( 0.611111111111111, X(2,2), 1e-13);
    EXPECT_NEAR(-0.388888888888889, X(2,3), 1e-13);
}

TEST_F (LinearAlgebraDenseTest, Square_SPD_batchCholeskySolve_Larger_LDB_NDB_Larger_NRHS_NBT) {
    int M=3, N=3, NRHS=4, num_matrices=2;
    int lda=3, nda=3;
    int ldb=5, ndb=9;
    // A = [[4, 1, 0], [1, 3, 1], [0, 1, 2]]
    std::vector<double> a = {4, 1, 0, 1, 3, 1, 0, 1, 2};
    // B = [[1, 2, 0, 1], [0, 1, 0, 2], [3, 0, 1, 0]]
    std::vector<double> b = {1, 2, 0, 1, 0, 1, 0, 2, 3, 0, 1, 0};
    SetUpDense(lda, nda, ldb, ndb, M, N, NRHS, num_matrices, a, false, b, false, false);
    GMLS_LinearAlgebra::batchSolve(pm, DenseSolverType::CHOL, A_d.data(), lda, nda, false, B_d.data(), ldb, ndb, false, M, N, NRHS, num_matrices, false);
    Kokkos::deep_copy(B, B_d);
    Kokkos::fence();
    // X = [[0.444444, 0.444444, 0.0555556, 0.0555556], [-0.777778, 0.222222, -0.222222, 0.777778], [1.88889, -0.111111, 0.611111, -0.388889]]
    host_unmanaged_matrix_right_type X(B.data() + 1*ldb*ndb, ldb, ndb);
    EXPECT_NEAR( 0.444444444444444, X(0,0), 1e-13);
    EXPECT_NEAR( 0.444444444444444, X(0,1), 1e-13);
    EXPECT_NEAR( 0.055555555555556, X(0,2), 1e-13);
    EXPECT_NEAR( 0.055555555555556, X(0,3), 1e-13);
    EXPECT_NEAR(-0.777777777777778, X(1,0), 1e-13);
    EXPECT_NEAR( 0.222222222222222, X(1,1), 1e-13);
    EXPECT_NEAR(-0.222222222222222, X(1,2), 1e-13);
    EXPECT_NEAR( 0.777777777777778, X(1,3), 1e-13);
    EXPECT_NEAR( 1.888888888888889, X(2,0), 1e-13);
    EXPECT_NEAR(-0.111111111111111, X(2,1), 1e-13);
    EXPECT_NEAR( 0.611111111111111, X(2,2), 1e-13);
    EXPECT_NEAR(-0.388888888888889, X(2,3), 1e-13);
}

#endif
