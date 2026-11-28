// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef TEST_TARGETS
#define TEST_TARGETS

#include "Compadre_GMLS.hpp"
#include "Compadre_Functors.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace Compadre;

class TargetTest: public ::testing::Test {
public:

    Kokkos::View<double*, Kokkos::DefaultExecutionSpace> epsilon_device;
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> source_coords_device;
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> target_coords_device;
    Kokkos::View<int**, Kokkos::DefaultExecutionSpace> neighbor_lists_device;
    Kokkos::View<double**, Kokkos::DefaultExecutionSpace> additional_target_coords_device;
    Kokkos::View<int**, Kokkos::DefaultExecutionSpace> additional_target_indices_device;

    TargetTest() {
    }

    // silence warning, letting compiler know we mean
    // to override
    using ::testing::Test::SetUp; 
    virtual void SetUp(const int dimension) {
        // initialization
        int num_source_sites = 3;
        if (dimension==2) num_source_sites = 9;

        source_coords_device = Kokkos::View<double**, Kokkos::DefaultExecutionSpace>("source coordinates", 
                num_source_sites, dimension);
        Kokkos::View<double**>::host_mirror_type source_coords = Kokkos::create_mirror_view(source_coords_device);
        if (dimension==1) {
            source_coords(0,0) = 0.0;
            source_coords(1,0) = 1.0;
            source_coords(2,0) = 2.0;
        } else if (dimension==2) {
            source_coords(0,0) = -1.0;
            source_coords(0,1) = -1.0;
            source_coords(1,0) =  0.0;
            source_coords(1,1) = -1.0;
            source_coords(2,0) =  1.0;
            source_coords(2,1) = -1.0;
            source_coords(3,0) = -1.0;
            source_coords(3,1) =  0.0;
            source_coords(4,0) =  0.0;
            source_coords(4,1) =  0.0;
            source_coords(5,0) =  1.0;
            source_coords(5,1) =  0.0;
            source_coords(6,0) = -1.0;
            source_coords(6,1) =  1.0;
            source_coords(7,0) =  0.0;
            source_coords(7,1) =  1.0;
            source_coords(8,0) =  1.0;
            source_coords(8,1) =  1.0;
        }
        Kokkos::fence();
        Kokkos::deep_copy(source_coords_device, source_coords);

        target_coords_device = Kokkos::View<double**, Kokkos::DefaultExecutionSpace>("target coordinates", 
                1, dimension);
        Kokkos::View<double**>::host_mirror_type target_coords = Kokkos::create_mirror_view(target_coords_device);
        if (dimension==1) {
            target_coords(0,0) = 1.0;
        } else if (dimension==2) {
            target_coords(0,0) = 0.0;
            target_coords(0,1) = 0.0;
        }
        Kokkos::fence();
        Kokkos::deep_copy(target_coords_device, target_coords);

        epsilon_device = Kokkos::View<double*, Kokkos::DefaultExecutionSpace>("h supports", 1);
        Kokkos::View<double*>::host_mirror_type epsilon = Kokkos::create_mirror_view(epsilon_device);
        if (dimension==1) {
            epsilon(0) = 1.0;
        } else if (dimension==2) {
            epsilon(0) = std::sqrt(2)+1e-10;
        }
        Kokkos::fence();
        Kokkos::deep_copy(epsilon_device, epsilon);

        neighbor_lists_device = Kokkos::View<int**, Kokkos::DefaultExecutionSpace>("neighbor lists", 
                1, num_source_sites+1); // first column is # of neighbors
        Kokkos::View<int**>::host_mirror_type neighbor_lists = Kokkos::create_mirror_view(neighbor_lists_device);
        if (dimension==1) {
            neighbor_lists(0,0) = 3;
            neighbor_lists(0,1) = 0;
            neighbor_lists(0,2) = 1;
            neighbor_lists(0,3) = 2;
        } else if (dimension==2) {
            neighbor_lists(0,0) = 9;
            neighbor_lists(0,1) = 0;
            neighbor_lists(0,2) = 1;
            neighbor_lists(0,3) = 2;
            neighbor_lists(0,4) = 3;
            neighbor_lists(0,5) = 4;
            neighbor_lists(0,6) = 5;
            neighbor_lists(0,7) = 6;
            neighbor_lists(0,8) = 7;
            neighbor_lists(0,9) = 8;
        }
        Kokkos::fence();
        Kokkos::deep_copy(neighbor_lists_device, neighbor_lists);

        additional_target_coords_device = Kokkos::View<double**, Kokkos::DefaultExecutionSpace>("additional target coordinates", 
                6, dimension);
        Kokkos::View<double**>::host_mirror_type additional_target_coords = Kokkos::create_mirror_view(additional_target_coords_device);
        if (dimension==1) {
            // out of bounds
            additional_target_coords(0,0) = -1.0;
            // out of bounds
            additional_target_coords(1,0) = 3.0;
            // in bounds
            additional_target_coords(2,0) = 0.5;
            // in bounds
            additional_target_coords(3,0) = 1.5;
            // on boundary
            additional_target_coords(4,0) = 0.0;
            // on boundary
            additional_target_coords(5,0) = 2.0;
        } else if (dimension==2) {
            // out of bounds
            additional_target_coords(0,0) =  -2.0;
            additional_target_coords(0,1) =  -2.0;
            // out of bounds
            additional_target_coords(1,0) =   3.0;
            additional_target_coords(1,1) =   3.0;
            // in bounds
            additional_target_coords(2,0) =   0.25;
            additional_target_coords(2,1) =  -0.25;
            // in bounds
            additional_target_coords(3,0) =  -0.5;
            additional_target_coords(3,1) =   0.5;
            // on boundary
            additional_target_coords(4,0) =   1.0;
            additional_target_coords(4,1) =   1.0;
            // on boundary
            additional_target_coords(5,0) =   1.0;
            additional_target_coords(5,1) =  -1.0;
        }
        Kokkos::fence();
        Kokkos::deep_copy(additional_target_coords_device, additional_target_coords);

        additional_target_indices_device = Kokkos::View<int**, Kokkos::DefaultExecutionSpace>("additional target indices", 
                1, 7 /* # of extra evaluation sites plus index for each */);
        Kokkos::View<int**>::host_mirror_type additional_target_indices = Kokkos::create_mirror_view(additional_target_indices_device);
        if (dimension==1) {
            additional_target_indices(0,0)=6;
            additional_target_indices(0,1)=0;
            additional_target_indices(0,2)=1;
            additional_target_indices(0,3)=2;
            additional_target_indices(0,4)=3;
            additional_target_indices(0,5)=4;
            additional_target_indices(0,6)=5;
        } else if (dimension==2) {
            additional_target_indices(0,0)=6;
            additional_target_indices(0,1)=0;
            additional_target_indices(0,2)=1;
            additional_target_indices(0,3)=2;
            additional_target_indices(0,4)=3;
            additional_target_indices(0,5)=4;
            additional_target_indices(0,6)=5;
        }
        Kokkos::fence();
        Kokkos::deep_copy(additional_target_indices_device, additional_target_indices);
    }

    void TearDown() {
    }
};

TEST_F (TargetTest, SupportedOrders) {
    // orders are checked (in every build type) when set, since basis evaluation does not check them
    EXPECT_NO_THROW(GMLS(ScalarTaylorPolynomial, PointSample, 10, 2, "QR"));
    EXPECT_THROW(GMLS(ScalarTaylorPolynomial, PointSample, 11, 2, "QR"), std::logic_error);
    EXPECT_THROW(GMLS(ScalarTaylorPolynomial, PointSample, PointSample, 2, 3, "QR", "MANIFOLD", "NO_CONSTRAINT", 11), std::logic_error);
    // StaggeredEdgeAnalyticGradientIntegralSample uses one order higher internally
    EXPECT_NO_THROW(GMLS(ScalarTaylorPolynomial, StaggeredEdgeAnalyticGradientIntegralSample, 9, 2, "QR"));
    EXPECT_THROW(GMLS(ScalarTaylorPolynomial, StaggeredEdgeAnalyticGradientIntegralSample, 10, 2, "QR"), std::logic_error);
    // setters check before changing anything
    GMLS gmls(ScalarTaylorPolynomial, PointSample, 2, 2, "QR");
    EXPECT_THROW(gmls.setPolynomialOrder(11), std::logic_error);
    EXPECT_EQ(gmls.getPolynomialOrder(), 2);
    EXPECT_THROW(gmls.setCurvaturePolynomialOrder(11), std::logic_error);
    GMLS gmls_staggered(ScalarTaylorPolynomial, StaggeredEdgeAnalyticGradientIntegralSample, 2, 2, "QR");
    EXPECT_THROW(gmls_staggered.setPolynomialOrder(10), std::logic_error);
}

TEST_F (TargetTest, RankDeficientQR) {
    // the 1D fixture has 3 neighbors (enough for order 2), but the outer two lie on the support boundary
    // and get zero weight, so P^T*W*P is singular; QR solves it (CHOL aborts, which cannot be tested here)
    int order = 2;
    int dimension = 1;
    SetUp(dimension);
    GMLS gmls(ScalarTaylorPolynomial, PointSample, order, dimension, "QR");
    gmls.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);
    gmls.addTargets(ScalarPointEvaluation);
    EXPECT_NO_THROW(gmls.generateAlphas());
}

TEST_F (TargetTest, NeighborCountCheck) {
    // the 1D fixture gives each target 3 neighbors: enough for an order 2 basis (3 unknowns), not order 3
    int dimension = 1;
    SetUp(dimension);
    // neighbor lists set after construction
    GMLS gmls(ScalarTaylorPolynomial, PointSample, 2, dimension, "QR");
    EXPECT_NO_THROW(gmls.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device));
    // order raised after the neighbor lists are set
    EXPECT_THROW(gmls.setPolynomialOrder(3), std::logic_error);
    // neighbor lists set when the order is already too high
    GMLS gmls_3(ScalarTaylorPolynomial, PointSample, 3, dimension, "QR");
    EXPECT_THROW(gmls_3.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device), std::logic_error);
    // an explicit opt-out accepts the lists, here and in later order changes
    GMLS gmls_opt_out(ScalarTaylorPolynomial, PointSample, 3, dimension, "QR");
    EXPECT_NO_THROW(gmls_opt_out.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device, false /*verify_unisolvency*/));
    EXPECT_NO_THROW(gmls_opt_out.setPolynomialOrder(4, false /*verify_unisolvency*/));
    EXPECT_THROW(gmls_opt_out.setPolynomialOrder(4), std::logic_error);
}

TEST_F (TargetTest, EvalBernsteinBasis1D) {

    int order = 2;
    int dimension = 1;

    SetUp(dimension);

    // this test checks basis evaluation only; neighbors on the support boundary have zero weight,
    // leaving the problem rank deficient, so use QR (CHOL requires a full rank problem)
    GMLS gmls(BernsteinPolynomial, VectorPointSample,
                 order, dimension, "QR");

    gmls.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);

    // set up additional sites to evaluate target operators
    gmls.setAdditionalEvaluationSitesData(additional_target_indices_device, additional_target_coords_device);
    
    // create a vector of target operations
    std::vector<TargetOperation> lro(1);
    lro[0] = ScalarPointEvaluation;
    
    // and then pass them to the GMLS class
    gmls.addTargets(lro);

    // solve to populate values
    // 1 batch, keep coefficients, don't clear cache
    gmls.generateAlphas(1, true, false);

    auto gmls_basis_data = gmls.extractBasisData();

    int local_index = 0;
    device_unmanaged_matrix_right_type P_target_row_device_unmanaged(gmls_basis_data.P_target_row_data 
            + TO_GLOBAL(local_index)*TO_GLOBAL(gmls_basis_data.P_target_row_dim_0*gmls_basis_data.P_target_row_dim_1), 
                gmls_basis_data.P_target_row_dim_0, gmls_basis_data.P_target_row_dim_1);
    auto P_target_row = Kokkos::create_mirror_view(P_target_row_device_unmanaged);
    Kokkos::deep_copy(P_target_row, P_target_row_device_unmanaged);

    int j=0; // operation
    int m=0; // tilesize
    int k=0; // tilesize
    SolutionSet<host_memory_space> h_ss(gmls_basis_data._d_ss);
    h_ss.copyAlphas(gmls_basis_data._d_ss);
    for (int e=0; e<7; ++e) {
        const int offset_index_jmke = h_ss.getTargetOffsetIndex(j,m,k,e);
        //printf("offset %d: %d\n", e, offset_index_jmke);
        //for (int i=0; i<order+1; ++i) {
        //    printf("P(%d,%d): %.16g\n", e, i, P_target_row(offset_index_jmke, i));
        //}
        // e==0 is target site
        if (e==0) {
            ASSERT_DOUBLE_EQ (P_target_row(offset_index_jmke, 0), 0.25);
            ASSERT_DOUBLE_EQ (P_target_row(offset_index_jmke, 1), 0.5);
            ASSERT_DOUBLE_EQ (P_target_row(offset_index_jmke, 2), 0.25);
        // e==1 is additional site out of bounds left, should have atleast one coefficient negative
        } else if (e==1) {
            bool found_negative = false;
            for (int i=0; i<order+1; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
            }
            ASSERT_TRUE(found_negative);
        // e==2 is additional site out of bounds right, should have atleast one coefficient negative
        } else if (e==2) {
            bool found_negative = false;
            for (int i=0; i<order+1; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
            }
            ASSERT_TRUE(found_negative);
        // e==3 is additional site in bounds left of target, should sum to 1.0 and have no negatives
        } else if (e==3) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<order+1; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        // e==4 is additional site in bounds right of target, should sum to 1.0 and have no negatives
        } else if (e==4) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<order+1; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        // e==5 is additional site at left bounds (given r=1), should sum to 1.0 and have no negatives
        } else if (e==5) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<order+1; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        // e==6 is additional site at right bounds (given r=1), should sum to 1.0 and have no negatives
        } else if (e==6) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<order+1; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        }
    }
}

TEST_F (TargetTest, EvalBernsteinBasis2D) {

    int order = 2;
    int dimension = 2;

    SetUp(dimension);

    // this test checks basis evaluation only; neighbors on the support boundary have zero weight,
    // leaving the problem rank deficient, so use QR (CHOL requires a full rank problem)
    GMLS gmls(BernsteinPolynomial, VectorPointSample,
                 order, dimension, "QR");

    gmls.setProblemData(neighbor_lists_device, source_coords_device, target_coords_device, epsilon_device);

    // set up additional sites to evaluate target operators
    gmls.setAdditionalEvaluationSitesData(additional_target_indices_device, additional_target_coords_device);
    
    // create a vector of target operations
    std::vector<TargetOperation> lro(1);
    lro[0] = ScalarPointEvaluation;
    
    // and then pass them to the GMLS class
    gmls.addTargets(lro);

    // solve to populate values
    // 1 batch, keep coefficients, don't clear cache
    gmls.generateAlphas(1, true, false);

    auto gmls_basis_data = gmls.extractBasisData();

    int local_index = 0;
    device_unmanaged_matrix_right_type P_target_row_device_unmanaged(gmls_basis_data.P_target_row_data 
            + TO_GLOBAL(local_index)*TO_GLOBAL(gmls_basis_data.P_target_row_dim_0*gmls_basis_data.P_target_row_dim_1), 
                gmls_basis_data.P_target_row_dim_0, gmls_basis_data.P_target_row_dim_1);
    auto P_target_row = Kokkos::create_mirror_view(P_target_row_device_unmanaged);
    Kokkos::deep_copy(P_target_row, P_target_row_device_unmanaged);

    int j=0; // operation
    int m=0; // tilesize
    int k=0; // tilesize
    SolutionSet<host_memory_space> h_ss(gmls_basis_data._d_ss);
    h_ss.copyAlphas(gmls_basis_data._d_ss);
    int basis_size = (order+1)*(order+1);
    for (int e=0; e<7; ++e) {
        const int offset_index_jmke = h_ss.getTargetOffsetIndex(j,m,k,e);
        //printf("offset %d: %d\n", e, offset_index_jmke);
        //for (int i=0; i<basis_size; ++i) {
        //    printf("P(%d,%d): %.16g\n", e, i, P_target_row(offset_index_jmke, i));
        //}
        // e==0 is target site
        if (e==0) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<basis_size; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ(sum, 1.0);
        // e==1 is additional site out of bounds left, should have atleast one coefficient negative
        } else if (e==1) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<basis_size; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(found_negative);
            ASSERT_DOUBLE_EQ(sum, 1.0);
        // e==2 is additional site out of bounds right, should have atleast one coefficient negative
        } else if (e==2) {
            bool found_negative = false;
            for (int i=0; i<basis_size; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
            }
            ASSERT_TRUE(found_negative);
        // e==3 is additional site in bounds left of target, should sum to 1.0 and have no negatives
        } else if (e==3) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<basis_size; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        // e==4 is additional site in bounds right of target, should sum to 1.0 and have no negatives
        } else if (e==4) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<basis_size; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        // e==5 is additional site (1,1), inside the square [-h,h]^2 the basis is defined on (h=sqrt(2)+1e-10), should sum to 1.0 and have no negatives
        } else if (e==5) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<basis_size; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        // e==6 is additional site (1,-1), inside the square [-h,h]^2 the basis is defined on (h=sqrt(2)+1e-10), should sum to 1.0 and have no negatives
        } else if (e==6) {
            bool found_negative = false;
            double sum = 0.0;
            for (int i=0; i<basis_size; ++i) {
                if (P_target_row(offset_index_jmke, i) < 0.0) found_negative = true;
                sum += P_target_row(offset_index_jmke, i);
            }
            ASSERT_TRUE(!found_negative);
            ASSERT_DOUBLE_EQ (sum, 1.0);
        }
    }
}

#endif
