#ifndef TEST_QUADRATURE
#define TEST_QUADRATURE

#include "Compadre_Quadrature.hpp"
#include <gtest/gtest.h>
#include <stdexcept>

using namespace Compadre;

// The quadrature tables cover orders 0-9 in 1D and 0-19 in 2D; requests outside those
// ranges must raise an error in every build type rather than read past the tables.

TEST (QuadratureTest, OneDimensionalOrderLimits) {
    for (int order = 0; order <= 9; ++order) {
        EXPECT_NO_THROW(Quadrature(order, 1, "LINE")) << "order " << order;
    }
    EXPECT_THROW(Quadrature(10, 1, "LINE"), std::logic_error);
    EXPECT_THROW(Quadrature(-1, 1, "LINE"), std::logic_error);
}

TEST (QuadratureTest, TwoDimensionalOrderLimits) {
    for (int order = 0; order <= 19; ++order) {
        EXPECT_NO_THROW(Quadrature(order, 2, "TRI")) << "order " << order;
    }
    EXPECT_THROW(Quadrature(20, 2, "TRI"), std::logic_error);
    EXPECT_THROW(Quadrature(-1, 2, "TRI"), std::logic_error);
}

TEST (QuadratureTest, OneDimensionalWeightsSumToOne) {
    // weights of each 1D rule are normalized to the unit interval
    for (int order = 0; order <= 9; ++order) {
        Quadrature q(order, 1, "LINE");
        double sum = 0.0;
        auto w = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), q.getWeights());
        for (int i = 0; i < q.getNumberOfQuadraturePoints(); ++i) sum += w(i);
        EXPECT_NEAR(1.0, sum, 1e-14) << "order " << order;
    }
}

#endif
