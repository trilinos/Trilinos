// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef TEST_XYZ
#define TEST_XYZ

#include "Compadre_Misc.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace Compadre;

class XyzTest: public ::testing::Test {
public:
    XYZ a;
    XYZ b;
    XYZ c;

    XyzTest( ) {
        // initialization
    }

    void SetUp( ) {
        // before test
        a = XYZ(1,1,1);
        b = XYZ(2,2,2);
        c = XYZ(1,2,3);
    }

    void TearDown( ) {
        // after test completes
    }
};

TEST_F (XyzTest, XYZMultiplyScalar) {
    ASSERT_DOUBLE_EQ (3.0, (a*3.0).x);
    ASSERT_DOUBLE_EQ (4.0, (a*4.0).y);
    ASSERT_DOUBLE_EQ (5.0, (a*5.0).y);
}
TEST_F (XyzTest, ScalarMultiplyXYZ) {
    ASSERT_DOUBLE_EQ (3.0, (3.0*a).x);
    ASSERT_DOUBLE_EQ (4.0, (4.0*a).y);
    ASSERT_DOUBLE_EQ (5.0, (5.0*a).y);
}
TEST_F (XyzTest, XYZAddXYZ) {
    auto ans = a+c;
    ASSERT_DOUBLE_EQ (2.0, ans.x);
    ASSERT_DOUBLE_EQ (3.0, ans.y);
    ASSERT_DOUBLE_EQ (4.0, ans.z);
}
TEST_F (XyzTest, XYZAddScalar) {
    ASSERT_DOUBLE_EQ (2.0, (a+1).x);
    ASSERT_DOUBLE_EQ (3.0, (a+2).y);
    ASSERT_DOUBLE_EQ (4.0, (a+3).z);
}
TEST_F (XyzTest, ScalarAddXYZ) {
    auto ans = 1 + c;
    ASSERT_DOUBLE_EQ (2.0, ans.x);
    ASSERT_DOUBLE_EQ (3.0, ans.y);
    ASSERT_DOUBLE_EQ (4.0, ans.z);
}
TEST_F (XyzTest, XYZSubtractXYZ) {
    auto ans = a - c;
    ASSERT_DOUBLE_EQ ( 0.0, ans.x);
    ASSERT_DOUBLE_EQ (-1.0, ans.y);
    ASSERT_DOUBLE_EQ (-2.0, ans.z);
}
TEST_F (XyzTest, XYZSubtractScalar) {
    auto ans = c - 1;
    ASSERT_DOUBLE_EQ (0.0, ans.x);
    ASSERT_DOUBLE_EQ (1.0, ans.y);
    ASSERT_DOUBLE_EQ (2.0, ans.z);
}
TEST_F (XyzTest, ScalarSubtractXYZ) {
    auto ans = 1 - c;
    ASSERT_DOUBLE_EQ ( 0.0, ans.x);
    ASSERT_DOUBLE_EQ (-1.0, ans.y);
    ASSERT_DOUBLE_EQ (-2.0, ans.z);
}
TEST_F (XyzTest, XYZDivideScalar) {
    auto ans = c / 4.0;
    ASSERT_DOUBLE_EQ ( 0.25, ans.x);
    ASSERT_DOUBLE_EQ ( 0.50, ans.y);
    ASSERT_DOUBLE_EQ (3/4.0, ans.z);
}

#endif
