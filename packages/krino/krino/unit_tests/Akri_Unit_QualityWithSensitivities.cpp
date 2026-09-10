// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_QualityMetricWithSensitivities.hpp>
#include <Akri_UnitTestUtils.hpp>
#include <gtest/gtest.h>

#include <stk_math/StkVector.hpp>

namespace krino {

void test_tet_mean_ratio_and_sensitivities(const std::array<stk::math::Vector3d,4> & tetCoords)
{
  const auto [qualUsingFAD,sensUsingFAD] = MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivities_using_FAD(tetCoords);
  const auto [qualAnalytic,sensAnalytic] = MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivities(tetCoords);

  const double tol = 1.e-6;
  EXPECT_NEAR(qualUsingFAD, qualAnalytic, tol);
  for (int i=0; i<4; ++i)
    expect_near(sensUsingFAD[i], sensAnalytic[i], tol);
}

TEST(MeanRatioQualityWithSensitivities, testValueAndSensitivities3d)
{
  const std::array<stk::math::Vector3d,4> rightTet
  {{
      { 0.0, 0.0, 0.0 },
      { 1.0, 0.0, 0.0 },
      { 0.0, 1.0, 0.0 },
      { 0.0, 0.0, 1.0 },
  }};

  test_tet_mean_ratio_and_sensitivities(rightTet);
}

void test_tri_mean_ratio_and_sensitivities(const std::array<stk::math::Vector2d,3> & triCoords)
{
  const auto [qualUsingFAD,sensUsingFAD] = MeanRatioQualityMetricWithSensitivities::tri2d_quality_and_sensitivities_using_FAD(triCoords[0].data(), triCoords[1].data(), triCoords[2].data());
  const auto [qualAnalytic,sensAnalytic] = MeanRatioQualityMetricWithSensitivities::tri2d_quality_and_sensitivities(triCoords[0].data(), triCoords[1].data(), triCoords[2].data());

  const double tol = 1.e-6;
  EXPECT_NEAR(qualUsingFAD, qualAnalytic, tol);
  for (int i=0; i<6; ++i)
    EXPECT_NEAR(sensUsingFAD[i], sensAnalytic[i], tol);
}

TEST(MeanRatioQualityWithSensitivities, testValueAndSensitivities2d)
{
  const std::array<stk::math::Vector2d,3> rightTri
  {{
      { 0.0, 0.0, },
      { 1.0, 0.0 },
      { 0.0, 1.0 }
  }};

  test_tri_mean_ratio_and_sensitivities(rightTri);
}

}
