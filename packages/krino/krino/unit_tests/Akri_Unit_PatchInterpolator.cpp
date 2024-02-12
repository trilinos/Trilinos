#include <gtest/gtest.h>
#include <memory>

#include <Akri_PatchInterpolator.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

class PatchInterpolator2D : public ::testing::Test
{
protected:
    std::array<stk::math::Vector3d,2> copy_to_3d_array(const std::array<stk::math::Vector2d,2> & vec2d)
    {
      return std::array<stk::math::Vector3d,2>{ stk::math::Vector3d(vec2d[0][0], vec2d[0][1], 0.), stk::math::Vector3d(vec2d[1][0], vec2d[1][1], 0.) };
    }
    std::unique_ptr<CubicSplineInterpolator> build_patch_interpolator(const std::array<stk::math::Vector2d,2> & patchNodeCoords, const std::array<stk::math::Vector2d,2> & patchNodeNormalDirs)
    {
      const auto nodeCoords = copy_to_3d_array(patchNodeCoords);
      auto nodeNormals = copy_to_3d_array(patchNodeNormalDirs);
      nodeNormals[0].unitize();
      nodeNormals[1].unitize();

      return std::make_unique<CubicSplineInterpolator>(nodeCoords, nodeNormals);
    }

    void test_interpolated_location(const std::array<stk::math::Vector2d,2> & patchNodeCoords, const std::array<stk::math::Vector2d,2> & patchNodeNormalDirs, const double paramX, const stk::math::Vector2d & goldLoc, const double tol=0.0001)
    {
      auto interpolator = build_patch_interpolator(patchNodeCoords, patchNodeNormalDirs);
      const stk::math::Vector3d interpolatedLoc = interpolator->evaluate(copy_to_3d_array(patchNodeCoords), paramX);
      const stk::math::Vector3d gold3D{goldLoc[0], goldLoc[1], 0.};
      EXPECT_GT(tol, (gold3D-interpolatedLoc).length()) << "Mismatch between gold location " << goldLoc << " and interpolated location" << interpolatedLoc;
    }
};

TEST_F(PatchInterpolator2D, constant)
{
  const std::array<stk::math::Vector2d,2> patchNodeCoords{ stk::math::Vector2d{0.,0.}, stk::math::Vector2d{1.,0.} };

  const std::array<stk::math::Vector2d,2> patchNodeNormalDirs{ stk::math::Vector2d{0.,1.}, stk::math::Vector2d{0.,1.} };
  test_interpolated_location(patchNodeCoords, patchNodeNormalDirs, 0.5, stk::math::Vector2d(0.5,0.0));
}

TEST_F(PatchInterpolator2D, quadratic)
{
  const std::array<stk::math::Vector2d,2> patchNodeCoords{ stk::math::Vector2d{0.,0.}, stk::math::Vector2d{0.,1.} };
  const std::array<stk::math::Vector2d,2> patchNodeNormalDirs{ stk::math::Vector2d{1.,-1.}, stk::math::Vector2d{1.,1.} };
  for (int i=0; i<5; ++i)
  {
    const double pos = i/5.;
    const double func = 0.25 - (pos-0.5)*(pos-0.5);
    test_interpolated_location(patchNodeCoords, patchNodeNormalDirs, pos, stk::math::Vector2d(func, pos));
  }
}

TEST_F(PatchInterpolator2D, cubic)
{
  const std::array<stk::math::Vector2d,2> patchNodeCoords{ stk::math::Vector2d{0.,0.}, stk::math::Vector2d{0.,1.} };
  const std::array<stk::math::Vector2d,2> patchNodeNormalDirs{ stk::math::Vector2d{1.,1.}, stk::math::Vector2d{1.,1.} };
  for (int i=0; i<5; ++i)
  {
    const double pos = i/5.;
    const double func = -pos + 3.*pos*pos - 2.*pos*pos*pos;
    test_interpolated_location(patchNodeCoords, patchNodeNormalDirs, pos, stk::math::Vector2d(func, pos));
  }
}

TEST_F(PatchInterpolator2D, approximateCircle)
{
  const std::array<stk::math::Vector2d,2> patchNodeCoords{ stk::math::Vector2d{1.,0.}, stk::math::Vector2d{0.,1.} };
  const std::array<stk::math::Vector2d,2> patchNodeNormalDirs{ stk::math::Vector2d{1.,0.}, stk::math::Vector2d{0.,1.} };

  for (int i=0; i<5; ++i)
  {
    const double pos = i/5.;
    const stk::math::Vector2d locOnCircle = ((1.-pos)*patchNodeCoords[0] + pos*patchNodeCoords[1]).unit_vector();
    test_interpolated_location(patchNodeCoords, patchNodeNormalDirs, pos, locOnCircle, 0.15);
  }
}

}

