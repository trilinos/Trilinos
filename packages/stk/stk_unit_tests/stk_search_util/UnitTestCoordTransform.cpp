#include "gtest/gtest.h"
#include "stk_search_util/CoordTransform.hpp"

namespace {

double pi()
{
  return 2*std::atan2(1.0, 0.0);
}
}

TEST(CoordTransformIdentity, TwoD)
{
  stk::search::CoordTransformIdentity transform;
  std::vector<double> coords{1, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_EQ(coords[0], 1.0);
  EXPECT_EQ(coords[1], 2.0);
}

TEST(CoordTransformIdentity, ThreeD)
{
  stk::search::CoordTransformIdentity transform;
  std::vector<double> coords{1, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[0], 1.0);
  EXPECT_EQ(coords[1], 2.0);
  EXPECT_EQ(coords[2], 3.0);
}

TEST(CoordTransformAddZ, AddsZDefaultExpression)
{
  stk::search::CoordTransformAddZ transform;
  std::vector<double> coords{1, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[2], 0.0);
}

TEST(CoordTransformAddZ, AddsZConstantExpression)
{
  stk::search::CoordTransformAddZ transform("z=5.4321");
  std::vector<double> coords{1, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[2], 5.4321);
}

TEST(CoordTransformAddZ, AddsZConstantExpression2)
{
  stk::search::CoordTransformAddZ transform("3.5");
  std::vector<double> coords{1, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[2], 3.5);
}

TEST(CoordTransformAddZ, AddsZVariableExpression)
{
  stk::search::CoordTransformAddZ transform("z=x+y");
  std::vector<double> coords{1, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[2], 3);
}

TEST(CoordTransformAddZ, AddsZVariableExpression2)
{
  stk::search::CoordTransformAddZ transform("z=max(x,y)+4");
  std::vector<double> coords{1, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[2], 6);
}

TEST(CoordTransformAddZ, AddsZInvalidExpression)
{
  stk::search::CoordTransformAddZ transform("z=w");
  std::vector<double> coords{1, 2};
  EXPECT_ANY_THROW(transform.transform(stk::search::spmd::EntityKeyPair{}, coords));
}

TEST(CoordTransformAddZ, AddsZInvalidExpression2)
{
  stk::search::CoordTransformAddZ transform("z=w+x");
  std::vector<double> coords{1, 2};
  EXPECT_ANY_THROW(transform.transform(stk::search::spmd::EntityKeyPair{}, coords));
}

TEST(CoordTransformAddZ, AddsZInvalidExpression3)
{
  stk::search::CoordTransformAddZ transform("z=w+x*y");
  std::vector<double> coords{1, 2};
  EXPECT_ANY_THROW(transform.transform(stk::search::spmd::EntityKeyPair{}, coords));
}

TEST(CoordTransformAddZ, AddsZInvalidExpression4)
{
  stk::search::CoordTransformAddZ transform("w=x");
  std::vector<double> coords{1, 2};
  EXPECT_ANY_THROW(transform.transform(stk::search::spmd::EntityKeyPair{}, coords));
}

TEST(CoordTransformAddZ, AddsZInvalidExpression5)
{
  stk::search::CoordTransformAddZ transform("w=x;z=y");
  std::vector<double> coords{1, 2};
  EXPECT_ANY_THROW(transform.transform(stk::search::spmd::EntityKeyPair{}, coords));
}

TEST(CoordTransformRemoveZ, RemoveZ)
{
  stk::search::CoordTransformRemoveZ transform;
  std::vector<double> coords{1, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_EQ(coords[0], 1.0);
  EXPECT_EQ(coords[1], 2.0);
}

TEST(CoordTransformTranslate, TwoD)
{
  std::vector<double> coords{1, 2};
  std::array<double, 3> shift = {4, 5, 0};

  stk::search::CoordTransformTranslate transform(shift);
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_EQ(coords[0], 5.0);
  EXPECT_EQ(coords[1], 7.0);
}

TEST(CoordTransformTranslate, ThreeD)
{
  std::vector<double> coords{1, 2, 3};
  std::array<double, 3> shift = {4, 5, 6};

  stk::search::CoordTransformTranslate transform(shift);
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[0], 5.0);
  EXPECT_EQ(coords[1], 7.0);
  EXPECT_EQ(coords[2], 9.0);
}

TEST(CoordTransformScale, TwoD)
{
  std::vector<double> coords{1, 2};
  std::array<double, 3> scale = {4, 5, 0};

  stk::search::CoordTransformScale transform(scale);
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_EQ(coords[0], 4.0);
  EXPECT_EQ(coords[1], 10.0);
}

TEST(CoordTransformScale, ThreeD)
{
  std::vector<double> coords{1, 2, 3};
  std::array<double, 3> scale = {4, 5, 6};

  stk::search::CoordTransformScale transform(scale);
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[0], 4.0);
  EXPECT_EQ(coords[1], 10.0);
  EXPECT_EQ(coords[2], 18.0);
}

TEST(CoordTransformPermute, Shift3D)
{
  std::array<unsigned, 3> dests = {1, 2, 0};
  std::vector<double> coords = {1, 2, 3};

  stk::search::CoordTransformPermute transform(dests);
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_EQ(coords[0], 3);
  EXPECT_EQ(coords[1], 1);
  EXPECT_EQ(coords[2], 2);
}

TEST(CoordTransformPermute, Shift2D)
{
  std::array<unsigned, 3> dests = {1, 0, static_cast<unsigned>(-1)};
  std::vector<double> coords = {1, 2};

  stk::search::CoordTransformPermute transform(dests);
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_EQ(coords[0], 2);
  EXPECT_EQ(coords[1], 1);
}

TEST(CoordTransformPermute, ThrowOnDupliateDest)
{
  std::array<unsigned, 3> dests = {1, 1, 2};
  EXPECT_ANY_THROW(stk::search::CoordTransformPermute transform(dests));
}

TEST(CoordTransformPermute, ThrowOnOutOfRangeIndex)
{
  std::array<unsigned, 3> dests = {1, 2, 3};
  EXPECT_ANY_THROW(stk::search::CoordTransformPermute transform(dests));
}

TEST(CoordTransformPermute, ThrowOnOutOfRangeIndex2ndIndex)
{
  std::array<unsigned, 3> dests = {1, 3, 2};
  EXPECT_ANY_THROW(stk::search::CoordTransformPermute transform(dests));
}

TEST(CoordTransformRotateZ, ZtoX)
{
  std::array<double, 3> new_z{1, 0, 0};
  std::vector<double> coords{1, 0, 0};

  stk::search::CoordTransformRotateZ transform(new_z);
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);

  EXPECT_NEAR(coords[0], 0, 1e-13);
  EXPECT_NEAR(coords[1], 0, 1e-13);
  EXPECT_NEAR(coords[2], 1, 1e-13);

  stk::search::CoordTransformRotateZ inverse_transform(new_z, false);
  inverse_transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_NEAR(coords[0], 1, 1e-13);
  EXPECT_NEAR(coords[1], 0, 1e-13);
  EXPECT_NEAR(coords[2], 0, 1e-13);
}

TEST(CoordTransformCylindrical, TwoD)
{
  std::vector<double> coords{std::sqrt(3)/2, 0.5};

  stk::search::CoordTransformCylindricalCoordinates transform;
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_DOUBLE_EQ(coords[0], 1);
  EXPECT_DOUBLE_EQ(coords[1], pi()/6);
}

TEST(CoordTransformCylindrical, ThreeD)
{
  std::vector<double> coords{std::sqrt(3)/2, 0.5, 10};

  stk::search::CoordTransformCylindricalCoordinates transform;
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_DOUBLE_EQ(coords[0], 1);
  EXPECT_DOUBLE_EQ(coords[1], pi()/6);
  EXPECT_EQ(coords[2], 10.0);
}

TEST(CoordTransformCylindrical, Quadrant2)
{
  std::vector<double> coords{-std::sqrt(3)/2, 0.5};

  stk::search::CoordTransformCylindricalCoordinates transform;
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_DOUBLE_EQ(coords[0], 1);
  EXPECT_DOUBLE_EQ(coords[1], 5*pi()/6);
}

TEST(CoordTransformCylindrical, Quadrant3)
{
  std::vector<double> coords{-std::sqrt(3)/2, -0.5};

  stk::search::CoordTransformCylindricalCoordinates transform;
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_DOUBLE_EQ(coords[0], 1);
  EXPECT_DOUBLE_EQ(coords[1], pi()/6 - pi());
}

TEST(CoordTransformCylindrical, Quadrant4)
{
  std::vector<double> coords{std::sqrt(3)/2, -0.5};

  stk::search::CoordTransformCylindricalCoordinates transform;
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_DOUBLE_EQ(coords[0], 1);
  EXPECT_DOUBLE_EQ(coords[1], -pi()/6);
}

TEST(CoordTransformCartesian, TwoD)
{
  std::vector<double> coords{2, pi()/6};

  stk::search::CoordTransformCartesianCoordinates transform;
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 2u);
  EXPECT_DOUBLE_EQ(coords[0], std::sqrt(3));
  EXPECT_DOUBLE_EQ(coords[1], 1.0);
}

TEST(CoordTransformCartesian, ThreeD)
{
  std::vector<double> coords{2, pi()/6, 10};

  stk::search::CoordTransformCartesianCoordinates transform;
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_DOUBLE_EQ(coords[0], std::sqrt(3));
  EXPECT_DOUBLE_EQ(coords[1], 1.0);
  EXPECT_EQ(coords[2], 10.0);
}

TEST(CoordTransformPeriodicSymmetry, UnitRange)
{
  int dim = 0;
  double xmin = 0;
  double xmax = 1;
  stk::search::CoordTransformPeriodicSymmetry transform(dim, xmin, xmax);

  std::vector<double> coords{0.5, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_DOUBLE_EQ(coords[0], 0.5);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);

  coords = {1.2, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_DOUBLE_EQ(coords[0], 0.8);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);

  coords = {2.2, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_NEAR(coords[0], 0.2, 1e-13);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);

  coords = {-0.2, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_NEAR(coords[0], 0.2, 1e-13);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);

  coords = {-1.2, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_NEAR(coords[0], 0.8, 1e-13);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);

  coords = {-2.2, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_NEAR(coords[0], 0.2, 1e-13);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);
}

TEST(CoordTransformPeriodicSymmetry, NonUnitRange)
{
  int dim = 0;
  double xmin = 1;
  double xmax = 3;
  stk::search::CoordTransformPeriodicSymmetry transform(dim, xmin, xmax);

  std::vector<double> coords{1.5, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_DOUBLE_EQ(coords[0], 1 + 0.25*2);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);

  coords = {3.5, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_DOUBLE_EQ(coords[0], 3 - 0.25*2);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);

  coords = {5.5, 2, 3};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);
  EXPECT_EQ(coords.size(), 3u);
  EXPECT_DOUBLE_EQ(coords[0], 1 + 0.25*2);
  EXPECT_DOUBLE_EQ(coords[1], 2);
  EXPECT_DOUBLE_EQ(coords[2], 3);
}

TEST(CoordTransformAxisymmetry, DefaultZAxis)
{
  stk::search::CoordTransformAxisymmetric2D transform;

  std::vector<double> coords{std::sqrt(3.0)/2, 0.5, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);

  EXPECT_DOUBLE_EQ(coords[0], 1);  // r
  EXPECT_DOUBLE_EQ(coords[1], 2);  // z

  coords = {-std::sqrt(3.0)/2, -0.5, 2};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);

  EXPECT_DOUBLE_EQ(coords[0], 1);  // r
  EXPECT_DOUBLE_EQ(coords[1], 2);  // z
}

TEST(CoordTransformAxisymmetry, XAxisWithTranslation)
{
  std::array<double, 3> axis = {1, 0, 0};
  std::array<double, 3> pt = {0, 1, 0};
  stk::search::CoordTransformAxisymmetric2D transform(axis, pt);

  std::vector<double> coords{2, std::sqrt(3.0)/2 - 1, 0.5};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);

  EXPECT_DOUBLE_EQ(coords[0], 1);  // r
  EXPECT_DOUBLE_EQ(coords[1], 2);  // z

  coords = {2, -std::sqrt(3.0)/2 - 1, -0.5};
  transform.transform(stk::search::spmd::EntityKeyPair{}, coords);

  EXPECT_DOUBLE_EQ(coords[0], 1);  // r
  EXPECT_DOUBLE_EQ(coords[1], 2);  // z
}