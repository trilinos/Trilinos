#include "stk_middle_mesh/bspline_patch.hpp"
#include "gtest/gtest.h"
#include <fstream>

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

std::vector<utils::Point> get_points_uniform(double xmin, double xmax, double ymin, double ymax, int nptsX, int nptsY)
{
  std::vector<utils::Point> xVals(nptsX * nptsY);
  double deltaX = (xmax - xmin) / (nptsX - 1);
  double deltaY = (ymax - ymin) / (nptsY - 1);

  for (int i = 0; i < nptsX; ++i)
    for (int j = 0; j < nptsY; ++j)
      xVals[i * nptsY + j] = utils::Point(xmin + deltaX * i, ymin + deltaY * j);

  return xVals;
}

template <typename T>
void get_f_values(std::vector<utils::Point>& values, T func)
{
  for (size_t i = 0; i < values.size(); ++i)
    values[i].z = func(values[i]);
}

/*
void writeFile(const std::string& fname, double xmin, double xmax, double ymin, double ymax,
               int npts_x, int npts_y, const std::vector<utils::Point>& values)
{
  std::ofstream xfile;
  xfile.open(fname + "_x.txt");
  xfile << std::setprecision(16);

  double delta_x = (xmax - xmin)/(npts_x - 1);
  for (int i=0; i < npts_x; ++i)
    xfile << xmin + delta_x * i << std::endl;

  xfile.close();

  std::ofstream yfile;
  yfile.open(fname + "_y.txt");
  yfile << std::setprecision(16);

  double delta_y = (ymax - ymin)/(npts_y - 1);
  for (int i=0; i < npts_y; ++i)
    yfile << ymin + delta_y * i << std::endl;

  yfile.close();


  std::ofstream valsfile;
  valsfile << std::setprecision(16);

  valsfile.open(fname + "_vals.txt");

  for (int j=0; j < npts_y; ++j)
  {
    for (int i=0; i < npts_x; ++i)
    {
      int idx = i * npts_y + j;
      valsfile << values[idx].z << " ";
    }

    valsfile << std::endl;
  }

  valsfile.close();
}
*/

template <typename T>
void test_bsplines(T f)
{
  auto vals = get_points_uniform(0, 2, 0, 2, 6, 6);
  get_f_values(vals, f);

  utils::impl::BSplinePatch patch(4);
  patch.construct_patch(vals);

  auto testVals = get_points_uniform(0, 2, 0, 2, 10, 10);
  for (size_t i = 0; i < testVals.size(); ++i)
  {
    auto testPt  = testVals[i];
    auto patchPt = patch.eval_point(testPt.x, testPt.y);
    EXPECT_NEAR(patchPt.z, f(testPt), 1e-13);
  }
}

} // namespace

// 1 dimensional polynomials
// 2 dimensional polynomials
// 2 dimensional polynomials specifying pt0 to be outside the range of the other data

TEST(BSplinePatch, Zero)
{
  auto f = [](const utils::Point& /*pt*/) { return 0; };
  test_bsplines(f);
}

TEST(BSplinePatch, Constant)
{
  auto f = [](const utils::Point& /*pt*/) { return 2; };
  test_bsplines(f);
}

TEST(BSplinePatch, Linearx)
{
  auto f = [](const utils::Point& pt) { return 2 * pt.x; };
  test_bsplines(f);
}

TEST(BSplinePatch, Quadraticx)
{
  auto f = [](const utils::Point& pt) { return 2 * pt.x * pt.x; };
  test_bsplines(f);
}

TEST(BSplinePatch, Cubicx)
{
  auto f = [](const utils::Point& pt) { return 2 * std::pow(pt.x, 3); };
  test_bsplines(f);
}

TEST(BSplinePatch, Lineary)
{
  auto f = [](const utils::Point& pt) { return 2 * pt.y; };
  test_bsplines(f);
}

TEST(BSplinePatch, Quadraticy)
{
  auto f = [](const utils::Point& pt) { return 2 * pt.y * pt.y; };
  test_bsplines(f);
}

TEST(BSplinePatch, Cubicy)
{
  auto f = [](const utils::Point& pt) { return 2 * std::pow(pt.y, 3); };
  test_bsplines(f);
}

TEST(BSplinePatch, Linearxy)
{
  auto f = [](const utils::Point& pt) { return 2 * pt.x * pt.y; };
  test_bsplines(f);
}

TEST(BSplinePatch, Quadraticxy)
{
  auto f = [](const utils::Point& pt) { return 2 * pt.x * pt.x * pt.y * pt.y; };
  test_bsplines(f);
}

TEST(BSplinePatch, Cubicxy)
{
  auto f = [](const utils::Point& pt) { return 2 * std::pow(pt.x, 3) * std::pow(pt.y, 3); };
  test_bsplines(f);
}

TEST(BSplinePatch, CubicxyUsingPt0)
{
  auto f    = [](const utils::Point& pt) { return 2 * std::pow(pt.x, 3) * std::pow(pt.y, 3); };
  auto vals = get_points_uniform(0, 2, 0, 2, 6, 6);
  get_f_values(vals, f);

  utils::Point pt0(3, 3);
  utils::impl::BSplinePatch patch(2);
  patch.construct_patch(vals, pt0);

  auto testVals = get_points_uniform(0, 3, 0, 3, 10, 10);
  for (size_t i = 0; i < testVals.size(); ++i)
  {
    auto testPt  = testVals[i];
    auto patchPt = patch.eval_point(testPt.x, testPt.y);
    EXPECT_NEAR(patchPt.z, f(testPt), 1e-11);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
