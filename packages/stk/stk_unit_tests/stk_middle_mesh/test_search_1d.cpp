#include "stk_middle_mesh/search_1d.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

struct TestData
{
    double xCoord;
    int elnum;
};
} // namespace

TEST(Search1D, Element)
{
  double xmin = 0;
  double xmax = 6;
  int numel   = 3;
  disc1d::impl::MeshSpec1D spec{numel, xmin, xmax};
  int basisDegree = 1, quadDegree = 3;

  auto disc = disc1d::impl::make_disc1_d(spec, basisDegree, quadDegree);

  disc1d::impl::Search1D search(disc);

  std::vector<TestData> testData = {{xmin, 0}, {1, 0},   {1.5, 0}, {2.4, 1}, {3.6, 1},
                                    {4.4, 2},  {5.6, 2}, {7.0, 2}, {-1, 0}};

  for (auto& data : testData)
  {
    std::cout << "\ntesting x coordinate " << data.xCoord << ", which should be in element " << data.elnum << std::endl;
    auto el = search.find_containing_element(data.xCoord);
    EXPECT_EQ(get_type_dimension(el->get_type()), 1);
    EXPECT_EQ(el->get_id(), data.elnum);

    double xiCoord = search.get_xi_coordinate(el, data.xCoord);

    double xiCoordZeroToOne = (xiCoord + 1) / 2;
    double xL               = el->get_down(0)->get_point_orig(0).get_x();
    double xR               = el->get_down(1)->get_point_orig(0).get_x();
    double xFromXi          = xL + (xR - xL) * xiCoordZeroToOne;
    EXPECT_NEAR(data.xCoord, xFromXi, 1e-13);
  }
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
