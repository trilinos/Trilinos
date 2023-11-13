#include "stk_middle_mesh/discretization_1d.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(BasisFunctions, OtherDegree1)
{
  disc1d::impl::BasisFunctions basis(1);
  EXPECT_EQ(basis.get_num_nodes(), 2);
  EXPECT_EQ(basis.get_xi_min(), -1);
  EXPECT_EQ(basis.get_xi_max(), 1);

  std::vector<double> nodeXi = {-1, 1};
  for (size_t i = 0; i < nodeXi.size(); ++i)
    EXPECT_NEAR(basis.get_node_xi()[i], nodeXi[i], 1e-13);
}

TEST(BasisFunctions, OtherDegree2)
{
  disc1d::impl::BasisFunctions basis(2);
  EXPECT_EQ(basis.get_num_nodes(), 3);
  EXPECT_EQ(basis.get_xi_min(), -1);
  EXPECT_EQ(basis.get_xi_max(), 1);

  std::vector<double> nodeXi = {-1, 0.0, 1};
  for (size_t i = 0; i < nodeXi.size(); ++i)
    EXPECT_NEAR(basis.get_node_xi()[i], nodeXi[i], 1e-13);
}

TEST(BasisFunctions, ValuesDegree1)
{
  disc1d::impl::BasisFunctions basis(1);

  EXPECT_NEAR(basis.get_value(-1, 0), 1, 1e-13);
  EXPECT_NEAR(basis.get_value(1, 0), 0, 1e-13);

  EXPECT_NEAR(basis.get_value(-1, 1), 0, 1e-13);
  EXPECT_NEAR(basis.get_value(1, 1), 1, 1e-13);

  EXPECT_NEAR(basis.get_value(0, 0), 0.5, 1e-13);
  EXPECT_NEAR(basis.get_value(0, 0), 0.5, 1e-13);
}

TEST(BasisFunctions, DerivsDegree1)
{
  disc1d::impl::BasisFunctions basis(1);

  EXPECT_NEAR(basis.get_derivative(-1, 0), -0.5, 1e-13);
  EXPECT_NEAR(basis.get_derivative(1, 0), -0.5, 1e-13);

  EXPECT_NEAR(basis.get_derivative(-1, 1), 0.5, 1e-13);
  EXPECT_NEAR(basis.get_derivative(1, 1), 0.5, 1e-13);
}

TEST(BasisFunctions, ValuesDegree2)
{
  disc1d::impl::BasisFunctions basis(2);

  EXPECT_NEAR(basis.get_value(-1, 0), 1, 1e-13);
  EXPECT_NEAR(basis.get_value(0, 0), 0, 1e-13);
  EXPECT_NEAR(basis.get_value(1, 0), 0, 1e-13);

  EXPECT_NEAR(basis.get_value(-1, 1), 0, 1e-13);
  EXPECT_NEAR(basis.get_value(0, 1), 1, 1e-13);
  EXPECT_NEAR(basis.get_value(1, 1), 0, 1e-13);

  EXPECT_NEAR(basis.get_value(-1, 2), 0, 1e-13);
  EXPECT_NEAR(basis.get_value(0, 2), 0, 1e-13);
  EXPECT_NEAR(basis.get_value(1, 2), 1, 1e-13);
}

TEST(BasisFunctions, DerivsDegree2)
{
  disc1d::impl::BasisFunctions basis(2);

  EXPECT_NEAR(basis.get_derivative(-1, 0), -1.5, 1e-13);
  EXPECT_NEAR(basis.get_derivative(0, 0), -0.5, 1e-13);
  EXPECT_NEAR(basis.get_derivative(1, 0), 0.5, 1e-13);

  EXPECT_NEAR(basis.get_derivative(-1, 1), 2, 1e-13);
  EXPECT_NEAR(basis.get_derivative(0, 1), 0, 1e-13);
  EXPECT_NEAR(basis.get_derivative(1, 1), -2, 1e-13);

  EXPECT_NEAR(basis.get_derivative(-1, 2), -0.5, 1e-13);
  EXPECT_NEAR(basis.get_derivative(0, 2), 0.5, 1e-13);
  EXPECT_NEAR(basis.get_derivative(1, 2), 1.5, 1e-13);
}
} // namespace impl
} // namespace middle_mesh
} // namespace stk
