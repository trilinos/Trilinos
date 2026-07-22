#include "stk_middle_mesh/discretization_1d.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

TEST(Discretization1D, MeshEntityCounts)
{
  double xmin = 0;
  double xmax = 4;
  int numel   = 10;
  disc1d::impl::MeshSpec1D spec{numel, xmin, xmax};
  int basisDegree = 1, quadDegree = 3;

  auto disc = disc1d::impl::make_disc1_d(spec, basisDegree, quadDegree);

  EXPECT_EQ(count_valid(disc->get_mesh()->get_vertices()), numel + 1);
  EXPECT_EQ(count_valid(disc->get_mesh()->get_edges()), numel);

  for (auto& edge : disc->get_mesh()->get_edges())
    EXPECT_EQ(edge->count_down(), 2);
}

TEST(Discretization1D, DofNums)
{
  double xmin = 0;
  double xmax = 4;
  int numel   = 10;
  disc1d::impl::MeshSpec1D spec{numel, xmin, xmax};
  int basisDegree = 1, quadDegree = 3;

  auto disc = disc1d::impl::make_disc1_d(spec, basisDegree, quadDegree);

  EXPECT_EQ(disc->get_num_dofs(), numel + 1);
  std::vector<int> usedDofs(disc->get_num_dofs(), 0);
  auto& dofs = *(disc->get_dof_nums());
  for (auto& vert : disc->get_mesh()->get_vertices())
  {
    int dof = dofs(vert, 0, 0);
    EXPECT_GE(dof, 0);
    EXPECT_LT(dof, disc->get_num_dofs());
    usedDofs[dof] += 1;
  }

  for (auto dofCount : usedDofs)
    EXPECT_EQ(dofCount, 1);
}

TEST(Discretization1D, IntegrateConstant)
{
  double xmin = 0;
  double xmax = 4;
  disc1d::impl::MeshSpec1D spec{10, xmin, xmax};

  for (int basisDegree = 1; basisDegree <= 2; ++basisDegree)
  {
    int quadDegree = 3;
    auto disc      = disc1d::impl::make_disc1_d(spec, basisDegree, quadDegree);
    auto f         = [](double /*x*/) { return 1; };
    auto func      = disc1d::impl::make_grid_function(disc, f);
    auto val       = disc1d::impl::integrate_function(disc, func);

    EXPECT_NEAR(val, xmax - xmin, 1e-13);
  }
}

TEST(Discretization1D, IntegrateLinear)
{
  double xmin = 0;
  double xmax = 4;
  disc1d::impl::MeshSpec1D spec{10, xmin, xmax};

  for (int basisDegree = 1; basisDegree <= 2; ++basisDegree)
  {
    int quadDegree = 3;
    auto disc      = disc1d::impl::make_disc1_d(spec, basisDegree, quadDegree);
    auto f         = [](double x) { return x; };
    auto func      = disc1d::impl::make_grid_function(disc, f);
    auto val       = disc1d::impl::integrate_function(disc, func);

    EXPECT_NEAR(val, 0.5 * std::pow(xmax, 2) - 0.5 * std::pow(xmin, 2), 1e-13);
  }
}

TEST(Discretization1D, SummedMassMatrix)
{
  double xmin = 0;
  double xmax = 4;
  disc1d::impl::MeshSpec1D spec{10, xmin, xmax};

  for (int basisDegree = 1; basisDegree <= 2; ++basisDegree)
  {
    int quadDegree = 3;
    auto disc      = disc1d::impl::make_disc1_d(spec, basisDegree, quadDegree);
    auto f         = [](double x) { return x; };
    auto func      = disc1d::impl::make_grid_function(disc, f);

    auto summedMassMatrix = disc1d::impl::compute_summed_mass_matrix(disc);
    // EXPECT_EQ(summed_mass_matrix.size(), disc->get_num_dofs());
    EXPECT_EQ((int)func.size(), disc->get_num_dofs());

    double val = 0.0;
    for (int i = 0; i < disc->get_num_dofs(); ++i)
      val += summedMassMatrix[i] * func[i];

    EXPECT_NEAR(val, 0.5 * std::pow(xmax, 2) - 0.5 * std::pow(xmin, 2), 1e-13);
  }
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
