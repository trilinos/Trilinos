#include "simple_interpolator.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

void apply_interpolation_from_matrix(InterpolatorBase& interp, const std::vector<double>& srcFunc,
                                     std::vector<double>& destFunc)
{
  std::vector<std::vector<double>> coeffs;
  std::vector<std::vector<int>> indices;
  interp.get_interpolation_matrix(coeffs, indices);

  for (size_t iDest = 0; iDest < coeffs.size(); ++iDest)
  {
    auto& coeffsI  = coeffs[iDest];
    auto& indicesI = indices[iDest];

    destFunc[iDest] = 0;
    for (size_t jSrc = 0; jSrc < coeffsI.size(); ++jSrc)
      destFunc[iDest] += coeffsI[jSrc] * srcFunc[indicesI[jSrc]];
  }
}

void SimpleInterpolator::interpolate(const std::vector<double>& srcFunc, std::vector<double>& destFunc)
{
  assert(srcFunc.size() == size_t(m_discSrc->get_num_dofs()));
  assert(destFunc.size() == size_t(m_discDest->get_num_dofs()));

  // TODO: make Fields non-copyable
  auto& dofnums = *(m_discDest->get_dof_nums());
  auto& coords  = *(m_discDest->get_coords());
  for (auto vert : m_discDest->get_mesh()->get_vertices())
  {
    double valInterp              = interpolate_point(srcFunc, coords(vert, 0, 0));
    destFunc[dofnums(vert, 0, 0)] = valInterp;
  }

  if (m_discDest->get_field_shape().count[1] == 1)
  {
    for (auto edge : m_discDest->get_mesh()->get_edges())
    {
      double valInterp              = interpolate_point(srcFunc, coords(edge, 0, 0));
      destFunc[dofnums(edge, 0, 0)] = valInterp;
    }
  }
}

double SimpleInterpolator::interpolate_point(const std::vector<double>& srcFunc, double x)
{
  auto elSrc   = m_search.find_containing_element(x);
  double xiSrc = m_search.get_xi_coordinate(elSrc, x);

  std::vector<double> elVals;
  get_element_values(m_discSrc, elSrc, srcFunc, elVals);
  const auto& basis = m_discSrc->get_basis();

  double valInterp = 0;
  for (int i = 0; i < basis.get_num_nodes(); ++i)
    valInterp += basis.get_value(xiSrc, i) * elVals[i];

  return valInterp;
}

void SimpleInterpolator::get_interpolation_data(double x, std::vector<double>& coeffs, std::vector<int>& indices)
{
  auto elSrc   = m_search.find_containing_element(x);
  double xiSrc = m_search.get_xi_coordinate(elSrc, x);

  const auto& basis = m_discSrc->get_basis();

  coeffs.resize(basis.get_num_nodes());
  indices.resize(basis.get_num_nodes());

  auto& dofsSrc = *(m_discSrc->get_dof_nums());
  if (basis.get_degree() == 1)
  {
    indices[0] = dofsSrc(elSrc->get_down(0), 0, 0);
    indices[1] = dofsSrc(elSrc->get_down(1), 0, 0);
  } else if (basis.get_degree() == 2)
  {
    indices[0] = dofsSrc(elSrc->get_down(0), 0, 0);
    indices[1] = dofsSrc(elSrc, 0, 0);
    indices[2] = dofsSrc(elSrc->get_down(1), 0, 0);
  }

  for (int i = 0; i < basis.get_num_nodes(); ++i)
    coeffs[i] = basis.get_value(xiSrc, i);
}

void SimpleInterpolator::get_interpolation_matrix(std::vector<std::vector<double>>& coeffs,
                                                  std::vector<std::vector<int>>& indices)
{
  coeffs.resize(m_discDest->get_num_dofs());
  indices.resize(m_discDest->get_num_dofs());

  auto& dofnums = *(m_discDest->get_dof_nums());
  auto& coords  = *(m_discDest->get_coords());
  for (auto vert : m_discDest->get_mesh()->get_vertices())
    get_interpolation_data(coords(vert, 0, 0), coeffs[dofnums(vert, 0, 0)], indices[dofnums(vert, 0, 0)]);

  if (m_discDest->get_field_shape().count[1] == 1)
    for (auto edge : m_discDest->get_mesh()->get_edges())
      get_interpolation_data(coords(edge, 0, 0), coeffs[dofnums(edge, 0, 0)], indices[dofnums(edge, 0, 0)]);
}

} // namespace impl
} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
