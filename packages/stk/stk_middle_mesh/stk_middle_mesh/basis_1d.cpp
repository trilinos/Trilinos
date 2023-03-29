#include "basis_1d.hpp"
#include <iostream>
#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

BasisFunctions::BasisFunctions(int degree)
  : m_degree(degree)
{
  if (degree == 1)
    m_nodeXi = {-1, 1};
  else if (degree == 2)
    m_nodeXi = {-1, 0, 1};
  else
    throw std::runtime_error("only degree 1 and 2 basis functions are supported");
}

double BasisFunctions::get_value(double xi, int node) const
{
  check_xi(xi);
  if (m_degree == 1)
    return get_value_degree1(xi, node);
  else if (m_degree == 2)
    return get_value_degree2(xi, node);

  throw std::runtime_error("unsupported degree");
}

double BasisFunctions::get_derivative(double xi, int node) const
{
  check_xi(xi);
  if (m_degree == 1)
    return get_deriv_degree1(xi, node);
  else if (m_degree == 2)
    return get_deriv_degree2(xi, node);

  throw std::runtime_error("unsupported degree");
}

double BasisFunctions::get_value_degree1(double xi, int node) const
{
  if (node == 0)
    return (1 - xi) / 2;
  else if (node == 1)
    return (1 + xi) / 2;

  throw std::runtime_error("unsupported node");
}

double BasisFunctions::get_deriv_degree1(double /*xi*/, int node) const
{
  if (node == 0)
    return -1.0 / 2;
  else if (node == 1)
    return 1.0 / 2;

  throw std::runtime_error("unsupported node");
}

double BasisFunctions::get_value_degree2(double xi, int node) const
{
  if (node == 0)
    return xi * (xi - 1) / 2;
  else if (node == 1)
    return 1 - xi * xi;
  else if (node == 2)
    return xi * (xi + 1) / 2;

  throw std::runtime_error("unsupported node");
}

double BasisFunctions::get_deriv_degree2(double xi, int node) const
{
  if (node == 0)
    return (2 * xi - 1) / 2;
  else if (node == 1)
    return -2 * xi;
  else if (node == 2)
    return (2 * xi + 1) / 2;

  throw std::runtime_error("unsupported node");
}

void BasisFunctions::check_xi(double xi) const
{
  if (xi < get_xi_min() - 1e-8 || xi > get_xi_max() + 1e-8)
    throw std::runtime_error("xi out of range");
}

BasisFunctionValues::BasisFunctionValues(BasisFunctions& basis, const std::vector<double>& xiPoints)
  : m_vals(basis.get_num_nodes() * xiPoints.size())
  , m_derivs(basis.get_num_nodes() * xiPoints.size())
  , m_numNodes(basis.get_num_nodes())
  , m_numPoints(xiPoints.size())
{
  for (int node = 0; node < basis.get_num_nodes(); ++node)
    for (size_t i = 0; i < xiPoints.size(); ++i)
    {
      m_vals[get_idx(node, i)]   = basis.get_value(xiPoints[i], node);
      m_derivs[get_idx(node, i)] = basis.get_derivative(xiPoints[i], node);
    }
}

} // namespace impl
} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
