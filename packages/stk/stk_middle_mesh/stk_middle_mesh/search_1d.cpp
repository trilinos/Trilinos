#include "search_1d.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

mesh::MeshEntityPtr Search1D::find_containing_element(double x)
{
  auto it = std::lower_bound(m_vertCoords.begin(), m_vertCoords.end(), x);
  int idx = std::distance(m_vertCoords.begin(), it);
  // int elnum = idx == 0 ? 0 : idx - 1;
  int elnum = clamp(0, idx - 1, m_disc->get_mesh()->get_edges().size() - 1);
  if (elnum < 0 || size_t(elnum) >= m_disc->get_mesh()->get_edges().size())
    throw std::runtime_error("elnum is out of range");

  return m_disc->get_mesh()->get_edges()[elnum];
}

// given a point x that lies within el, compute its parametric coordinate
double Search1D::get_xi_coordinate(mesh::MeshEntityPtr el, double x)
{
  assert(get_type_dimension(el->get_type()) == 1);
  double xL = el->get_down(0)->get_point_orig(0).get_x();
  double xR = el->get_down(1)->get_point_orig(0).get_x();

  double xiZeroToOne = (x - xL) / (xR - xL);
  return 2 * xiZeroToOne - 1;
}

void Search1D::setup_for_search()
{
  int idx = 0;
  for (auto& vert : m_disc->get_mesh()->get_vertices())
    m_vertCoords[idx++] = vert->get_point_orig(0).get_x();
}

int Search1D::clamp(int min, int val, int max)
{
  int val1 = std::max(min, val);
  return std::min(val1, max);
}

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
