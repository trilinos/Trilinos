#include "stencil_interpolator.hpp"
#include <algorithm>

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

void StencilInterpolator::interpolate(const std::vector<double>& srcFunc, std::vector<double>& destFunc)
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

void StencilInterpolator::get_interpolation_matrix(std::vector<std::vector<double>>& coeffs,
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

namespace {

mesh::MeshEntityPtr get_adjacent_element(mesh::MeshEntityPtr el, int vertIdx)
{
  auto vert = el->get_down(vertIdx);
  if (vert->count_up() == 1)
    return nullptr;
  else
  {
    auto el1 = vert->get_up(0);
    auto el2 = vert->get_up(1);
    return el1 == el ? el2 : el1;
  }
}

mesh::MeshEntityPtr get_left_adjacent_element(mesh::MeshEntityPtr el)
{
  return get_adjacent_element(el, 0);
}

mesh::MeshEntityPtr get_right_adjacent_element(mesh::MeshEntityPtr el)
{
  return get_adjacent_element(el, 1);
}

} // namespace

double StencilInterpolator::interpolate_point(const std::vector<double>& srcFunc, double x)
{
  // double xi_src = m_search.getXiCoordinate(el_src, x);
  auto elements = get_mesh_entities_in_stencil(x);
  auto coords   = get_node_coords(elements);
  auto dofs     = get_node_dofs(elements);
  auto coeffs   = compute_coefficients(coords, x);

  int npts         = coeffs.size();
  double valInterp = 0;
  for (int i = 0; i < npts; ++i)
    valInterp += coeffs[i] * srcFunc[dofs[i]];

  return valInterp;
}

void StencilInterpolator::get_interpolation_data(double x, std::vector<double>& coeffs, std::vector<int>& indices)
{
  auto elements = get_mesh_entities_in_stencil(x);
  auto coords   = get_node_coords(elements);
  indices       = get_node_dofs(elements);
  coeffs        = compute_coefficients(coords, x);
}

std::vector<mesh::MeshEntityPtr> StencilInterpolator::get_mesh_entities_in_stencil(double x)
{
  auto elSrc = m_search.find_containing_element(x);

  int nptsInBasis = m_discSrc->get_basis().get_num_nodes();
  int npts        = nptsInBasis;

  mesh::MeshEntityPtr leftEl = elSrc, rightEl = elSrc;
  std::vector<mesh::MeshEntityPtr> elements = {elSrc};
  while (npts < m_stencilsize)
  {
    auto leftAdjacentEl = get_left_adjacent_element(leftEl);
    if (leftAdjacentEl)
    {
      elements.push_back(leftAdjacentEl);
      leftEl = leftAdjacentEl;
      npts += nptsInBasis - 1;
    }

    if (npts >= m_stencilsize)
      break;

    auto rightAdjacentEl = get_right_adjacent_element(rightEl);
    if (rightAdjacentEl)
    {
      elements.push_back(rightAdjacentEl);
      rightEl = rightAdjacentEl;
      npts += nptsInBasis - 1;
    }
  }

  auto comp = [](const mesh::MeshEntityPtr& e1, const mesh::MeshEntityPtr& e2) { return e1->get_id() < e2->get_id(); };
  std::sort(elements.begin(), elements.end(), comp);

  return elements;
}

std::vector<double> StencilInterpolator::get_node_coords(const std::vector<mesh::MeshEntityPtr>& elements)
{
  std::vector<double> nodeCoords;
  auto& coords = *(m_discSrc->get_coords());

  nodeCoords.push_back(coords(elements[0]->get_down(0), 0, 0));
  for (auto el : elements)
  {
    if (m_discSrc->get_field_shape().count[1] == 1)
      nodeCoords.push_back(coords(el, 0, 0));

    nodeCoords.push_back(coords(el->get_down(1), 0, 0));
  }

  return nodeCoords;
}

std::vector<int> StencilInterpolator::get_node_dofs(const std::vector<mesh::MeshEntityPtr>& elements)
{
  std::vector<int> nodeDofs;
  auto& dofs = *(m_discSrc->get_dof_nums());

  nodeDofs.push_back(dofs(elements[0]->get_down(0), 0, 0));
  for (auto el : elements)
  {
    if (m_discSrc->get_field_shape().count[1] == 1)
      nodeDofs.push_back(dofs(el, 0, 0));

    nodeDofs.push_back(dofs(el->get_down(1), 0, 0));
  }

  return nodeDofs;
}

std::vector<double> StencilInterpolator::compute_coefficients(const std::vector<double>& nodeCoords, double x)
{
  int npts = nodeCoords.size();
  utils::impl::Matrix<double> a(npts, npts);
  std::vector<double> b(npts);

  double xmin = nodeCoords[0];
  double xmax = nodeCoords[1];

  double xXi = 2 * (x - xmin) / (xmax - xmin) - 1;

  LegendrePolynomials legendre;
  for (int polyDegree = 0; polyDegree < npts; ++polyDegree)
  {
    b[polyDegree] = legendre.compute(polyDegree, xXi);

    for (int pt = 0; pt < npts; ++pt)
    {
      double xPt        = nodeCoords[pt];
      double xiPt       = 2 * (xPt - xmin) / (xmax - xmin) - 1;
      a(polyDegree, pt) = legendre.compute(polyDegree, xiPt);
    }
  }

  std::vector<int> ipiv(npts);
  solve_linear_system(a, ipiv.data(), b.data());

  return b;
}

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
