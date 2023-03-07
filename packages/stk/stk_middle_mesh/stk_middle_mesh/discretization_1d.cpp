#include "discretization_1d.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

mesh::FieldShape get_field_shape(int degree)
{
  if (degree == 1)
    return mesh::FieldShape(1, 0, 0);
  else if (degree == 2)
    return mesh::FieldShape(1, 1, 0);
  else
    throw std::runtime_error("only degree 1 and 2 fieldshapes supported");
}

void Discretization1D::set_dof_nums()
{
  const int notAssigned = -1;
  m_dofnums->set(notAssigned);

  auto& fieldRef = *m_dofnums;
  int dofnum     = 0;
  for (auto& edge : m_mesh->get_edges())
  {
    mesh::MeshEntityPtr vertL = edge->get_down(0);
    mesh::MeshEntityPtr vertR = edge->get_down(1);
    if (fieldRef(vertL, 0, 0) == notAssigned)
      fieldRef(vertL, 0, 0) = dofnum++;

    if (m_fieldshape.count[1] == 1)
      fieldRef(edge, 0, 0) = dofnum++;

    if (fieldRef(vertR, 0, 0) == notAssigned)
      fieldRef(vertR, 0, 0) = dofnum++;
  }

  int numExpectedDofs = m_fieldshape.count[0] * count_valid(m_mesh->get_vertices()) +
                        m_fieldshape.count[1] * count_valid(m_mesh->get_edges());

  if (dofnum != numExpectedDofs)
    throw std::runtime_error("number of dofs is incorrect");

  m_numDofs = dofnum;
}

void Discretization1D::compute_coordinates()
{
  auto& coords = *m_coords;
  for (auto& vert : get_mesh()->get_vertices())
    coords(vert, 0, 0) = vert->get_point_orig(0).get_x();

  if (m_fieldshape.count[1] == 1)
  {
    for (auto& edge : get_mesh()->get_edges())
    {
      double xL          = edge->get_down(0)->get_point_orig(0).get_x();
      double xR          = edge->get_down(1)->get_point_orig(0).get_x();
      coords(edge, 0, 0) = (xL + xR) / 2;
    }
  }
}

void Discretization1D::compute_metrics()
{
  // for linear elements, dx/dxi = delta_x / delta_xi
  auto& dxdxiRef  = *(m_dxdxi);
  auto& coordsRef = *(m_coords);
  for (auto edge : m_mesh->get_edges())
  {
    mesh::MeshEntityPtr vertL = edge->get_down(0);
    mesh::MeshEntityPtr vertR = edge->get_down(1);
    double deltaX             = coordsRef(vertR, 0, 0) - coordsRef(vertL, 0, 0);
    double deltaXi            = m_basis.get_xi_max() - m_basis.get_xi_min();

    for (int i = 0; i < m_quad.get_num_points(); ++i)
      dxdxiRef(edge, 0, i) = deltaX / deltaXi;
  }
}

Discretization1DPtr make_disc1_d(MeshSpec1D meshspec, int basisDegree, int quadDegree)
{
  return std::shared_ptr<Discretization1D>(new Discretization1D(meshspec, basisDegree, quadDegree));
}

std::shared_ptr<mesh::Mesh> make1_d_mesh(MeshSpec1D meshspec)
{
  auto mesh = mesh::make_empty_mesh();

  double deltaX             = (meshspec.xmax - meshspec.xmin) / meshspec.numel;
  double xI                 = meshspec.xmin;
  mesh::MeshEntityPtr vertL = mesh->create_vertex(xI, 0, 0);
  mesh::MeshEntityPtr vertR = nullptr;

  for (int i = 0; i < meshspec.numel; ++i)
  {
    xI += deltaX;
    vertR = mesh->create_vertex(xI, 0, 0);
    mesh->create_edge(vertL, vertR);
    vertL = vertR;
  }

  return mesh;
}

void get_element_values(Discretization1DPtr disc, mesh::MeshEntityPtr el, const std::vector<double>& valsMesh,
                        std::vector<double>& elVals)
{
  assert(get_type_dimension(el->get_type()) == 1);
  assert(valsMesh.size() == size_t(disc->get_num_dofs()));

  auto& dofsRef = *(disc->get_dof_nums());

  auto vert1 = el->get_down(0);
  auto vert2 = el->get_down(1);

  if (disc->get_basis_degree() == 1)
  {
    elVals.resize(2);
    elVals[0] = valsMesh[dofsRef(vert1, 0, 0)];
    elVals[1] = valsMesh[dofsRef(vert2, 0, 0)];
  } else if (disc->get_basis_degree() == 2)
  {
    elVals.resize(3);
    elVals[0] = valsMesh[dofsRef(vert1, 0, 0)];
    elVals[1] = valsMesh[dofsRef(el, 0, 0)];
    elVals[2] = valsMesh[dofsRef(vert2, 0, 0)];
  } else
    throw std::runtime_error("unsupported degree");
}

void add_element_values(Discretization1DPtr disc, mesh::MeshEntityPtr el, const std::vector<double>& elVals,
                        std::vector<double>& valsMesh)
{
  assert(get_type_dimension(el->get_type()) == 1);
  assert(valsMesh.size() == size_t(disc->get_num_dofs()));

  auto& dofsRef = *(disc->get_dof_nums());

  auto vert1 = el->get_down(0);
  auto vert2 = el->get_down(1);

  if (disc->get_basis_degree() == 1)
  {
    valsMesh[dofsRef(vert1, 0, 0)] += elVals[0];
    valsMesh[dofsRef(vert2, 0, 0)] += elVals[1];
  } else if (disc->get_basis_degree() == 2)
  {
    valsMesh[dofsRef(vert1, 0, 0)] += elVals[0];
    valsMesh[dofsRef(el, 0, 0)] += elVals[1];
    valsMesh[dofsRef(vert2, 0, 0)] += elVals[2];
  } else
    throw std::runtime_error("unsupported degree");
}

double integrate_function(Discretization1DPtr disc, const std::vector<double>& uVals)
{
  int numEdges    = count_valid(disc->get_mesh()->get_edges());
  auto& quad      = disc->get_quadrature();
  auto& basisVals = disc->get_quad_basis_values();
  auto& dxdxi     = *(disc->getd_xdxi());

  double val = 0;
  std::vector<double> elVals;
  for (int i = 0; i < numEdges; ++i)
  {
    auto edge = disc->get_mesh()->get_edges()[i];
    get_element_values(disc, edge, uVals, elVals);

    for (int k = 0; k < quad.get_num_points(); ++k)
    {
      double uValInterp = 0;
      for (int j = 0; j < basisVals.get_num_nodes(); ++j)
        uValInterp += basisVals.get_value(j, k) * elVals[j]; // TODO: reorder basis_vals to cache efficiency

      val += quad.get_weights()[k] * uValInterp * dxdxi(edge, 0, k);
    }
  }

  return val;
}

// TODO: untested
/*
void computeMassMatrix(Discretization1DPtr disc)
{
  int num_edges = count_valid(disc->get_mesh()->get_edges());
  auto& quad = disc->get_quadrature();
  auto& basis_vals = disc->get_quad_basis_values();
  auto& dxdxi = *(disc->getd_xdxi());

  utils::impl::Matrix<double> mass_el(basis_vals.get_num_nodes(), basis_vals.get_num_nodes());
  for (int i=0; i < num_edges; ++i)
  {
    auto edge = disc->get_mesh()->get_edges()[i];
    mass_el.fill(0);

    for (int k=0; k < quad.get_num_points(); ++k)
    {
      double u_val_interp = 0;
      for (int j1=0; j1 < basis_vals.get_num_nodes(); ++j1)
        for (int j2=0; j2 < basis_vals.get_num_nodes(); ++j2)
          mass_el(j1, j2) += basis_vals.get_value(j1, k) * quad.getWeights()[k] * dxdxi(edge, 0, k) *
basis_vals.get_value(j2, k);
    }
  }
}
*/

std::vector<double> compute_summed_mass_matrix(Discretization1DPtr disc)
{
  int numEdges    = count_valid(disc->get_mesh()->get_edges());
  auto& quad      = disc->get_quadrature();
  auto& basisVals = disc->get_quad_basis_values();
  auto& dxdxi     = *(disc->getd_xdxi());

  std::vector<double> summedMassMatrix(disc->get_num_dofs(), 0), elContrib(basisVals.get_num_nodes());
  for (int i = 0; i < numEdges; ++i)
  {
    auto edge = disc->get_mesh()->get_edges()[i];
    for (int j = 0; j < basisVals.get_num_nodes(); ++j)
      elContrib[j] = 0;

    for (int k = 0; k < quad.get_num_points(); ++k)
    {
      // TODO: there is a more efficient way to do this: the basis functions sum to one everywhere, so
      //       the j1 loop is unnecessary
      for (int j1 = 0; j1 < basisVals.get_num_nodes(); ++j1)
        for (int j2 = 0; j2 < basisVals.get_num_nodes(); ++j2)
          elContrib[j2] +=
              basisVals.get_value(j1, k) * quad.get_weights()[k] * dxdxi(edge, 0, k) * basisVals.get_value(j2, k);
    }

    add_element_values(disc, edge, elContrib, summedMassMatrix);
  }

  return summedMassMatrix;
}

std::vector<double> apply_mass_matrix(Discretization1DPtr disc, const std::vector<double>& valsMesh)
{
  int numEdges    = count_valid(disc->get_mesh()->get_edges());
  auto& quad      = disc->get_quadrature();
  auto& basisVals = disc->get_quad_basis_values();
  auto& dxdxi     = *(disc->getd_xdxi());

  std::vector<double> result(disc->get_num_dofs(), 0), valsEl(basisVals.get_num_nodes()),
      resultEl(basisVals.get_num_nodes());
  for (int i = 0; i < numEdges; ++i)
  {
    auto edge = disc->get_mesh()->get_edges()[i];
    for (int j = 0; j < basisVals.get_num_nodes(); ++j)
      resultEl[j] = 0;

    get_element_values(disc, edge, valsMesh, valsEl);

    for (int k = 0; k < quad.get_num_points(); ++k)
      for (int j1 = 0; j1 < basisVals.get_num_nodes(); ++j1)
        for (int j2 = 0; j2 < basisVals.get_num_nodes(); ++j2)
          resultEl[j1] += basisVals.get_value(j1, k) * quad.get_weights()[k] * dxdxi(edge, 0, k) *
                          basisVals.get_value(j2, k) * valsEl[j2];

    add_element_values(disc, edge, valsEl, result);
  }

  return result;
}

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
