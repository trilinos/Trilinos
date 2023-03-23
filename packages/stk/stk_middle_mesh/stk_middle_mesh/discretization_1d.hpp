#ifndef DISCRETIZATION_1D
#define DISCRETIZATION_1D

#include "basis_1d.hpp"
#include "field.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "quadrature.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

mesh::FieldShape get_field_shape(int degree);

struct MeshSpec1D
{
    int numel;
    double xmin;
    double xmax;
};

std::shared_ptr<mesh::Mesh> make1_d_mesh(MeshSpec1D meshspec);

class Discretization1D
{
  public:
    int get_num_dofs() const { return m_numDofs; }

    int get_basis_degree() const { return m_basis.get_degree(); }

    const mesh::FieldPtr<int> get_dof_nums() const { return m_dofnums; }

    const mesh::FieldPtr<double> get_coords() const { return m_coords; }

    const mesh::FieldPtr<double> getd_xdxi() const { return m_dxdxi; }

    const Quadrature& get_quadrature() const { return m_quad; }

    const BasisFunctionValues& get_quad_basis_values() const { return m_quadValues; }

    const BasisFunctions& get_basis() const { return m_basis; }

    mesh::FieldShape get_field_shape() const { return m_fieldshape; }

    std::shared_ptr<mesh::Mesh> get_mesh() const { return m_mesh; }

  private:
    Discretization1D(MeshSpec1D meshspec, int basisDegree, int quadDegree)
      : m_mesh(make1_d_mesh(meshspec))
      , m_fieldshape(disc1d::impl::get_field_shape(basisDegree))
      , m_dofnums(mesh::create_field<int>(m_mesh, m_fieldshape, 1, 0))
      , m_coords(mesh::create_field<double>(m_mesh, m_fieldshape, 1, 0))
      , m_basis(basisDegree)
      , m_quad(quadDegree)
      , m_quadValues(m_basis, m_quad.get_points())
    {
      m_dxdxi = mesh::create_field<double>(m_mesh, mesh::FieldShape(0, 1, 0), m_quad.get_num_points(), 0);
      set_dof_nums();
      compute_coordinates();
      compute_metrics();
    }

    void set_dof_nums();

    void compute_coordinates();

    void compute_metrics();

    friend std::shared_ptr<Discretization1D> make_disc1_d(MeshSpec1D meshspec, int basisDegree, int quadDegree);

    std::shared_ptr<mesh::Mesh> m_mesh;
    mesh::FieldShape m_fieldshape;
    mesh::FieldPtr<int> m_dofnums;
    mesh::FieldPtr<double> m_coords;
    mesh::FieldPtr<double> m_dxdxi; // dxidx at the quadrature points
    int m_numDofs;
    BasisFunctions m_basis;
    Quadrature m_quad;
    BasisFunctionValues m_quadValues;
};

using Discretization1DPtr = std::shared_ptr<Discretization1D>;

Discretization1DPtr make_disc1_d(MeshSpec1D meshspec, int basisDegree, int quadDegree);

std::shared_ptr<mesh::Mesh> make1_d_mesh(MeshSpec1D meshspec);

// T must be callable as T(x)
template <typename T>
std::vector<double> make_grid_function(Discretization1DPtr disc, T func)
{
  auto& coords = *(disc->get_coords());
  std::vector<double> vals(disc->get_num_dofs());
  auto& dofsRef = *(disc->get_dof_nums());

  for (auto vert : disc->get_mesh()->get_vertices())
    vals[dofsRef(vert, 0, 0)] = func(coords(vert, 0, 0));

  if (disc->get_field_shape().count[1] == 1)
    for (auto edge : disc->get_mesh()->get_edges())
      vals[dofsRef(edge, 0, 0)] = func(coords(edge, 0, 0));

  return vals;
}

void get_element_values(Discretization1DPtr disc, mesh::MeshEntityPtr el, const std::vector<double>& valsMesh,
                        std::vector<double>& elVals);

void add_element_values(Discretization1DPtr disc, mesh::MeshEntityPtr el, const std::vector<double>& elVals,
                        std::vector<double>& valsMesh);

double integrate_function(Discretization1DPtr disc, const std::vector<double>& uVals);

// TODO: untested
/*
void computeMassMatrix(Discretization1DPtr disc);
*/

// Compute 1^T M, store result in a vector
std::vector<double> compute_summed_mass_matrix(Discretization1DPtr disc);

std::vector<double> apply_mass_matrix(Discretization1DPtr disc, const std::vector<double>& valsMesh);

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
#endif
