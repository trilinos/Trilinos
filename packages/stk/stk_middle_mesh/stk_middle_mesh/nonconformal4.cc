#include "nonconformal4.h"
#include "mesh_projection_calculator.h"
#include "middle_grid_constraint_generator.h"
#include "middle_grid_triangulator.h"

#include "mesh_io.h" //TODO: DEBUGGING

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

// inverts a the field field_ptr, defined on mesh_in, creating a new field on mesh_other.
// field_ptr must give the elements on mesh_other
mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr>
invert_classification_field(std::shared_ptr<mesh::Mesh> meshIn, std::shared_ptr<mesh::Mesh> meshOther,
                            mesh::FieldPtr<mesh::MeshEntityPtr> fieldPtr)
{
  mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> fieldInvPtr =
      mesh::create_variable_size_field<mesh::MeshEntityPtr>(meshOther, mesh::impl::FieldShape(0, 0, 1));
  auto& field    = *fieldPtr;
  auto& fieldInv = *fieldInvPtr;
  for (auto& el : meshIn->get_elements())
    if (el)
    {
      mesh::MeshEntityPtr elOther = field(el, 0, 0);
      fieldInv.insert(elOther, 0, el);
    }

  return fieldInvPtr;
}

std::shared_ptr<mesh::Mesh> Nonconformal4::create()
{
  {
    MeshProjectionCalculator meshProjection(m_mesh1, m_mesh2, m_meshIn, m_relationalData, m_classifier,
                                            m_edgeTracerTolerances);
    meshProjection.project();
  }

  {
    MiddleGridConstraintGenerator generator(m_mesh1, m_mesh2, m_meshIn, m_relationalData);
    generator.generate();
  }

  {
    MiddleGridTriangulator triangulator(m_mesh1, m_mesh2, m_meshIn, m_relationalData, m_classifier);
    triangulator.triangulate();
  }

  apply_geometry_improvers();

  return m_meshIn;
}

void Nonconformal4::apply_geometry_improvers()
{
  for (auto enumVal : m_geometryImprovers)
  {
    std::shared_ptr<GeometryImprover> improver = geometry_improver_factory(
        enumVal, m_mesh1, m_mesh2, m_meshIn, m_relationalData, m_classifier->get_normal_field());
    improver->run();
  }
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
