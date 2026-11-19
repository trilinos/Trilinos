#include <Akri_InterfaceInverseGradientArea.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_SimplexGradient.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <set>
#include <string>
#include <vector>

namespace krino
{

static std::string interface_inverse_gradient_area_field_name(const std::string & lsFieldName)
{
  return "SURFINVGRAD_" + lsFieldName;
}

FieldRef
register_interface_inverse_gradient_area_field(stk::mesh::MetaData & meta, const LS_Field & lsField)
{
  return AuxMetaData::get(meta).register_stk_field(
      interface_inverse_gradient_area_field_name(lsField.isovar.name()),
      FieldType::REAL,
      stk::topology::NODE_RANK,
      1,
      1,
      Phase_Support::get(meta).get_negative_levelset_interface_selector(lsField.identifier));
}

FieldRef get_interface_inverse_gradient_area_field(
    const stk::mesh::MetaData & meta, const LS_Field & lsField)
{
  const FieldRef field = meta.get_field(
      stk::topology::NODE_RANK, interface_inverse_gradient_area_field_name(lsField.isovar.name()));
  STK_ThrowRequireMsg(field.valid(),
      "The interface inverse gradient area of "
          << lsField.isovar.name()
          << " was not registered; construct a ConformingDecompositionWithSensitivities before "
             "committing the mesh.");
  return field;
}

static stk::math::Vector3d get_unsnapped_coordinates(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef snapDisplacements,
    const stk::mesh::Entity node)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  stk::math::Vector3d coords = get_vector_field(mesh, coordsField, node, dim);
  if (snapDisplacements.valid()) coords -= get_vector_field(mesh, snapDisplacements, node, dim);
  return coords;
}

// The gradient of the level set over the background (parent) element of a decomposed element, from
// the unsnapped coordinates and level set of its parent nodes, where the level set is linear.
static stk::math::Vector3d
compute_parent_element_levelset_gradient(const stk::mesh::BulkData & mesh,
    const CDFEM_Support & cdfemSupport,
    const FieldRef levelSetField,
    const stk::mesh::Entity elem)
{
  std::set<stk::mesh::Entity> parentNodes;
  for (auto node : StkMeshEntities{mesh.begin_nodes(elem), mesh.end_nodes(elem)})
    recursively_fill_parent_nodes(
        mesh, node, cdfemSupport.get_parent_node_ids_field(), parentNodes);
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  STK_ThrowRequireMsg(parentNodes.size() == dim + 1, "Expected a linear simplex parent element.");

  const FieldRef snapDisplacements = cdfemSupport.get_cdfem_snap_displacements_field();
  auto levelSetData = levelSetField.data<double>();
  std::vector<stk::math::Vector3d> x;
  std::vector<double> ls;
  for (auto node : parentNodes)
  {
    x.push_back(
        get_unsnapped_coordinates(mesh, cdfemSupport.get_coords_field(), snapDisplacements, node));
    ls.push_back(levelSetData.entity_values(node)());
  }
  return calculate_simplex_gradient(x, ls);
}

FieldRef
update_interface_inverse_gradient_area(const stk::mesh::BulkData & mesh, const LS_Field & lsField)
{
  return update_interface_inverse_gradient_area(mesh,
      lsField,
      Phase_Support::get(mesh.mesh_meta_data())
          .get_negative_levelset_interface_selector(lsField.identifier));
}

FieldRef update_interface_inverse_gradient_area(const stk::mesh::BulkData & mesh,
    const LS_Field & lsField,
    const stk::mesh::Selector & interfaceSideSelector)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(meta);
  const Phase_Support & phaseSupport = Phase_Support::get(meta);
  const FieldRef invGradAreaField = get_interface_inverse_gradient_area_field(meta, lsField);
  const FieldRef coordsField = cdfemSupport.get_coords_field();
  const stk::mesh::Selector & sideSelector = interfaceSideSelector;
  const stk::mesh::Selector elementSelector =
      phaseSupport.get_negative_levelset_block_selector(lsField.identifier) &
      AuxMetaData::get(meta).active_part() & meta.locally_owned_part();
  const unsigned dim = meta.spatial_dimension();

  stk::mesh::field_fill(0.0, invGradAreaField);
  auto invGradAreaData = invGradAreaField.data<double, stk::mesh::ReadWrite>();
  for (auto * bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector))
  {
    for (auto elem : *bucketPtr)
    {
      bool hasInterfaceSide = false;
      for (auto side :
          StkMeshEntities{mesh.begin(elem, meta.side_rank()), mesh.end(elem, meta.side_rank())})
        if (sideSelector(mesh.bucket(side))) hasInterfaceSide = true;
      if (!hasInterfaceSide) continue;

      const double gradMag =
          compute_parent_element_levelset_gradient(mesh, cdfemSupport, lsField.isovar, elem)
              .length();
      if (gradMag == 0.) continue;

      for (auto side :
          StkMeshEntities{mesh.begin(elem, meta.side_rank()), mesh.end(elem, meta.side_rank())})
      {
        if (!sideSelector(mesh.bucket(side))) continue;
        const StkMeshEntities sideNodes{mesh.begin_nodes(side), mesh.end_nodes(side)};
        const stk::math::Vector3d x0 = get_vector_field(mesh, coordsField, sideNodes[0], dim);
        const stk::math::Vector3d x1 = get_vector_field(mesh, coordsField, sideNodes[1], dim);
        const double area = (sideNodes.size() == 2) ? (x1 - x0).length() :
                                                      0.5 *
                Cross(x1 - x0, get_vector_field(mesh, coordsField, sideNodes[2], dim) - x0)
                    .length();
        for (auto node : sideNodes)
          invGradAreaData.entity_values(node)() += area / gradMag;
      }
    }
  }

  stk::mesh::parallel_sum(mesh, {&invGradAreaField.field()});
  stk::mesh::communicate_field_data(mesh, {&invGradAreaField.field()});
  return invGradAreaField;
}

} // namespace krino
