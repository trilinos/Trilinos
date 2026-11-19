#include <Akri_AuxMetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_SideAreaVector.hpp>
#include <Akri_SurfaceNodalNormal.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <algorithm>

namespace krino
{

FieldRef register_surface_nodal_normal_field(stk::mesh::MetaData & meta,
    const std::string & name,
    const stk::mesh::Selector & surfaceSelector)
{
  const FieldType & vecType =
      (meta.spatial_dimension() == 3) ? FieldType::VECTOR_3D : FieldType::VECTOR_2D;
  return AuxMetaData::get(meta).register_stk_field(
      name, vecType, stk::topology::NODE_RANK, 1, 1, surfaceSelector);
}

void update_surface_nodal_normal(const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const FieldRef normalField,
    const stk::mesh::Selector & sideSelector,
    const stk::mesh::Selector & elementSelector)
{
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  const unsigned dim = meta.spatial_dimension();

  stk::mesh::field_fill(0.0, normalField);
  auto normalData = normalField.data<double, stk::mesh::ReadWrite>();

  // Each (side, element) pair is added once, by the element's owner, so the parallel sum gives the
  // full sum without aura.
  for (auto * bucketPtr :
      mesh.get_buckets(stk::topology::ELEMENT_RANK, elementSelector & meta.locally_owned_part()))
  {
    for (auto elem : *bucketPtr)
    {
      const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
      STK_ThrowRequireMsg(
          elemNodes.size() == dim + 1, "Surface nodal normal expects linear tri or tet elements.");
      for (auto side :
          StkMeshEntities{mesh.begin(elem, meta.side_rank()), mesh.end(elem, meta.side_rank())})
      {
        if (!sideSelector(mesh.bucket(side))) continue;
        const StkMeshEntities sideNodes{mesh.begin_nodes(side), mesh.end_nodes(side)};
        const stk::math::Vector3d areaVec =
            compute_line_or_tri_side_area_vector_out_of_element(mesh, coordsField, side, elem);
        for (auto node : sideNodes)
        {
          STK_ThrowRequireMsg(has_field_data(normalField, node),
              "Surface nodal normal " << normalField.name()
                                      << " is not defined on the nodes of all the selected sides.");
          auto normalVals = normalData.entity_values(node);
          for (unsigned d = 0; d < dim; ++d)
            normalVals(stk::mesh::ComponentIdx(d)) += areaVec[d];
        }
      }
    }
  }

  // Shared nodes hold the full sums after parallel_sum; ghosted nodes get the owner's value from
  // communicate_field_data.
  stk::mesh::parallel_sum(mesh, {&normalField.field()});
  stk::mesh::communicate_field_data(mesh, {&normalField.field()});
}

} // namespace krino
