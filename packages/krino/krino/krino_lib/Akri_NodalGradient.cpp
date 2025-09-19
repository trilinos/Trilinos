/*
 * Akri_NodalGradient.cpp
 *
 *  Created on: Sep 11, 2025
 *      Author: drnoble
 */
#include <string>
#include <stk_mesh/base/MetaData.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_AuxMetaData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalGradient.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <Akri_SimplexGradient.hpp>

namespace krino {

std::string get_nodal_gradient_field_name(const std::string & scalarFieldName)
{
  std::string gradientFieldName = "GRAD_" + scalarFieldName;
  return gradientFieldName;
}

FieldRef get_nodal_gradient_for_scalar_field(const stk::mesh::MetaData & meta, const FieldRef scalarField)
{
  return meta.get_field(stk::topology::NODE_RANK, get_nodal_gradient_field_name(scalarField.name()));
}

FieldRef register_nodal_gradient_for_scalar_field(stk::mesh::MetaData & meta, const FieldRef scalarField)
{
  const FieldType & vecType = (meta.spatial_dimension() == 3) ? FieldType::VECTOR_3D : FieldType::VECTOR_2D;
  AuxMetaData & auxMeta = krino::AuxMetaData::get(meta);
  FieldRef gradField = auxMeta.register_stk_field(get_nodal_gradient_field_name(scalarField.name()), vecType, stk::topology::NODE_RANK, 2, 1, stk::mesh::selectField(scalarField));
  return gradField;
}

static stk::math::Vector3d calculate_tri_or_tet_volume_weighted_gradient(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef scalarField, const StkMeshEntities & elemNodes)
{
  STK_ThrowRequireMsg(elemNodes.size() == 3 || elemNodes.size() == 4, "Current expect either linear tri or tet");

  stk::math::Vector3d grad;
  double volOrArea;
  if (elemNodes.size() == 3)
  {
    const std::array<stk::math::Vector3d,3> nodeCoords = get_triangle_vector(mesh, coordsField, elemNodes, 2);
    const std::array<double,3> nodeScalar = get_triangle_scalar(mesh, scalarField, elemNodes);
    calculate_triangle2d_gradient_and_area(nodeCoords, nodeScalar, grad, volOrArea);
    return volOrArea*grad;
  }

  const std::array<stk::math::Vector3d,4> nodeCoords = get_tetrahedron_vector(mesh, coordsField, elemNodes);
  const std::array<double,4> nodeScalar = get_tetrahedron_scalar(mesh, scalarField, elemNodes);
  calculate_tetrahedron_gradient_and_volume(nodeCoords, nodeScalar, grad, volOrArea);
  return volOrArea*grad;
}

double compute_tri_or_tet_volume(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Entity elem)
{
  const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
  STK_ThrowAssert(elemNodes.size() == 4 || elemNodes.size() == 3);
  if (elemNodes.size() == 4)
    return compute_tet_volume(get_tetrahedron_vector(mesh, coordsField, elemNodes));
  return compute_tri_volume(get_triangle_vector(mesh, coordsField, elemNodes, 2));
}

static double calculate_nodal_volume(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const stk::mesh::Selector & elemSelector, const stk::mesh::Entity node)
{
  double nodalVol = 0.;
  for (auto elem : StkMeshEntities{mesh.begin_elements(node), mesh.end_elements(node)})
    if (is_entity_selected(mesh, elemSelector, elem))
      nodalVol += compute_tri_or_tet_volume(mesh, coordsField, elem);
  return nodalVol;
}

void update_nodal_gradient(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef scalarField)
{
  update_nodal_gradient(mesh, coordsField, scalarField, get_nodal_gradient_for_scalar_field(mesh.mesh_meta_data(), scalarField));
}

void update_nodal_gradient(const stk::mesh::BulkData & mesh, const FieldRef coordsField, const FieldRef scalarField, const FieldRef gradientField)
{
  STK_ThrowRequire((mesh.is_automatic_aura_on() || mesh.parallel_size() == 1)); // Could be easily generalized, but currently uses aura during averaging

  // Compute lumped nodal gradient
  // Could be made cheaper by storing volume to use during averaging instead of recomputing it and this would also eliminate need for aura

  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const AuxMetaData & auxMeta = krino::AuxMetaData::get(mesh.mesh_meta_data());

  stk::mesh::field_fill(0.0, gradientField);

  stk::mesh::Selector activeFieldSelector = stk::mesh::selectField(gradientField) & auxMeta.active_part();

  for ( auto * bucketPtr : mesh.get_buckets(stk::topology::ELEMENT_RANK, activeFieldSelector) )
  {
    for (auto elem : *bucketPtr)
    {
      const StkMeshEntities elemNodes{mesh.begin_nodes(elem), mesh.end_nodes(elem)};
      const stk::math::Vector3d contrib = calculate_tri_or_tet_volume_weighted_gradient(mesh, coordsField, scalarField, elemNodes);
      for (auto node : elemNodes)
      {
        double * gradData = field_data<double>(gradientField, node);
        for (unsigned i=0; i<dim; ++i)
          gradData[i] += contrib[i];
      }
    }
  }

  stk::mesh::Selector activeFieldNotGhostSelector = stk::mesh::selectField(gradientField) & auxMeta.active_not_ghost_selector();

  for ( auto * bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, activeFieldNotGhostSelector) )
  {
    const auto & bucket = *bucketPtr;
    double * gradData = field_data<double>(gradientField, bucket);

    for (size_t i = 0; i < bucket.size(); ++i)
    {
      const double nodalVolume = calculate_nodal_volume(mesh, coordsField, activeFieldSelector, bucket[i]); // uses aura
      for ( unsigned d = 0; d < dim; d++ )
      {
        gradData[dim*i+d] /= nodalVolume;
      }
    }
  }

  stk::mesh::communicate_field_data(mesh, {&gradientField.field()});
}

}



