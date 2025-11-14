#include "stk_field_copier.hpp"

namespace stk {
namespace middle_mesh {
namespace stk_interface {



StkFieldCopier::StkFieldCopier(std::shared_ptr<stk::mesh::BulkData> bulkDataPtr, stk::mesh::Part* part,
                               std::shared_ptr<mesh::Mesh> middleMesh,
                               stk::mesh::Field<VertIdType>* stkVertField) :
  m_bulkDataPtr(bulkDataPtr),
  m_part(part),
  m_middleMesh(middleMesh),
  m_stkVertField(stkVertField)
{}

mesh::FieldPtr<double> StkFieldCopier::create_middle_mesh_field(const stk::mesh::Field<double>& stkField)
{
  std::pair<mesh::FieldShape, int> field_shape_and_comp = get_field_shape_and_num_components(stkField);
  auto meshField = mesh::create_field<double>(m_middleMesh, field_shape_and_comp.first, field_shape_and_comp.second);
  check_field_shapes(stkField, meshField);

  return meshField;
}

stk::mesh::Field<double>* StkFieldCopier::create_stk_field(mesh::FieldPtr<double> middleMeshField, const std::string& name)
{
  auto metaData = m_bulkDataPtr->mesh_meta_data_ptr();
  stk::mesh::Field<double>* stkField = &(metaData->declare_field<double>(stk::topology::NODE_RANK, name));
  stk::mesh::put_field_on_mesh(*stkField, *m_part, middleMeshField->get_num_comp(),
                                middleMeshField->get_field_shape().get_num_nodes(0), 0);

  return stkField;
}

void StkFieldCopier::copy(const stk::mesh::Field<double>& stkField, mesh::FieldPtr<double>& middleMeshFieldPtr)
{
  check_field_shapes(stkField, middleMeshFieldPtr);

  auto meshMetaDataPtr = m_bulkDataPtr->mesh_meta_data_ptr();
  stk::mesh::Selector selector(stkField & (meshMetaDataPtr->locally_owned_part() | meshMetaDataPtr->globally_shared_part()));
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::NODE_RANK, selector);

  int numNodesPerEntity = middleMeshFieldPtr->get_field_shape().get_num_nodes(0);
  int numCompPerNode = middleMeshFieldPtr->get_num_comp();
  auto& middleMeshField = *middleMeshFieldPtr;

  auto stkVertFieldData = m_stkVertField->data();
  auto stkFieldData     = stkField.data();
  for (const stk::mesh::Bucket* bucket : buckets)
    for (stk::mesh::Entity node : *bucket)
    {
      VertIdType meshVertId = stkVertFieldData.entity_values(node)(0_comp);
      mesh::MeshEntityPtr vert = m_middleMesh->get_vertices()[meshVertId];

      auto stkFieldDataForNode = stkFieldData.entity_values(node);
      for (stk::mesh::CopyIdx i=0_copy; i < numNodesPerEntity; ++i)
        for (stk::mesh::ComponentIdx j=0_comp; j < numCompPerNode; ++j)
          middleMeshField(vert, i, j) = stkFieldDataForNode(i, j);
    }
}

void StkFieldCopier::copy(const mesh::FieldPtr<double> middleMeshFieldPtr, stk::mesh::Field<double>& stkField)
{
  check_field_shapes(stkField, middleMeshFieldPtr);

  auto meshMetaDataPtr = m_bulkDataPtr->mesh_meta_data_ptr();
  stk::mesh::Selector selector(stkField & (meshMetaDataPtr->locally_owned_part() | meshMetaDataPtr->globally_shared_part()));
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::NODE_RANK, selector);

  int numNodesPerEntity = middleMeshFieldPtr->get_field_shape().get_num_nodes(0);
  int numCompPerNode = middleMeshFieldPtr->get_num_comp();
  auto& middleMeshField = *middleMeshFieldPtr;

  auto stkVertFieldData = m_stkVertField->data();
  auto stkFieldData = stkField.data<stk::mesh::ReadWrite>();
  for (const stk::mesh::Bucket* bucket : buckets)
    for (stk::mesh::Entity node : *bucket)
    {
      VertIdType meshVertId = stkVertFieldData.entity_values(node)(0_comp);
      mesh::MeshEntityPtr vert = m_middleMesh->get_vertices()[meshVertId];

      auto stkFieldDataForNode = stkFieldData.entity_values(node);
      for (stk::mesh::CopyIdx i=0_copy; i < numNodesPerEntity; ++i)
        for (stk::mesh::ComponentIdx j=0_comp; j < numCompPerNode; ++j)
          stkFieldDataForNode(i, j) = middleMeshField(vert, i, j);
    }
}


void StkFieldCopier::check_field_shapes(const stk::mesh::Field<double>& stkField, const mesh::FieldPtr<double> meshField)
{
  std::pair<int, int> stk_field_dims = get_field_shape(stkField);
  mesh::FieldShape fshape = meshField->get_field_shape();

  if (stkField.entity_rank() != stk::topology::NODE_RANK)
  {
    std::stringstream ss;
    ss << "Field is defined on " << stkField.entity_rank() << ", but only NODE_RANK fields are supported";
    throw std::runtime_error(ss.str());
  }

  if (fshape.get_num_nodes(0) != stk_field_dims.second)
  {
    throw std::runtime_error(
      std::string("Field shapes not compatible: stk field has ") + std::to_string(stk_field_dims.second) +
      " nodes per entity, while the middle mesh field has " + std::to_string(fshape.get_num_nodes(0))
    );
  }

  if (meshField->get_num_comp() != stk_field_dims.first)
  {
    throw std::runtime_error(
      std::string("Field shapes not compatible: stk field has ") + std::to_string(stk_field_dims.second) +
      " components per node, while the middle mesh field has " + std::to_string(meshField->get_num_comp())
    );
  }
}

std::pair<mesh::FieldShape, int> StkFieldCopier::get_field_shape_and_num_components(const stk::mesh::Field<double>& stkField)
{
  std::pair<int, int> dims = get_field_shape(stkField);

  // stk fields are column major for some strange reason (ie. the first index is contiguous in memory).
  // Make that data be contiguous in memory in the mesh field as well.
  // This requires the first dimension of the stk field to be the second dimension of the mesh field
  return std::make_pair(mesh::FieldShape(dims.second, 0, 0), dims.first);
}

std::pair<int, int> StkFieldCopier::get_field_shape(const stk::mesh::Field<double>& stkField)
{
  stk::mesh::Selector selector(stkField);
  const stk::mesh::BucketVector& buckets = m_bulkDataPtr->get_buckets(stk::topology::NODE_RANK, selector);
  if (buckets.size() == 0)
    throw std::runtime_error("Field has zero node-rank buckets.  Only node fields are supported");


  int dim0 = stk::mesh::field_extent0_per_entity(stkField, *(buckets[0]));
  int dim1 = stk::mesh::field_extent1_per_entity(stkField, *(buckets[0]));

  return std::make_pair(dim0, dim1);
}


}
}
}
