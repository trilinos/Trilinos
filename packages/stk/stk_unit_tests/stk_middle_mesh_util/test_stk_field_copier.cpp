#include "gtest/gtest.h"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_middle_mesh_util/stk_field_copier.hpp"

using namespace stk::middle_mesh;
namespace {

double getCoordValue(const utils::Point& pt)
{
  return pt.x + 2*pt.y + 3*pt.z;
}

void set_field(mesh::FieldPtr<double> meshFieldPtr)
{
  auto mesh = meshFieldPtr->get_mesh();
  int numNodesPerEntity = meshFieldPtr->get_field_shape().get_num_nodes(0);
  int numCompPerNode = meshFieldPtr->get_num_comp();
  auto& field = *meshFieldPtr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      double valStart = getCoordValue(pt);

      for (int i=0; i < numNodesPerEntity; ++i)
        for (int j=0; j < numCompPerNode; ++j)
          field(vert, i, j) = valStart + i * numCompPerNode + j;
    }
}

void set_field(std::shared_ptr<stk::mesh::BulkData> bulkDataPtr, stk::mesh::Field<double>& field)
{
  const stk::mesh::FieldBase& coordField = *(bulkDataPtr->mesh_meta_data_ptr()->coordinate_field());

  stk::mesh::Selector selector(field);
  const stk::mesh::BucketVector& buckets = bulkDataPtr->get_buckets(stk::topology::NODE_RANK, selector);

  auto coordFieldData = coordField.data<double>();
  auto fieldData = field.data<stk::mesh::OverwriteAll>();
  for (stk::mesh::Bucket* bucket : buckets)
    for (stk::mesh::Entity node : *bucket)
    {
      auto coordData = coordFieldData.entity_values(node);
      auto nodeFieldData = fieldData.entity_values(node);
      double valStart = getCoordValue({coordData(0_comp), coordData(1_comp), coordData(2_comp)});

      for (stk::mesh::CopyIdx i : nodeFieldData.copies())
        for (stk::mesh::ComponentIdx j : nodeFieldData.components())
          nodeFieldData(i, j) = valStart + static_cast<int>(i * nodeFieldData.num_components()) + static_cast<int>(j);
    }
}

void check_field(mesh::FieldPtr<double> meshFieldPtr)
{
  auto mesh = meshFieldPtr->get_mesh();
  int numNodesPerEntity = meshFieldPtr->get_field_shape().get_num_nodes(0);
  int numCompPerNode = meshFieldPtr->get_num_comp();
  auto& field = *meshFieldPtr;
  for (auto& vert : mesh->get_vertices())
    if (vert)
    {
      utils::Point pt = vert->get_point_orig(0);
      double valStart = getCoordValue(pt);

      for (int i=0; i < numNodesPerEntity; ++i)
        for (int j=0; j < numCompPerNode; ++j)
          EXPECT_EQ(field(vert, i, j), valStart + i *numCompPerNode + j);
    }
}

void check_field(std::shared_ptr<stk::mesh::BulkData> bulkDataPtr, stk::mesh::Field<double>& field)
{
  const stk::mesh::FieldBase& coordField = *(bulkDataPtr->mesh_meta_data_ptr()->coordinate_field());

  auto meshMetaDataPtr = bulkDataPtr->mesh_meta_data_ptr();
  stk::mesh::Selector selector(field & (meshMetaDataPtr->locally_owned_part() | meshMetaDataPtr->globally_shared_part()));
  const stk::mesh::BucketVector& buckets = bulkDataPtr->get_buckets(stk::topology::NODE_RANK, selector);

  auto coordFieldData = coordField.data<double>();
  auto fieldData      = field.data();
  for (stk::mesh::Bucket* bucket : buckets)
    for (stk::mesh::Entity node : *bucket)
    {
      auto coordData = coordFieldData.entity_values(node);
      auto nodeFieldData = fieldData.entity_values(node);
      double valStart = getCoordValue({coordData(0_comp), coordData(1_comp), coordData(2_comp)});

      for (stk::mesh::CopyIdx i : nodeFieldData.copies())
        for (stk::mesh::ComponentIdx j : nodeFieldData.components())
          EXPECT_EQ(nodeFieldData(i, j), valStart + static_cast<int>(i * nodeFieldData.num_components()) + static_cast<int>(j));

    }
}

}

TEST(StkFieldCopier, MiddleMeshToStk)
{

  std::string meshFileName1 = "generated:3x3x2|sideset:Z|bbox:0,0,0,1,1,1";
  std::string partName1 = "surface_1";
  stk_interface::StkMeshCreator creator1(meshFileName1, "RCB", MPI_COMM_WORLD);
  stk_interface::MeshPart meshPart = creator1.create_mesh_from_part(partName1);

  mesh::FieldPtr<double> meshField = mesh::create_field<double>(meshPart.mesh, mesh::FieldShape(2, 0, 0), 3);
  set_field(meshField);

  stk_interface::StkFieldCopier copier(creator1.get_bulk_data_ptr(), meshPart.stkPart, meshPart.mesh, meshPart.stkVertField);
  stk::mesh::Field<double>* stkField = copier.create_stk_field(meshField, "stkfield");
  copier.copy(meshField, *stkField);
  check_field(creator1.get_bulk_data_ptr(), *stkField);

  meshField->set(0);
  copier.copy(*stkField, meshField);
  check_field(meshField);
}

TEST(StkFieldCopier, StkToMiddleMesh)
{

  std::string meshFileName1 = "generated:3x3x2|sideset:Z|bbox:0,0,0,1,1,1";
  std::string partName1 = "surface_1";
  stk_interface::StkMeshCreator creator1(meshFileName1, "NONE", MPI_COMM_WORLD);
  stk_interface::MeshPart meshPart = creator1.create_mesh_from_part(partName1);

  stk::mesh::Field<double>* stkField = &(creator1.get_meta_data_ptr()->declare_field<double>(stk::topology::NODE_RANK, "stkfield"));
  stk::mesh::put_field_on_mesh(*stkField, *(meshPart.stkPart), 3, 2, 0);
  //mesh::FieldPtr<double> meshField = mesh::create_field<double>(meshPart.mesh, mesh::FieldShape(2, 0, 0), 3);
  set_field(creator1.get_bulk_data_ptr(), *stkField);

  stk_interface::StkFieldCopier copier(creator1.get_bulk_data_ptr(), meshPart.stkPart, meshPart.mesh, meshPart.stkVertField);
  //stk::mesh::Field<double>* stkField = copier.create_stk_field(meshField, "stkfield");
  mesh::FieldPtr<double> meshField = copier.create_middle_mesh_field(*stkField);
  copier.copy(*stkField, meshField);
  check_field(meshField);

  copier.copy(meshField, *stkField);
  check_field(creator1.get_bulk_data_ptr(), *stkField);
}
