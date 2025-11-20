#include "stk_middle_mesh_util/field_output_adaptor.hpp"
#include "gtest/gtest.h"
#include "stk_middle_mesh_util/exodus_writer.hpp"
#include "stk_middle_mesh_util/stk_interface.hpp"
#include "stk_middle_mesh/create_mesh.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_topology/topology.hpp"


namespace {

using namespace stk::middle_mesh;

std::shared_ptr<mesh::Mesh> create_simple_mesh(bool makeTriangles)
{
  mesh::impl::MeshSpec spec;
  spec.xmin = 0;   spec.ymin = 0;
  spec.xmax = 1;   spec.ymax = 1;
  spec.numelX = 3; spec.numelY = 3;
  auto f = [](const utils::Point& pt) { return pt; };
  return mesh::impl::create_mesh(spec, f, MPI_COMM_WORLD, makeTriangles);
}

void read_stk_mesh(const std::string& fname, stk::mesh::BulkData& bulkData)
{
  stk::io::StkMeshIoBroker reader(bulkData.parallel());
  reader.set_bulk_data(bulkData);
  reader.add_mesh_database(fname, stk::io::READ_MESH);
  reader.create_input_mesh();
  reader.add_all_mesh_fields_as_input_fields();
  // reader.populate_bulk_data();
  reader.populate_mesh();
  reader.populate_field_data();
  reader.read_defined_input_fields(1);
}

int count_entities(stk::mesh::BulkData& bulkData, stk::mesh::EntityRank rank)
{
  stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  return stk::mesh::count_entities(bulkData, rank, metaData.universal_part());
  /*
  int nentities = 0;
  for (stk::mesh::Bucket* bucket : bulkData.get_buckets(rank, metaData.universal_part()))
    nentities += bucket->size();

  return nentities;
  */
}

void check_vertex_coordinates(std::unique_ptr<stk::mesh::BulkData>& bulkDataPtr,
                              std::shared_ptr<mesh::Mesh> mesh)
{
  stk::mesh::MetaData& metaData = bulkDataPtr->mesh_meta_data();
  int meshElIdx = 0;
  for (stk::mesh::Bucket* bucket : bulkDataPtr->get_buckets(stk::topology::ELEM_RANK, metaData.universal_part()))
    for (size_t i=0; i < bucket->size(); ++i)
    {
      stk::mesh::Entity stkEl = (*bucket)[i];
      mesh::MeshEntityPtr meshEl = mesh->get_elements()[meshElIdx];

      std::vector<stk::mesh::Entity> stkVerts(bulkDataPtr->begin_nodes(stkEl), bulkDataPtr->end_nodes(stkEl));

      std::vector<mesh::MeshEntityPtr> meshVerts(4);
      int nverts = mesh::get_downward(meshEl, 0, meshVerts.data());
      meshVerts.resize(nverts);

      EXPECT_EQ(meshVerts.size(), stkVerts.size());
      const stk::mesh::FieldBase& stkCoordField = *(metaData.coordinate_field());
      auto stkCoordFieldData = stkCoordField.data<double>();
      for (int j=0; j < nverts; ++j)
      {
        auto stkVertCoords = stkCoordFieldData.entity_values(stkVerts[j]);
        utils::Point meshVertCoords = meshVerts[j]->get_point_orig(0);

        for (stk::mesh::ComponentIdx d=0_comp; d < 3_comp; ++d)
        {
          EXPECT_NEAR(stkVertCoords(d), meshVertCoords[d], 1e-13);
        }
      }

      meshElIdx++;
    }
}

mesh::FieldPtr<double> make_simple_field(std::shared_ptr<mesh::Mesh> mesh, int dim)
{
  mesh::FieldShape fshape;
  fshape.count[dim] = 2;
  auto field = mesh::create_field<double>(mesh, fshape, 4);
  int idx = 0;
  for (auto& entity : mesh->get_mesh_entities(dim))
    if (entity)
    {
      double valBase = idx;

      for (int i=0; i < 2; ++i)
        for (int j=0; j < 4; ++j)
        {
          int valOffset = i * 4 + j;
          (*field)(entity, i, j) = valBase + valOffset;
        }

      idx++;
    }

  return field;
}

void check_simple_field(std::unique_ptr<stk::mesh::BulkData>& bulkDataPtr, int dim, stk::mesh::FieldBase& field)
{
  int idx = 0;
  stk::mesh::EntityRank rank = dim == 0 ? stk::topology::NODE_RANK : stk::topology::ELEM_RANK;
  auto fieldData = field.data<double>();
  for (stk::mesh::Bucket* bucket : bulkDataPtr->buckets(rank))
    for (const stk::mesh::Entity& entity : *bucket)
    {
      auto entityData = fieldData.entity_values(entity);

      double valBase = idx;
      for (stk::mesh::CopyIdx i=0_copy; i < 2; ++i)
        for (stk::mesh::ComponentIdx j=0_comp; j < 4; ++j)
        {
          int valOffset = static_cast<int>(i*4) + j;
          EXPECT_NEAR(entityData(i, j), valBase + valOffset, 1e-13);
        }

      idx++;
    }
}

} // namespace

TEST(ExodusWriter, EntityCounts)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  for (bool makeTriangles : {true, false})
  {
    auto mesh = create_simple_mesh(makeTriangles);

    std::string fname = "tmp.exo";
    stk_interface::impl::ExodusWriter writer(mesh);
    writer.write(fname);

    auto bulkDataPtr  = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
    read_stk_mesh(fname, *bulkDataPtr);

    EXPECT_EQ(count_entities(*bulkDataPtr, stk::topology::NODE_RANK), 16);
    EXPECT_EQ(count_entities(*bulkDataPtr, stk::topology::ELEM_RANK), makeTriangles ? 18 : 9);
    check_vertex_coordinates(bulkDataPtr, mesh);
  }
}

TEST(ExodusWriter, FieldValuesVert)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  for (bool makeTriangles : {true})
  {
    auto mesh = create_simple_mesh(makeTriangles);
    auto field = make_simple_field(mesh, 0);

    std::string fieldName = "test_field";
    auto fieldOutputAdaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field, fieldName);

    std::string fname = "tmp.exo";
    stk_interface::impl::ExodusWriter writer(mesh, {fieldOutputAdaptor});
    writer.write(fname);

    auto bulkDataPtr  = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
    read_stk_mesh(fname, *bulkDataPtr);

    stk::mesh::MetaData& metaData = bulkDataPtr->mesh_meta_data();
    stk::mesh::FieldBase* stkField = metaData.get_field(stk::topology::NODE_RANK, fieldName + "_verts");

    EXPECT_NE(stkField, nullptr);
    check_simple_field(bulkDataPtr, 0, *stkField);
  }
}


TEST(ExodusWriter, FieldValuesElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  for (bool makeTriangles : {true, false})
  {
    auto mesh = create_simple_mesh(makeTriangles);
    auto field = make_simple_field(mesh, 2);

    std::string fieldName = "test_field";
    auto fieldOutputAdaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field, fieldName);

    std::string fname = "tmp.exo";
    stk_interface::impl::ExodusWriter writer(mesh, {fieldOutputAdaptor});
    writer.write(fname);

    auto bulkDataPtr  = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
    read_stk_mesh(fname, *bulkDataPtr);

    stk::mesh::MetaData& metaData = bulkDataPtr->mesh_meta_data();
    stk::mesh::FieldBase* stkField = metaData.get_field(stk::topology::ELEM_RANK, fieldName + "_elems");

    EXPECT_NE(stkField, nullptr);
    check_simple_field(bulkDataPtr, 2, *stkField);
  }
}


TEST(ExodusWriter, FieldValuesDeletedElement)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 1)
    GTEST_SKIP();

  for (bool makeTriangles : {true, false})
  {
    auto mesh = create_simple_mesh(makeTriangles);
    mesh->delete_face(mesh->get_elements()[5]);
    auto field = make_simple_field(mesh, 2);

    std::string fieldName = "test_field";
    auto fieldOutputAdaptor = std::make_shared<stk_interface::impl::FieldOutputAdaptorDouble>(field, fieldName);

    std::string fname = "tmp.exo";
    stk_interface::impl::ExodusWriter writer(mesh, {fieldOutputAdaptor});
    writer.write(fname);

    auto bulkDataPtr  = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
    read_stk_mesh(fname, *bulkDataPtr);

    stk::mesh::MetaData& metaData = bulkDataPtr->mesh_meta_data();
    stk::mesh::FieldBase* stkField = metaData.get_field(stk::topology::ELEM_RANK, fieldName + "_elems");

    EXPECT_NE(stkField, nullptr);
    check_simple_field(bulkDataPtr, 2, *stkField);
  }
}
