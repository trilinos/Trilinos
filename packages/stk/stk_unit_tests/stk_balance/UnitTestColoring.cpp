#include <gtest/gtest.h>

#include "stk_balance/balanceUtils.hpp"
#include "stk_io/IossBridge.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_io/WriteMesh.hpp"
#include "stk_mesh/base/FindRestriction.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_unit_test_utils/BuildMesh.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "stk_unit_test_utils/stk_mesh_fixtures/heterogeneous_mesh.hpp"
#include "stk_util/diag/StringUtil.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include <stk_balance/balance.hpp>

namespace {
using stk::unit_test_util::build_mesh;

class BasicColoring : public stk::unit_test_util::MeshFixture {};

void test_adjacent_elements_have_different_coloring(const stk::mesh::BulkData& bulk)
{
  stk::mesh::EntityVector nodes;
  const stk::mesh::MetaData& metaData = bulk.mesh_meta_data();
  stk::mesh::get_selected_entities(metaData.locally_owned_part(), bulk.buckets(stk::topology::NODE_RANK), nodes);

  for (stk::mesh::Entity& node : nodes)
  {
    std::set<stk::mesh::Part*> coloringParts;

    unsigned numElems = bulk.num_elements(node);
    unsigned numLocallyOwnedElems = 0;
    const stk::mesh::Entity* elems = bulk.begin_elements(node);
    for (unsigned i = 0; i < numElems; ++i)
    {
      stk::mesh::Entity elem = elems[i];
      if (bulk.bucket(elem).owned()) {
        ++numLocallyOwnedElems;
        stk::mesh::Part* colorPart = stk::balance::get_coloring_part(bulk, elem);
        if (nullptr != colorPart) {
          coloringParts.insert(colorPart);
        }
      }
    }
    EXPECT_EQ(numLocallyOwnedElems, coloringParts.size());
  }
}

TEST_F(BasicColoring, checkAdjacentElementsHaveDifferentColors)
{
  if (stk::parallel_machine_size(get_comm()) > 3) return;
  setup_mesh("generated:3x3x3", stk::mesh::BulkData::AUTO_AURA);

  stk::balance::BasicColoringSettings coloringSettings;
  bool meshIsColored = stk::balance::colorStkMesh(coloringSettings, *bulkData);
  EXPECT_TRUE(meshIsColored);
  test_adjacent_elements_have_different_coloring(*bulkData);
}

TEST_F(BasicColoring, checkForCorrectColors)
{
  if (stk::parallel_machine_size(get_comm()) > 1) return;
  setup_mesh("generated:3x3x3", stk::mesh::BulkData::AUTO_AURA);

  stk::balance::BasicColoringSettings coloringSettings;
  bool meshIsColored = stk::balance::colorStkMesh(coloringSettings, *bulkData);
  EXPECT_TRUE(meshIsColored);

  stk::mesh::PartVector coloringParts;
  stk::balance::fill_coloring_parts(*metaData, coloringParts);

  size_t goldNumberOfColors = 9;
  EXPECT_EQ(goldNumberOfColors, coloringParts.size());

  for (size_t i = 1; i <= goldNumberOfColors; ++i)
  {
    std::string goldPartName = stk::balance::construct_coloring_part_name(i, get_meta().locally_owned_part());
    bool colorFound = false;
    for (stk::mesh::Part* part : coloringParts)
    {
      if (part->name() == goldPartName)
      {
        colorFound = true;
        break;
      }
    }
    EXPECT_TRUE(colorFound);
  }
}

class ColorMeshWithColoringFieldsSettings : public stk::balance::BalanceSettings
{
public:
  ColorMeshWithColoringFieldsSettings() {}

  virtual GraphOption getGraphOption() const
  {
    return BalanceSettings::COLOR_MESH_AND_OUTPUT_COLOR_FIELDS;
  }
};

void declare_color_fields(stk::mesh::MetaData& meta)
{
  const stk::mesh::PartVector& parts = meta.get_parts();
  for(stk::mesh::Part* part : parts)
  {
    if(part->primary_entity_rank() == stk::topology::ELEM_RANK && stk::io::is_part_io_part(*part))
    {
      stk::topology topo = part->topology();
      stk::mesh::Part& rootTopologyPart = meta.get_topology_root_part(topo);
      stk::mesh::Field<int>& colorField = meta.declare_field<int>(stk::topology::ELEMENT_RANK,
                                                                  rootTopologyPart.topology().name() + "coloring");
      stk::mesh::put_field_on_mesh(colorField, rootTopologyPart, nullptr);
    }
  }
}

void check_coloring_by_topology(const stk::mesh::MetaData& meta, std::map<stk::topology, unsigned>& goldTopologyCounts)
{
  unsigned goldNumberOfColors = 0;
  for (auto& topoCount : goldTopologyCounts)
  {
    goldNumberOfColors += topoCount.second;

    stk::mesh::PartVector topoColoringParts;
    stk::topology topo = topoCount.first;
    stk::balance::fill_coloring_parts_with_topology(meta, topo, topoColoringParts);
    EXPECT_EQ(goldTopologyCounts[topo], topoColoringParts.size()) << " incorrect size for topology " << topo;
  }

  stk::mesh::PartVector coloringParts;
  stk::balance::fill_coloring_parts(meta, coloringParts);
  EXPECT_EQ(goldNumberOfColors, coloringParts.size());

}

TEST(ColorByTopology, colorHeterogeneousMesh)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  stk::mesh::fixtures::VectorFieldType & node_coord = meta.declare_field<double>(stk::topology::NODE_RANK,
                                                                                                "coordinates");
  stk::mesh::put_field_on_mesh(node_coord, meta.universal_part(), 3, nullptr);

  stk::mesh::fixtures::heterogeneous_mesh_meta_data( meta , node_coord );
  declare_color_fields(meta);
  meta.commit();

  stk::mesh::fixtures::heterogeneous_mesh_bulk_data( *bulk , node_coord );

  stk::balance::BasicColoringByTopologySettings coloringSettings;
  //    ColorMeshWithColoringFieldsSettings coloringSettings;
  bool meshIsColored = stk::balance::colorStkMesh(coloringSettings, *bulk);
  EXPECT_TRUE(meshIsColored);

  std::map<stk::topology, unsigned> goldTopologyCounts;
  goldTopologyCounts[stk::topology::SHELL_TRIANGLE_3] = 3;
  goldTopologyCounts[stk::topology::SHELL_QUADRILATERAL_4] = 2;
  goldTopologyCounts[stk::topology::PYRAMID_5] = 2;
  goldTopologyCounts[stk::topology::TETRAHEDRON_4] = 3;
  goldTopologyCounts[stk::topology::HEXAHEDRON_8] = 2;
  goldTopologyCounts[stk::topology::WEDGE_6] = 3;

  check_coloring_by_topology(meta, goldTopologyCounts);
}

typedef stk::mesh::Field<double> Vector2dFieldType;

void quad_tri_mesh_meta_data(stk::mesh::MetaData & meta_data,
                             const Vector2dFieldType & node_coord)
{
  stk::mesh::Part & universal        = meta_data.universal_part();
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("quad4", stk::topology::QUADRILATERAL_4_2D));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("tri3", stk::topology::TRIANGLE_3_2D));

  const stk::mesh::FieldBase::Restriction & res =
      stk::mesh::find_restriction(node_coord, stk::topology::NODE_RANK , universal );

  if ( res.num_scalars_per_entity() != 2 ) {
    std::ostringstream msg ;
    msg << "stk_mesh/unit_tests/quad_tri_mesh_meta_data FAILED, coordinate dimension must be 3 != " << res.num_scalars_per_entity() ;
    throw std::runtime_error( msg.str() );
  }
}

class Color2DMesh : public stk::unit_test_util::MeshFixture2D {};

TEST_F(Color2DMesh, colorHeterogeneousMeshWithQuadsSurroundingTriangles)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) return;

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  Vector2dFieldType & node_coord = get_meta().declare_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field_on_mesh(node_coord, get_meta().universal_part(), 2, nullptr);

  quad_tri_mesh_meta_data(get_meta() , node_coord);

  declare_color_fields(get_meta());

  std::string meshDesc = "0,1,QUAD_4_2D,  1,  5,  6,  2,block_1\n"
                         "0,2,QUAD_4_2D,  2,  6,  7,  3,block_1\n"
                         "0,3,QUAD_4_2D,  3,  7,  8,  4,block_1\n"
                         "0,4,QUAD_4_2D,  5,  9, 10,  6,block_1\n"
                         "0,5,QUAD_4_2D,  7, 11, 12,  8,block_1\n"
                         "0,6,QUAD_4_2D,  9, 13, 14, 10,block_1\n"
                         "0,7,QUAD_4_2D, 10, 14, 15, 11,block_1\n"
                         "0,8,QUAD_4_2D, 11, 15, 16, 12,block_1\n"
                         "0,9, TRI_3_2D,  6, 10,  7    ,block_2\n"
                         "0,10,TRI_3_2D, 10, 11,  7    ,block_2\n";

  std::vector<double> coordinates = { 0,3, 1,3, 2,3, 3,3, 0,2, 1,2, 2,2, 3,2,
                                      0,1, 1,1, 2,1, 3,1, 0,0, 1,0, 2,0, 3,0 };

  stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

  ColorMeshWithColoringFieldsSettings coloringSettings;
  bool meshIsColored = stk::balance::colorStkMesh(coloringSettings, get_bulk());
  EXPECT_TRUE(meshIsColored);

  std::map<stk::topology, unsigned> goldTopologyCounts;
  goldTopologyCounts[stk::topology::QUADRILATERAL_4_2D] = 3;
  goldTopologyCounts[stk::topology::TRIANGLE_3_2D] = 2;

  check_coloring_by_topology(get_meta(), goldTopologyCounts);
}

}
