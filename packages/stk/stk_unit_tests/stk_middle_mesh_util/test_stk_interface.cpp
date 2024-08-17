#ifdef STK_BUILT_FOR_SIERRA

#include "stk_middle_mesh_util/constants.hpp"
#include "stk_middle_mesh_util/stk_interface.hpp"
#include "stk_middle_mesh/utils.hpp"
#include "gtest/gtest.h"

#include <ostream>

namespace stk {
namespace middle_mesh {
namespace impl {

namespace mesh {
std::ostream& printit(std::ostream& os, const stk::mesh::SideSetEntry& e, stk::mesh::BulkData& bulkData)
{
  os << bulkData.identifier(e.element) << ", " << e.side;
  return os;
}
} // namespace mesh

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
}

TEST(StkInterface, twoToThree)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) > 1)
    GTEST_SKIP();

  std::vector<middle_mesh::mesh::impl::MeshInput::NamePair> interfaces{
      {"left_block_interface", "right_block_interface"}};

  std::string fnameOut  = "nonconformal_cubes_out.exo";
  std::string fnameOut2 = "nonconformal_cubes_surface_out.exo";
  middle_mesh::mesh::impl::MeshInput input("nonconformal_cubes_ejoined.exo", fnameOut, fnameOut2, interfaces);

  stk_interface::impl::StkInterface stkInterface(input);

  NonconformalOpts opts;
  opts.enableSnapAndQuality = false;
  stkInterface.compute_interface(opts);
  stkInterface.write_output();

  // read output mesh and test it
  auto bulkDataPtr  = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();
  auto bulkData2Ptr = stk::mesh::MeshBuilder(MPI_COMM_WORLD).set_spatial_dimension(3).create();

  auto& bulkData  = *bulkDataPtr;
  auto& bulkData2 = *bulkData2Ptr;

  auto& metaData  = bulkDataPtr->mesh_meta_data();
  auto& metaData2 = bulkData2Ptr->mesh_meta_data();

  read_stk_mesh(fnameOut, bulkData);
  read_stk_mesh(fnameOut2, bulkData2);

  std::string partName = stk_interface::impl::NAME_PREFIX + "nonconformal_interface_" + interfaces[0].first + "_and_" +
                         interfaces[0].second;
  stk::mesh::Part* ncPart = metaData2.get_part(partName);

  // check number of entities
  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(*ncPart, bulkData2, entityCounts);
  EXPECT_EQ(entityCounts[stk::topology::NODE_RANK], 25u);
  EXPECT_EQ(entityCounts[stk::topology::ELEM_RANK], 32u);

  // check GID field
  std::map<stk::mesh::SideSetEntry, int> gidCountL;
  std::map<stk::mesh::SideSetEntry, int> gidCountR;

  using FType           = stk::mesh::Field<double>;
  std::string fieldName = stk_interface::impl::NAME_PREFIX + "nonconformal_interface_gid_field";

  FType* fPtr = metaData2.get_field<double>(stk::topology::ELEM_RANK, fieldName);

  ASSERT_NE(fPtr, nullptr);

  stk::mesh::Selector selPart            = *ncPart;
  stk::mesh::BucketVector const& buckets = bulkData2.get_buckets(stk::topology::ELEM_RANK, selPart);
  std::vector<stk::mesh::Entity> entities;
  stk::mesh::get_selected_entities(selPart, buckets, entities);

  for (auto& e : entities)
  {
    double* dataE           = stk::mesh::field_data(*fPtr, e);
    stk::mesh::Entity elemL = bulkData.get_entity(stk::topology::ELEM_RANK, dataE[0]);
    stk::mesh::Entity elemR = bulkData.get_entity(stk::topology::ELEM_RANK, dataE[2]);
    stk::mesh::SideSetEntry entryL(elemL, static_cast<int>(dataE[1] - 1));
    stk::mesh::SideSetEntry entryR(elemR, static_cast<int>(dataE[3] - 1));
    gidCountL[entryL] += 1;
    gidCountR[entryR] += 1;
  }

  for (auto& p : gidCountL)
    EXPECT_EQ(p.second, 8);

  for (auto& p : gidCountR)
    EXPECT_TRUE(p.second == 2 || p.second == 4 || p.second == 8);

  // check GID field contains only entities in the corresponding sidesets
  // and vice-versa
  stk::mesh::Part* partL      = metaData.get_part(interfaces[0].first);
  stk::mesh::SideSet sidesetL = bulkData.get_sideset(*partL);
  stk::mesh::Part* partR      = metaData.get_part(interfaces[0].second);
  stk::mesh::SideSet sidesetR = bulkData.get_sideset(*partR);

  for (auto& p : gidCountL)
    EXPECT_TRUE(sidesetL.contains(p.first));

  for (auto& p : gidCountR)
    EXPECT_TRUE(sidesetR.contains(p.first));

  for (auto& entry : sidesetL)
    EXPECT_EQ(gidCountL.count(entry), 1u);

  for (auto& entry : sidesetR)
    EXPECT_EQ(gidCountR.count(entry), 1u);
}

#endif
}
}
}
