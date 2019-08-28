#include <gtest/gtest.h>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_io/StkIoUtils.hpp>
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include "stk_unit_test_utils/FaceTestingUtils.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"
#include "stk_unit_test_utils/GenerateALefRAMesh.hpp"

namespace {

stk::mesh::Entity verify_and_get_face(stk::mesh::BulkData &bulk, stk::mesh::Part* ssPart)
{
    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(*ssPart, bulk.buckets(stk::topology::FACE_RANK), faces);
    EXPECT_EQ(1u, faces.size());
    return faces[0];
}

void check_polarity_of_modified_AA(stk::ParallelMachine pm,
                                   stk::unit_test_util::SidesetDirection direction,
                                   const std::string& name,
                                   bool isPositivePolarity)
{
    size_t spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData bulk(meta, pm);

    stk::mesh::Part* ssPart = stk::unit_test_util::create_AA_mesh_with_sideset(bulk, direction);
    stk::mesh::Entity  face = verify_and_get_face(bulk, ssPart);

    std::pair<bool,bool> sidesetPolarity = stk::io::is_positive_sideset_polarity(bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first);
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second);
}

void check_polarity_of_modified_AB(stk::ParallelMachine pm,
                                   stk::unit_test_util::SidesetDirection direction,
                                   const std::string& name,
                                   bool isPositivePolarity)
{
    size_t spatialDim = 3;
    stk::mesh::MetaData meta(spatialDim);
    stk::mesh::BulkData bulk(meta, pm);

    stk::mesh::Part* ssPart = stk::unit_test_util::create_AB_mesh_with_sideset(bulk, direction);
    stk::mesh::Entity  face = verify_and_get_face(bulk, ssPart);

    std::pair<bool,bool> sidesetPolarity = stk::io::is_positive_sideset_polarity(bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first);
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second);

    std::string file = name+".e";
    stk::io::write_mesh(file, bulk);

    sidesetPolarity = stk::io::is_positive_sideset_polarity(bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first);
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second);
    unlink(file.c_str());
}

TEST(StkIo, sideset_polarity_AA)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);

  if(p_size == 1)
  {
      check_polarity_of_modified_AA(pm, stk::unit_test_util::RIGHT, "ARA", false);
      check_polarity_of_modified_AA(pm, stk::unit_test_util::LEFT , "ALA", true);
  }
}

TEST(StkIo, sideset_polarity_AB)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  const unsigned p_size = stk::parallel_machine_size(pm);

  if(p_size == 1)
  {
      check_polarity_of_modified_AB(pm, stk::unit_test_util::RIGHT, "ARB", false);
      check_polarity_of_modified_AB(pm, stk::unit_test_util::LEFT , "ALB", true);
  }
}

}
