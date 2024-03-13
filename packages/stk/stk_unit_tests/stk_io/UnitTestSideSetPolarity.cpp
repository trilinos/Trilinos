#include <gtest/gtest.h>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SidesetUpdater.hpp>
#include <stk_io/StkIoUtils.hpp>
#include "stk_unit_test_utils/ReadWriteSidesetTester.hpp"
#include "stk_unit_test_utils/FaceTestingUtils.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"
#include "stk_unit_test_utils/GenerateALefRAMesh.hpp"
#include <stk_unit_test_utils/BuildMesh.hpp>

namespace {

using stk::unit_test_util::build_mesh;
using stk::unit_test_util::build_mesh_no_simple_fields;

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
    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, pm);

    if (!bulk->has_observer_type<stk::mesh::SidesetUpdater>()) {
      stk::mesh::Selector activeSelector = bulk->mesh_meta_data().universal_part();

      bulk->register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(*bulk, activeSelector));
    }

    bulk->initialize_face_adjacent_element_graph();

    stk::mesh::Part* ssPart = stk::unit_test_util::create_AA_mesh_with_sideset(*bulk, direction);
    stk::mesh::Entity  face = verify_and_get_face(*bulk, ssPart);

    std::pair<bool,bool> sidesetPolarity = stk::mesh::is_positive_sideset_polarity(*bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first) << "Failure for " << name << " on first";
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second) << "Failure for " << name << " on second";
}

void check_polarity_of_modified_AB(stk::ParallelMachine pm,
                                   stk::unit_test_util::SidesetDirection direction,
                                   const std::string& name,
                                   bool isPositivePolarity)
{
    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, pm);

    if (!bulk->has_observer_type<stk::mesh::SidesetUpdater>()) {
      stk::mesh::Selector activeSelector = bulk->mesh_meta_data().universal_part();

      bulk->register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(*bulk, activeSelector));
    }

    bulk->initialize_face_adjacent_element_graph();

    stk::mesh::Part* ssPart = stk::unit_test_util::create_AB_mesh_with_sideset(*bulk, direction);
    stk::mesh::Entity  face = verify_and_get_face(*bulk, ssPart);

    std::pair<bool,bool> sidesetPolarity = stk::mesh::is_positive_sideset_polarity(*bulk, *ssPart, face);
    EXPECT_TRUE(sidesetPolarity.first) << "Failure for " << name << " on first";
    EXPECT_EQ(isPositivePolarity, sidesetPolarity.second) << "Failure for " << name << " on second";
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

void update_sidesets(const stk::mesh::BulkData& bulk, std::ostream& os)
{
    std::vector<std::shared_ptr<stk::mesh::SidesetUpdater> > updaters = bulk.get_observer_type<stk::mesh::SidesetUpdater>();
    STK_ThrowRequireMsg(!updaters.empty(), "ERROR, no SidesetUpdater found on stk::mesh::BulkData");

    updaters[0]->set_output_stream(os);

    std::vector<size_t> values;
    updaters[0]->fill_values_to_reduce(values);
    std::vector<size_t> maxValues(values);

    if (stk::parallel_machine_size(bulk.parallel()) > 1) {
        stk::all_reduce_max(bulk.parallel(), values.data(), maxValues.data(), maxValues.size());
    }

    updaters[0]->set_reduced_values(maxValues);
}

TEST(StkIo, check_internal_sideset_warning_with_reconstruction)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(pm) != 1) {
      return;
    }

    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::unit_test_util::SidesetDirection direction = stk::unit_test_util::LEFT;
    stk::mesh::Part* ssPart = stk::unit_test_util::create_AA_mesh_with_sideset(*bulk, direction);
    EXPECT_TRUE(ssPart != nullptr);

    stk::mesh::Selector activeSelector(meta.universal_part());
    bulk->register_observer(std::make_shared<stk::mesh::ReconstructionSidesetUpdater>(*bulk, activeSelector));
    bulk->modification_end();

    std::vector<std::shared_ptr<stk::mesh::SidesetUpdater> > updaters = bulk->get_observer_type<stk::mesh::SidesetUpdater>();
    STK_ThrowRequireMsg(!updaters.empty(), "ERROR, no SidesetUpdater found on stk::mesh::BulkData");
    updaters[0]->set_warn_about_internal_sideset(true);

    {//mesh was just modified, including elements touching the sideset, so expect warning
      std::ostringstream os;
      update_sidesets(*bulk, os);

      std::string warningString = os.str();

      std::string expected = "WARNING, Internal sideset ("+ ssPart->name() +") detected. STK doesn't support internal sidesets\n(i.e., sidesets between elements where both elements are in the same block)\nExecution will continue but correct results are not guaranteed. Contact sierra-help@sandia.gov\n";

      EXPECT_EQ(expected, warningString);
    }

    bulk->modification_begin();
    bulk->modification_end();

    {//mesh now not modified, so expect no warning about internal sideset
      std::ostringstream os;
      update_sidesets(*bulk, os);

      std::string warningString = os.str();
      std::string expected = "";
      EXPECT_EQ(expected, warningString);
    }

    bulk->modification_begin();

    //add an element that doesn't touch the sideset:
    stk::mesh::Part& elemTopoPart = meta.get_topology_root_part(stk::topology::HEX_8);
    stk::mesh::EntityId elemId = 3;
    stk::mesh::EntityIdVector nodeIds = {9, 10, 11, 12, 13, 14, 15, 16};
    stk::mesh::declare_element(*bulk, elemTopoPart, elemId, nodeIds);

    bulk->modification_end();

    {//mesh modified, but not sideset, so still no warning about internal sideset
      std::ostringstream os;
      update_sidesets(*bulk, os);

      std::string warningString = os.str();
      std::string expected = "";
      EXPECT_EQ(expected, warningString);
    }
}

TEST(StkIo, check_internal_sideset_warning_with_incremental_update)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(pm) != 1) {
      return;
    }

    const unsigned spatialDim = 3;
    std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(spatialDim, pm);
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::mesh::Selector activeSelector(meta.universal_part());
    bulk->register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(*bulk, activeSelector));

    std::vector<std::shared_ptr<stk::mesh::SidesetUpdater> > updaters = bulk->get_observer_type<stk::mesh::SidesetUpdater>();
    STK_ThrowRequireMsg(!updaters.empty(), "ERROR, no SidesetUpdater found on stk::mesh::BulkData");
    std::ostringstream os;
    updaters[0]->set_output_stream(os);
    updaters[0]->set_warn_about_internal_sideset(true);

    stk::unit_test_util::SidesetDirection direction = stk::unit_test_util::LEFT;
    stk::mesh::Part* ssPart = stk::unit_test_util::create_AA_mesh_with_sideset(*bulk, direction);
    EXPECT_TRUE(ssPart != nullptr);

    bulk->modification_end();

    {//mesh was just modified, including elements touching the sideset, so expect warning
      std::string warningString = os.str();

      std::string expected = "WARNING, Internal sideset ("+ ssPart->name() +") detected. STK doesn't support internal sidesets\n(i.e., sidesets between elements where both elements are in the same block)\nExecution will continue but correct results are not guaranteed. Contact sierra-help@sandia.gov\n";

      EXPECT_EQ(expected, warningString);
    }

    os.str("");
    os.clear();
    bulk->modification_begin();
    bulk->modification_end();

    {//mesh now not modified, so expect no warning about internal sideset
      std::string warningString = os.str();
      std::string expected = "";
      EXPECT_EQ(expected, warningString);
    }

    os.str("");
    os.clear();
    bulk->modification_begin();

    //add an element that doesn't touch the sideset:
    stk::mesh::Part& elemTopoPart = meta.get_topology_root_part(stk::topology::HEX_8);
    stk::mesh::EntityId elemId = 3;
    stk::mesh::EntityIdVector nodeIds = {9, 10, 11, 12, 13, 14, 15, 16};
    stk::mesh::declare_element(*bulk, elemTopoPart, elemId, nodeIds);

    bulk->modification_end();

    {//mesh modified, but not sideset, so still no warning about internal sideset
      std::string warningString = os.str();
      std::string expected = "";
      EXPECT_EQ(expected, warningString);
    }
}

}
