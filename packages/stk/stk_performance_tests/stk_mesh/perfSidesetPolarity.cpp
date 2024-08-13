#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>

#include <stk_unit_test_utils/TextMesh.hpp>
#include "stk_mesh/base/Comm.hpp"
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SideSetUtil.hpp>
#include <stk_mesh/base/SidesetUpdater.hpp>
#include <stk_mesh/base/PolarityUtil.hpp>

#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp"
#include "stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp"

namespace
{

//==============================================================================
class PolarityPerformanceBase : public stk::unit_test_util::PerformanceTester
{
public:
  PolarityPerformanceBase(stk::mesh::BulkData &bulk)
    : stk::unit_test_util::PerformanceTester(bulk.parallel()),
      m_bulk(bulk)
  {
    create_mesh();
  }

  virtual ~PolarityPerformanceBase() = default;

protected:

  virtual size_t get_value_to_output_as_iteration_count() override
  {
    return 1;
  }

  std::string get_mesh_spec()
  {
    std::string meshDesc;

    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2) {
      meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
                 "1,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
    } else {
      meshDesc = "0,1,HEX_8,1,2,3,4,5, 6, 7, 8,block_1\n"
                 "0,2,HEX_8,5,6,7,8,9,10,11,12,block_2";
    }

    return meshDesc;
  }

  void create_mesh()
  {
    stk::mesh::MetaData &meta = m_bulk.mesh_meta_data();

    // Register sideset updater
    if (!m_bulk.has_observer_type<stk::mesh::SidesetUpdater>()) {
      stk::mesh::Selector activeSelector = meta.universal_part();

      m_bulk.register_observer(std::make_shared<stk::mesh::IncrementalSidesetUpdater>(m_bulk, activeSelector));
    }

    // Create sideset parts
    std::vector<std::string> sidesetNames{"surface_1", "surface_block_1_QUAD4_1"};
    int sidesetId = 1;
    for(const std::string &name : sidesetNames)
    {
      stk::mesh::Part &part = meta.declare_part_with_topology(name, stk::topology::QUAD_4);
      meta.set_part_id(part, sidesetId);
      m_sidesetParts.push_back(&part);
    }

    meta.declare_part_subset(*m_sidesetParts[0], *m_sidesetParts[1]);

    // Build mesh
    std::string meshDesc = get_mesh_spec();
    std::vector<double> coordinates = { 0,0,0, 1,0,0, 1,1,0, 0,1,0,
                                        0,0,1, 1,0,1, 1,1,1, 0,1,1,
                                        0,0,2, 1,0,2, 1,1,2, 0,1,2 };

    m_bulk.initialize_face_adjacent_element_graph();
    stk::unit_test_util::setup_text_mesh(m_bulk, stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));

    // Create and populate sideset between elements 1 & 2, based on element 1
    stk::mesh::SideSet &sideset = m_bulk.create_sideset(*m_sidesetParts[0]);
    sideset.set_accept_all_internal_non_coincident_entries(false);
    m_sideset = &sideset;

    std::vector<stk::mesh::EntityId> elementIds{1, 2};
    std::vector<stk::mesh::ConnectivityOrdinal> elementOrdinals{5, 4};
    std::vector<stk::mesh::Entity> elements(2);
    elements[0] = m_bulk.get_entity(stk::topology::ELEMENT_RANK, elementIds[0]);
    elements[1] = m_bulk.get_entity(stk::topology::ELEMENT_RANK, elementIds[1]);

    if (m_bulk.is_valid(elements[0])) {
      sideset.add({elements[0], elementOrdinals[0]});
    }

    // Setup surface to block mapping for the sideset updater
    stk::mesh::Part* block_1 = meta.get_part("block_1");
    EXPECT_TRUE(block_1 != nullptr);

    std::vector<const stk::mesh::Part*> touchingParts{block_1};
    meta.set_surface_to_block_mapping(m_sidesetParts[1], touchingParts);

    // Create faces
    m_bulk.modification_begin();
    int myProc = m_bulk.parallel_rank();
    m_face = m_bulk.declare_element_side(elements[myProc], elementOrdinals[myProc], m_sidesetParts);
    m_bulk.modification_end();

    // Make sure there is only 1 global face
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(m_bulk, counts);
    EXPECT_EQ(1u, counts[meta.side_rank()]);

    // Make sure it matches the one we created
    stk::mesh::EntityVector faces;
    stk::mesh::get_selected_entities(meta.universal_part(), m_bulk.buckets(meta.side_rank()), faces);
    EXPECT_EQ(1u, faces.size());
    EXPECT_EQ(m_face, faces[0]);
  }

  int get_num_loops()
  {
    return stk::unit_test_util::get_command_line_option("-l", m_defaultLoops);
  }

  stk::mesh::BulkData &m_bulk;
  stk::mesh::PartVector m_sidesetParts;
  stk::mesh::SideSet *m_sideset{nullptr};
  stk::mesh::Entity m_face;

  const int m_defaultLoops = 100000;
};

class PolarityUtilPerformance : public PolarityPerformanceBase
{
public:
  PolarityUtilPerformance(stk::mesh::BulkData &bulk)
    : PolarityPerformanceBase(bulk),
      m_helper(bulk, bulk.mesh_meta_data().universal_part())
  {
  }

  ~PolarityUtilPerformance() = default;

protected:
  virtual void run_algorithm_to_time() override
  {
    // Get polarity results with util wrapper code
    std::pair<bool,bool> expectedPolarity(true, true);
    int numLoops = get_num_loops();

    for(int i=0; i<numLoops; i++) {
      std::pair<bool,bool> polarity1 = m_helper.is_positive_sideset_polarity(*m_sidesetParts[1], m_face, m_sideset);
      std::pair<bool,bool> polarity2 = m_helper.is_positive_sideset_face_polarity(m_face);


      EXPECT_EQ(expectedPolarity, polarity1);
      EXPECT_EQ(expectedPolarity, polarity2);
    }
  }

  stk::mesh::PolarityUtil m_helper;
};

class PolarityFreePerformance : public PolarityPerformanceBase
{
public:
  PolarityFreePerformance(stk::mesh::BulkData &bulk)
    : PolarityPerformanceBase(bulk)
  {
  }

  ~PolarityFreePerformance() = default;

protected:
  virtual void run_algorithm_to_time() override
  {
    // Get polarity results with util wrapper code
    std::pair<bool,bool> expectedPolarity(true, true);
    int numLoops = get_num_loops();
    stk::mesh::Part& activePart = m_bulk.mesh_meta_data().universal_part();

    for(int i=0; i<numLoops; i++) {
      // Get polarity results with free function code
      std::pair<bool,bool> polarity1 = stk::mesh::is_positive_sideset_polarity(m_bulk, *m_sidesetParts[1], m_face,
                                                                               &activePart, m_sideset);
      std::pair<bool,bool> polarity2 = stk::mesh::is_positive_sideset_face_polarity(m_bulk, m_face, &activePart);

      EXPECT_EQ(expectedPolarity, polarity1);
      EXPECT_EQ(expectedPolarity, polarity2);
    }
  }
};

class PolarityPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
  void run_util_polarity_perf_test()
  {
    PolarityUtilPerformance perfTester(get_bulk());
    perfTester.run_performance_test();
  }

  void run_free_polarity_perf_test()
  {
    PolarityFreePerformance perfTester(get_bulk());
    perfTester.run_performance_test();
  }
};

TEST_F(PolarityPerformanceTest, utilPolarityTestWithAura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {GTEST_SKIP();}
  allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
  run_util_polarity_perf_test();
}

TEST_F(PolarityPerformanceTest, utilPolarityTestNoAura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {GTEST_SKIP();}
  allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
  run_util_polarity_perf_test();
}

TEST_F(PolarityPerformanceTest, freePolarityTestWithAura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {GTEST_SKIP();}
  allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
  run_free_polarity_perf_test();
}

TEST_F(PolarityPerformanceTest, freePolarityTestNoAura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) > 2) {GTEST_SKIP();}
  allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
  run_free_polarity_perf_test();
}

}




