#include "gtest/gtest.h"
#include "mpi.h"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphImpl.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "stk_unit_test_utils/TextMesh.hpp"
#include <stk_mesh/baseImpl/elementGraph/GraphEdgeData.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>
#include <stk_mesh/baseImpl/elementGraph/ParallelInfoForGraph.hpp>
#include <stk_io/IossBridge.hpp>

class ElemElemGraphTester : public stk::mesh::ElemElemGraph
{
public:
  ElemElemGraphTester(stk::mesh::BulkData& bulkData) :
    stk::mesh::ElemElemGraph(bulkData) { }

  stk::mesh::impl::ParallelGraphInfo& get_parallel_info() { return m_parallelInfoForGraphEdges.get_parallel_graph_info(); }
};

class ParallelGraphUpdate : public stk::unit_test_util::MeshFixture
{
public:
  ParallelGraphUpdate() : stk::unit_test_util::MeshFixture()
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    airPart = &get_meta().declare_part_with_topology("air", stk::topology::HEXAHEDRON_8);
    air = *airPart;

    skinPart = &get_meta().declare_part_with_topology("things to skin", stk::topology::HEXAHEDRON_8);
    skinSelector = *skinPart;

    make_part_an_element_block_by_setting_id();
  }

  void put_some_elements_into_some_parts()
  {
    stk::mesh::EntityVector elements = get_half_of_elements();
    add_part_to_elements(elements,airPart);
    add_part_to_elements(elements,skinPart);
  }

  void add_part_to_elements(const stk::mesh::EntityVector &elements, stk::mesh::Part* part)
  {
    get_bulk().modification_begin();
    for(stk::mesh::Entity element : elements )
      get_bulk().change_entity_parts(element, stk::mesh::ConstPartVector{part});
    get_bulk().modification_end();
  }

  void remove_part_from_elements(const stk::mesh::EntityVector &elements, stk::mesh::Part* part)
  {
    get_bulk().modification_begin();
    for(stk::mesh::Entity element : elements )
      get_bulk().change_entity_parts(element, stk::mesh::ConstPartVector{}, stk::mesh::ConstPartVector{part});
    get_bulk().modification_end();
  }

  stk::mesh::Selector get_air_selector() const
  {
    return air;
  }

  stk::mesh::Selector get_skin_selector() const
  {
    return skinSelector;
  }

  void move_elements_out_of_parts()
  {
    remove_part_from_elements(get_half_of_elements(), airPart);
    remove_part_from_elements(get_half_of_elements(), skinPart);
  }

  void move_elements_into_parts_and_verify()
  {
    expect_no_elements_selected(get_air_selector());
    expect_no_elements_selected(get_skin_selector());
    put_some_elements_into_some_parts();
    expect_some_elements_selected(get_air_selector());
    expect_some_elements_selected(get_skin_selector());
  }

  void verify_selectors_not_in_parallel_graph(const stk::mesh::impl::ParallelGraphInfo& parallel_info)
  {
    for(const auto& item : parallel_info)
    {
      stk::mesh::impl::LocalId elemOnOtherProc = item.first.elem2();
      EXPECT_FALSE(remoteAirSelector[elemOnOtherProc]);
      EXPECT_FALSE(remoteSkinSelector[elemOnOtherProc]);
    }
  }

  void verify_part_ordinals_not_in_parallel_graph(const stk::mesh::impl::ParallelPartInfo &parallelPartInfo) const
  {
    for(const auto& item : parallelPartInfo)
    {
      const stk::mesh::impl::PartOrdinals &partOrdinals = item.second.elementPartOrdinals;
      EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), airPart->mesh_meta_data_ordinal()));
      EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), skinPart->mesh_meta_data_ordinal()));
    }
  }

  void verify_selectors_in_parallel_graph(const stk::mesh::impl::ParallelGraphInfo& parallel_info)
  {
    for(const auto& item : parallel_info)
    {
      stk::mesh::impl::LocalId elemOnOtherProc = item.first.elem2();
      if(elemOnOtherProc %2 == 0)
      {
        EXPECT_TRUE(remoteAirSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
        EXPECT_TRUE(remoteSkinSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
      }
      else
      {
        EXPECT_FALSE(remoteAirSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
        EXPECT_FALSE(remoteSkinSelector[elemOnOtherProc]) << " for element " << -elemOnOtherProc;
      }
    }
  }

  void verify_part_ordinals_are_in_parallel_graph(const stk::mesh::impl::ParallelPartInfo &parallelPartInfo) const
  {
    for(const auto& item : parallelPartInfo)
    {
      const stk::mesh::impl::LocalId elemOnOtherProc = item.first;
      const stk::mesh::impl::PartOrdinals &partOrdinals = item.second.elementPartOrdinals;

      if(elemOnOtherProc%2 == 0)
      {
        EXPECT_TRUE(std::binary_search(partOrdinals.begin(), partOrdinals.end(), airPart->mesh_meta_data_ordinal()));
        EXPECT_TRUE(std::binary_search(partOrdinals.begin(), partOrdinals.end(), skinPart->mesh_meta_data_ordinal()));
      }
      else
      {
        EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), airPart->mesh_meta_data_ordinal()));
        EXPECT_TRUE(!std::binary_search(partOrdinals.begin(), partOrdinals.end(), skinPart->mesh_meta_data_ordinal()));
      }
    }
  }


  void make_part_an_element_block_by_setting_id()
  {
    get_meta().set_part_id(*airPart, 100);
    get_meta().set_part_id(*skinPart, 101);
  }

  stk::mesh::EntityVector get_half_of_elements()
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), elements);
    return get_elements_that_have_even_identifiers(elements);
  }

  stk::mesh::EntityVector get_elements_that_have_even_identifiers(const stk::mesh::EntityVector& elements)
  {
    stk::mesh::EntityVector someElements;
    for(stk::mesh::Entity element : elements)
      if(get_bulk().identifier(element)%2 == 0)
        someElements.push_back(element);
    return someElements;
  }

  void expect_no_elements_selected(stk::mesh::Selector selector)
  {
    size_t num_elements = stk::mesh::count_selected_entities(selector, get_bulk().buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(0u, num_elements);
  }

  void expect_some_elements_selected(stk::mesh::Selector selector)
  {
    size_t num_elements= stk::mesh::count_selected_entities(selector, get_bulk().buckets(stk::topology::ELEM_RANK));
    EXPECT_TRUE(num_elements != 0u);
  }

  stk::mesh::Part* airPart = nullptr;
  stk::mesh::Part* skinPart = nullptr;
  stk::mesh::Selector air;
  stk::mesh::Selector skinSelector;
  stk::mesh::impl::ParallelSelectedInfo remoteSkinSelector;
  stk::mesh::impl::ParallelSelectedInfo remoteAirSelector;
};

// Need updater for part ordinals

TEST_F(ParallelGraphUpdate, updateAirSelector)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD)==4)
  {
    stk::io::fill_mesh("generated:4x4x4", get_bulk());
    ElemElemGraphTester graph(get_bulk());

    stk::mesh::impl::ParallelPartInfo parallelPartInfo;
    stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
    stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_skin_selector(), remoteSkinSelector);

    verify_selectors_not_in_parallel_graph(graph.get_parallel_info());
    verify_part_ordinals_not_in_parallel_graph(parallelPartInfo);

    move_elements_into_parts_and_verify();
    stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_air_selector(), remoteAirSelector);
    stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_skin_selector(), remoteSkinSelector);
    stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
    verify_selectors_in_parallel_graph(graph.get_parallel_info());
    verify_part_ordinals_are_in_parallel_graph(parallelPartInfo);

    move_elements_out_of_parts();
    stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_air_selector(), remoteAirSelector);
    stk::mesh::impl::populate_selected_value_for_remote_elements(get_bulk(), graph, get_skin_selector(), remoteSkinSelector);
    stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo);
    verify_selectors_not_in_parallel_graph(graph.get_parallel_info());
    verify_part_ordinals_not_in_parallel_graph(parallelPartInfo);
  }
}

TEST_F(ParallelGraphUpdate, deleteBothSolidElementsOnParallelEdge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD)==2)
  {
    stk::io::fill_mesh("generated:1x1x4", get_bulk());
    get_bulk().initialize_face_adjacent_element_graph();

    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1u);
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2u);
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3u);
    stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4u);

    if(get_bulk().parallel_rank() == 0) {
      ASSERT_TRUE(get_bulk().is_valid(elem1));
      ASSERT_TRUE(get_bulk().is_valid(elem2));
      ASSERT_FALSE(get_bulk().is_valid(elem3));
      ASSERT_FALSE(get_bulk().is_valid(elem4));
    }

    if(get_bulk().parallel_rank() == 1) {
      ASSERT_FALSE(get_bulk().is_valid(elem1));
      ASSERT_FALSE(get_bulk().is_valid(elem2));
      ASSERT_TRUE(get_bulk().is_valid(elem3));
      ASSERT_TRUE(get_bulk().is_valid(elem4));
    }

    const stk::mesh::ElemElemGraph &graph = get_bulk().get_face_adjacent_element_graph();
    const stk::mesh::impl::ParallelGraphInfo& parallelGraphInfo = graph.get_parallel_info_for_graph_edges().get_parallel_graph_info();

    EXPECT_EQ(1u, parallelGraphInfo.size());
    const stk::mesh::GraphEdge& graphEdge = parallelGraphInfo.begin()->first;

    if(get_bulk().parallel_rank() == 0) {
      stk::mesh::impl::LocalId localId = graphEdge.elem1();
      stk::mesh::Entity localElem = graph.get_entity(localId);
      EXPECT_EQ(2u, get_bulk().identifier(localElem));
      EXPECT_EQ(5, graphEdge.side1());
      EXPECT_EQ(-3, graphEdge.elem2());
      EXPECT_EQ(4, graphEdge.side2());
    }
    if(get_bulk().parallel_rank() == 1) {
      stk::mesh::impl::LocalId localId = graphEdge.elem1();
      stk::mesh::Entity localElem = graph.get_entity(localId);
      EXPECT_EQ(3u, get_bulk().identifier(localElem));
      EXPECT_EQ(4, graphEdge.side1());
      EXPECT_EQ(-2, graphEdge.elem2());
      EXPECT_EQ(5, graphEdge.side2());
    }

    get_bulk().modification_begin();
    get_bulk().destroy_entity(elem2);
    get_bulk().destroy_entity(elem3);
    get_bulk().modification_end();

    EXPECT_EQ(0u, parallelGraphInfo.size());

    stk::mesh::impl::ParallelPartInfo parallelPartInfo;
    EXPECT_NO_THROW(stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo));
  }
}

TEST_F(ParallelGraphUpdate, deleteBothShellElementsOnParallelEdge)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD)==2)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8,block_2\n"
                           "1,3,SHELL_QUAD_4,5,6,7,8,block_3\n"
                           "1,4,HEX_8,5,6,7,8,9,10,11,12,block_1";

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 1,1,0, 0,1,0,
      0,0,1, 1,0,1, 1,1,1, 0,1,1,
      0,0,2, 1,0,2, 1,1,2, 0,1,2
    };

    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
    get_bulk().initialize_face_adjacent_element_graph();

    stk::mesh::Entity elem1 = get_bulk().get_entity(stk::topology::ELEM_RANK, 1u);
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2u);
    stk::mesh::Entity elem3 = get_bulk().get_entity(stk::topology::ELEM_RANK, 3u);
    stk::mesh::Entity elem4 = get_bulk().get_entity(stk::topology::ELEM_RANK, 4u);

    if(get_bulk().parallel_rank() == 0) {
      ASSERT_TRUE(get_bulk().is_valid(elem1));
      ASSERT_TRUE(get_bulk().is_valid(elem2));
      ASSERT_FALSE(get_bulk().is_valid(elem3));
      ASSERT_FALSE(get_bulk().is_valid(elem4));
    }

    if(get_bulk().parallel_rank() == 1) {
      ASSERT_FALSE(get_bulk().is_valid(elem1));
      ASSERT_FALSE(get_bulk().is_valid(elem2));
      ASSERT_TRUE(get_bulk().is_valid(elem3));
      ASSERT_TRUE(get_bulk().is_valid(elem4));
    }

    const stk::mesh::ElemElemGraph &graph = get_bulk().get_face_adjacent_element_graph();
    const stk::mesh::impl::ParallelGraphInfo& parallelGraphInfo = graph.get_parallel_info_for_graph_edges().get_parallel_graph_info();

    EXPECT_EQ(4u, parallelGraphInfo.size());

    get_bulk().modification_begin();
    get_bulk().destroy_entity(elem2);
    get_bulk().destroy_entity(elem3);
    get_bulk().modification_end();

    EXPECT_EQ(1u, parallelGraphInfo.size());

    stk::mesh::impl::ParallelPartInfo parallelPartInfo;
    EXPECT_NO_THROW(stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo));
  }
}

TEST_F(ParallelGraphUpdate, createAefA_FromScratch)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD)==2)
  {
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "0,2,SHELL_QUAD_4,5,6,7,8,block_2\n"
                           "1,3,SHELL_QUAD_4,5,6,7,8,block_3\n"
                           "1,4,HEX_8,5,6,7,8,9,10,11,12,block_1";

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 1,1,0, 0,1,0,
      0,0,1, 1,0,1, 1,1,1, 0,1,1,
      0,0,2, 1,0,2, 1,1,2, 0,1,2
    };

    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
    get_bulk().initialize_face_adjacent_element_graph();

    const stk::mesh::ElemElemGraph &graph = get_bulk().get_face_adjacent_element_graph();
    const stk::mesh::impl::ParallelGraphInfo& parallelGraphInfo = graph.get_parallel_info_for_graph_edges().get_parallel_graph_info();

    EXPECT_EQ(4u, parallelGraphInfo.size());

    stk::mesh::impl::ParallelPartInfo parallelPartInfo;
    EXPECT_NO_THROW(stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo));
  }
}

TEST_F(ParallelGraphUpdate, createAefA_FromAA)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD)==2)
  {
    stk::mesh::Part& block_2 = get_meta().declare_part_with_topology("block_2", stk::topology::SHELL_QUAD_4);
    stk::mesh::Part& block_3 = get_meta().declare_part_with_topology("block_3", stk::topology::SHELL_QUAD_4);

    stk::mesh::PartVector shellParts{&block_2, &block_3};

    stk::io::put_io_part_attribute(block_2);
    stk::io::put_io_part_attribute(block_3);

    get_meta().set_part_id(block_2, 2u);
    get_meta().set_part_id(block_3, 3u);

    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                           "1,4,HEX_8,5,6,7,8,9,10,11,12,block_1";

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 1,1,0, 0,1,0,
      0,0,1, 1,0,1, 1,1,1, 0,1,1,
      0,0,2, 1,0,2, 1,1,2, 0,1,2
    };

    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
    get_bulk().initialize_face_adjacent_element_graph();

    const stk::mesh::ElemElemGraph &graph = get_bulk().get_face_adjacent_element_graph();
    const stk::mesh::impl::ParallelGraphInfo& parallelGraphInfo = graph.get_parallel_info_for_graph_edges().get_parallel_graph_info();

    EXPECT_EQ(1u, parallelGraphInfo.size());

    stk::mesh::EntityIdVector shellNodeIds;
    for(unsigned i=5; i<=8; i++) {
      stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, i);
      ASSERT_TRUE(get_bulk().is_valid(node));
      shellNodeIds.push_back(i);
    }

    get_bulk().modification_begin();
    stk::mesh::EntityId shellId = get_bulk().parallel_rank() + 2;
    stk::mesh::Part* shellPart = shellParts[get_bulk().parallel_rank()];
    stk::mesh::Entity elem = stk::mesh::declare_element(get_bulk(), *shellPart, shellId, shellNodeIds);
    get_bulk().modification_end();

    EXPECT_EQ(4u, parallelGraphInfo.size());

    get_bulk().modification_begin();
    get_bulk().destroy_entity(elem);
    get_bulk().modification_end();

    EXPECT_EQ(1u, parallelGraphInfo.size());

    stk::mesh::impl::ParallelPartInfo parallelPartInfo;
    EXPECT_NO_THROW(stk::mesh::impl::populate_part_ordinals_for_remote_edges(get_bulk(), graph, parallelPartInfo));
  }
}
