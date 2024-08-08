#include "stk_balance/balance.hpp"
#include "stk_balance/internal/DetectAndFixMechanisms.hpp"
#include "stk_balance/internal/Zoltan2ParallelGraph.hpp"
#include "stk_balance/internal/privateDeclarations.hpp"
#include "stk_mesh/base/DestroyElements.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp"
#include "stk_tools/mesh_tools/CustomAura.hpp"
#include "stk_unit_test_utils/getOption.h"
#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace
{

class MechanismMesh2x2 : public stk::unit_test_util::MeshFixture
{
protected:
  void setup_mechanistic_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x2x2", auraOption);
    std::vector<std::vector<stk::mesh::EntityId>> decomp = {
      {1u, 4u},
      {2u, 3u}
    };

    std::vector<stk::mesh::EntityId> myElements = decomp[get_bulk().parallel_rank()];

    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), elements);

    stk::mesh::EntityProcVec elementMoves;

    int otherProc = 1 - get_bulk().parallel_rank();

    for(stk::mesh::Entity element : elements)
    {
      stk::mesh::EntityId id = get_bulk().identifier(element);
      if(std::find(myElements.begin(), myElements.end(), id) == myElements.end())
        elementMoves.push_back(std::make_pair(element, otherProc));
    }

    get_bulk().change_entity_owner(elementMoves);
  }
};

TEST_F(MechanismMesh2x2, detection_with_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_mechanistic_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::io::write_mesh("before_junk.g", get_bulk());
    stk::balance::GraphCreationSettings graphSettings;
    if(graphSettings.shouldFixMechanisms())
    {
      EXPECT_TRUE(stk::balance::internal::detectAndFixMechanisms(graphSettings, get_bulk()));
    }
    stk::io::write_mesh("after_junk.g", get_bulk());
  }
}

TEST_F(MechanismMesh2x2, detection_without_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_mechanistic_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::io::write_mesh("junk.g", get_bulk());
    stk::balance::GraphCreationSettings graphSettings;
    if(graphSettings.shouldFixMechanisms())
    {
      EXPECT_TRUE(stk::balance::internal::detectAndFixMechanisms(graphSettings, get_bulk()));
    }
  }
}

TEST_F(MechanismMesh2x2, move_components)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_mechanistic_mesh(stk::mesh::BulkData::AUTO_AURA);

    std::vector<stk::mesh::EntityVector> elementsToMove;
    elementsToMove.resize(1);

    stk::mesh::Entity element;
    if(get_bulk().parallel_rank()==0)
      element = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);
    else
      element = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);

    ASSERT_TRUE(get_bulk().is_valid(element) && get_bulk().bucket(element).owned());
    elementsToMove[0].push_back(element);

    stk::balance::GraphCreationSettings graphSettings;
    stk::mesh::impl::LocalIdMapper localIds(get_bulk(), stk::topology::ELEM_RANK);

    stk::tools::create_custom_aura(get_bulk(), get_bulk().mesh_meta_data().globally_shared_part(), "customAura");

    Zoltan2ParallelGraph zoltan2Graph;
    stk::balance::internal::fill_zoltan2_parallel_graph(get_bulk(), graphSettings, zoltan2Graph);

    std::vector<int> componentsToMove = {0};
    stk::balance::internal::move_components(zoltan2Graph, localIds, get_bulk(), elementsToMove, componentsToMove);

    if(get_bulk().parallel_rank()==0)
      element = get_bulk().get_entity(stk::topology::ELEM_RANK, 3);
    else
      element = get_bulk().get_entity(stk::topology::ELEM_RANK, 1);

    EXPECT_TRUE(get_bulk().is_valid(element) && get_bulk().bucket(element).owned());

    stk::io::write_mesh("junk.g", get_bulk());
  }
}


class LotsOfComponentsMesh : public stk::unit_test_util::MeshFixture
{
protected:
  void setup_mechanistic_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x10x10", auraOption);
    std::vector<std::vector<stk::mesh::EntityId>> decomp = {
      {2u, 3u, 12u, 13u, 24u, 25u, 34u, 35u, 46u, 27u, 28u, 37u, 38u, 9u, 10u, 19u, 20u},
      {}
    };

    std::vector<stk::mesh::EntityId> myElements = decomp[get_bulk().parallel_rank()];

    stk::mesh::EntityProcVec elementMoves;

    int otherProc = 1 - get_bulk().parallel_rank();

    for(stk::mesh::EntityId elementId : myElements)
    {
      stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elementId);
      EXPECT_TRUE(get_bulk().is_valid(element));
      elementMoves.push_back(std::make_pair(element, otherProc));
    }

    get_bulk().change_entity_owner(elementMoves);
  }
};

TEST_F(LotsOfComponentsMesh, detection_with_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_mechanistic_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::io::write_mesh("before_junk.g", get_bulk());
    stk::balance::GraphCreationSettings graphSettings;
    if(graphSettings.shouldFixMechanisms())
    {
      EXPECT_TRUE(stk::balance::internal::detectAndFixMechanisms(graphSettings, get_bulk()));
    }
    stk::io::write_mesh("after_junk.g", get_bulk());
  }
}

TEST_F(LotsOfComponentsMesh, detection_without_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_mechanistic_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::io::write_mesh("junk.g", get_bulk());
    stk::balance::GraphCreationSettings graphSettings;
    if(graphSettings.shouldFixMechanisms())
    {
      EXPECT_TRUE(stk::balance::internal::detectAndFixMechanisms(graphSettings, get_bulk()));
    }
  }
}


class GlobalMeshWithMechanism : public stk::unit_test_util::MeshFixture
{
protected:
  void setup_mechanistic_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_mesh("generated:1x3x3", auraOption);
    std::vector<stk::mesh::EntityId> deleteThese = {
      {3u, 6u, 7u, 8u}
    };

    stk::mesh::EntityVector localElementsToDelete;

    for(stk::mesh::EntityId elementId : deleteThese)
    {
      stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, elementId);
      if(get_bulk().is_valid(element) && get_bulk().bucket(element).owned())
        localElementsToDelete.push_back(element);
    }

    stk::mesh::destroy_elements(get_bulk(), localElementsToDelete);

    stk::balance::GraphCreationSettings graphSettings;
    stk::balance::balanceStkMesh(graphSettings, get_bulk());
  }
};

TEST_F(GlobalMeshWithMechanism, detection_with_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_mechanistic_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::io::write_mesh("junk.g", get_bulk());
    stk::balance::GraphCreationSettings graphSettings;
    if(graphSettings.shouldFixMechanisms())
    {
      EXPECT_TRUE(stk::balance::internal::detectAndFixMechanisms(graphSettings, get_bulk()));
    }
  }
}

TEST_F(GlobalMeshWithMechanism, detection_without_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    setup_mechanistic_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::io::write_mesh("junk.g", get_bulk());
    stk::balance::GraphCreationSettings graphSettings;
    if(graphSettings.shouldFixMechanisms())
    {
      EXPECT_TRUE(stk::balance::internal::detectAndFixMechanisms(graphSettings, get_bulk()));
    }
  }
}

class MeshTest : public stk::unit_test_util::MeshFixture {};

TEST_F(MeshTest, forExodusFile)
{
  std::string filename = stk::unit_test_util::get_option("-i", "generated:1x1x100");
  setup_mesh(filename, stk::mesh::BulkData::AUTO_AURA);

  stk::balance::GraphCreationSettings graphSettings;
  bool foundMechanisms = stk::balance::internal::detectAndFixMechanisms(graphSettings, get_bulk());

  if(foundMechanisms && get_bulk().parallel_rank()==0)
    std::cerr << "Found and fixed mechanisms\n";

  stk::io::write_mesh("junk.g", get_bulk());
}


} // namespace
