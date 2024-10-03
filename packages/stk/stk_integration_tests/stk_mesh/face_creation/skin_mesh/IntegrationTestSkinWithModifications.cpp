/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include "mpi.h"                        // for MPI_COMM_WORLD
#include <stddef.h>                     // for size_t, nullptr
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>
#include <stk_mesh/base/FaceCreator.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/SkinMeshUtil.hpp>
#include <string>                       // for string
#include <Ionit_Initializer.h>     // for Initializer
#include "stk_io/DatabasePurpose.hpp"
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>

namespace
{

class SkinWithModification : public stk::unit_test_util::MeshFixture
{
protected:
  SkinWithModification() : boundaryPart(nullptr)
  {
  }
  void setup_for_skinning()
  {
    boundaryPart = &get_meta().declare_part("boundary", get_meta().side_rank());

    thingsToSkin = get_meta().universal_part();
  }
  ~SkinWithModification()
  {
  }
  void test_skinning(const SideTestUtil::TestCase& exteriorTestCase,
                     const SideTestUtil::TestCase& interiorTestCase)
  {
    test_exposed_boundary(exteriorTestCase);
    test_interior_block_boundary(interiorTestCase);
  }
  void destroy_element(stk::mesh::Entity elem)
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank() == 0)
      get_bulk().destroy_entity(elem);
    get_bulk().modification_end();
  }

  void delete_shells_13_and_14()
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank() == 0)
    {
      get_bulk().destroy_entity(shell13);
      get_bulk().destroy_entity(shell14);
    }
    get_bulk().modification_end();
  }

  void create_shells_13_and_14()
  {
    create_shell_13();
    create_shell_14();
  }

  void create_shell_13()
  {
    shell13 = create_shell_with_id(shellId1);
  }

  void create_shell_14()
  {
    shell14 = create_shell_with_id(shellId2);
  }
  stk::mesh::Entity create_shell_with_id(const stk::mesh::EntityId shellId)
  {
    stk::mesh::Entity shell;
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank() == 0)
      shell = stk::mesh::declare_element(get_bulk(), get_meta().get_topology_root_part(stk::topology::SHELL_QUAD_4), shellId, {5, 6, 7, 8});
    get_bulk().modification_end();
    return shell;
  }

private:
  void test_exposed_boundary(const SideTestUtil::TestCase& testCase)
  {
    make_exterior_sides();
    SideTestUtil::expect_exposed_sides_connected_as_specified_in_test_case(get_bulk(), testCase, get_things_to_skin(), *boundaryPart);
    destroy_sides_in_boundary_part();
  }
  void make_exterior_sides()
  {
    stk::mesh::create_exposed_block_boundary_sides(get_bulk(), get_things_to_skin(), {boundaryPart});
  }

  void destroy_sides_in_boundary_part()
  {
    stk::mesh::EntityVector boundarySides;
    stk::mesh::get_selected_entities(*boundaryPart, get_bulk().buckets(get_meta().side_rank()), boundarySides);
    destroy_boundary_part_entities(boundarySides);
    SideTestUtil::expect_global_num_sides_in_part(get_bulk(), 0, get_meta().universal_part());
  }
  void destroy_entities(const stk::mesh::EntityVector& entities)
  {
    for(stk::mesh::Entity side : entities)
      destroy_side(side);
  }
  void destroy_side(stk::mesh::Entity side)
  {
    stk::mesh::EntityVector elements(get_bulk().begin_elements(side), get_bulk().end_elements(side));
    std::vector<stk::mesh::ConnectivityOrdinal> elementOrdinals(get_bulk().begin_element_ordinals(side), get_bulk().end_element_ordinals(side));
    for(size_t i = 0; i < elements.size(); i++)
      get_bulk().destroy_relation(elements[i], side, elementOrdinals[i]);
    get_bulk().destroy_entity(side);
  }
  void test_interior_block_boundary(const SideTestUtil::TestCase& testCase)
  {
    make_interior_block_boundary_sides();
    SideTestUtil::expect_interior_sides_connected_as_specified_in_test_case(get_bulk(), testCase, get_things_to_skin(), *boundaryPart);
    destroy_sides_in_boundary_part();
  }
  void make_interior_block_boundary_sides()
  {
    stk::mesh::create_interior_block_boundary_sides(get_bulk(), get_things_to_skin(), {boundaryPart});
  }

  const stk::mesh::Selector &get_things_to_skin()
  {
    return thingsToSkin;
  }

  void destroy_boundary_part_entities(stk::mesh::EntityVector &boundarySides)
  {
    get_bulk().modification_begin();
    destroy_entities(boundarySides);
    get_bulk().modification_end();
  }

protected:
  const stk::mesh::EntityId shellId1 = 13;
  const stk::mesh::EntityId shellId2 = 14;

  stk::mesh::Entity shell13;
  stk::mesh::Entity shell14;
private:
  stk::mesh::Part *boundaryPart;
  stk::mesh::Selector thingsToSkin;
};

class SkinFileWithModification : public SkinWithModification
{
protected:
  void setup_mesh_from_initial_configuration(const std::string &initialConfiguration,
                                             bool emptyMeshAlreadyCreated = false)
  {
    if (!emptyMeshAlreadyCreated) {
      setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    }
    SideTestUtil::read_and_decompose_mesh(initialConfiguration, get_bulk());
    setup_for_skinning();
  }

  void test_no_mods(const SideTestUtil::TestCase &exterior, const SideTestUtil::TestCase &interior)
  {
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
      setup_mesh_from_initial_configuration (get_filename());
      test_skinning(exterior, interior);
    }
  }

  void test_adding_one_shell(const SideTestUtil::TestCase &exterior, const SideTestUtil::TestCase &interior)
  {
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
      setup_mesh_from_initial_configuration(get_filename());
      create_shell_13();
      test_skinning(exterior, interior);
    }
  }

  void test_adding_two_shell(const SideTestUtil::TestCase &exterior, const SideTestUtil::TestCase &interior)
  {
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
      setup_mesh_from_initial_configuration(get_filename());
      create_shells_13_and_14();
      test_skinning(exterior, interior);
    }
  }

  void test_adding_two_shells_then_delete_one(const SideTestUtil::TestCase &exterior, const SideTestUtil::TestCase &interior)
  {
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
      setup_mesh_from_initial_configuration(get_filename());
      create_shells_13_and_14();
      destroy_element(shell14);
      test_skinning(exterior, interior);
    }
  }

  void test_adding_two_shells_then_delete_both(const SideTestUtil::TestCase &exterior, const SideTestUtil::TestCase &interior)
  {
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
      setup_mesh_from_initial_configuration(get_filename());
      create_shells_13_and_14();
      destroy_element(shell14);
      destroy_element(shell13);
      test_skinning(exterior, interior);
    }
  }

  void test_adding_two_shells_then_delete_both_in_same_mod_cycle(const SideTestUtil::TestCase &exterior, const SideTestUtil::TestCase &interior)
  {
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
      setup_mesh_from_initial_configuration(get_filename());
      create_shells_13_and_14();
      delete_shells_13_and_14();
      test_skinning(exterior, interior);
    }
  }

  void test_adding_partial_coincident_hex()
  {
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
      stk::mesh::EntityId elemId = 30;
      create_partial_coincident_hexes(elemId);
      test_skinning_partial_coincident(elemId);
    }
  }

  stk::mesh::EntityId create_partial_coincident_hexes(stk::mesh::EntityId elemId)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::mesh::Part& block2 = create_part_with_id(get_meta(), 2, stk::topology::HEX_8);
    const bool emptyMeshAlreadyCreated = true;
    setup_mesh_from_initial_configuration (get_filename(), emptyMeshAlreadyCreated);
    add_partially_coincident_hex(elemId, block2);
    return elemId;
  }

  void add_partially_coincident_hex(stk::mesh::EntityId elemId, stk::mesh::Part& block2)
  {
    get_bulk().modification_begin();
    if(get_bulk().parallel_rank() == 0)
      stk::mesh::declare_element(get_bulk(), block2, elemId, stk::mesh::EntityIdVector {31, 32, 33, 34, 5, 6, 7, 8});
    get_bulk().modification_end();
  }

  void test_skinning_partial_coincident(stk::mesh::EntityId elemId)
  {
    const SideTestUtil::TestCase exteriorCase = {"AA.e", 2, 15, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4},
                                                                 {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5},
                                                                 {elemId, 0}, {elemId, 1}, {elemId, 2}, {elemId, 3}, {elemId, 4}}};
    const SideTestUtil::TestCase interiorCase = {"AA.e", 2, 1, {{1, 5}, {elemId, 5}, {2, 4}}};
    test_skinning(exteriorCase, interiorCase);
  }


  void create_2_quads_with_partial_coincident_quad_mesh(stk::mesh::Part& block1, stk::mesh::Part& block2)
  {
    std::string meshDesc;
    if (get_parallel_size() == 1) {
      meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                 "0,2,QUAD_4_2D,4,3,5,6\n"
                 "0,3,QUAD_4_2D,21,22,3,4";
    }
    else {
      meshDesc = "0,1,QUAD_4_2D,1,2,3,4\n"
                 "1,2,QUAD_4_2D,4,3,5,6\n"
                 "0,3,QUAD_4_2D,21,22,3,4";
    }
    stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

    get_bulk().modification_begin();
    put_entity_into_part(get_bulk(), 1, block1);
    put_entity_into_part(get_bulk(), 2, block2);
    put_entity_into_part(get_bulk(), 3, block1);
    get_bulk().modification_end();
  }

  void put_entity_into_part(stk::mesh::BulkData &bulkData, stk::mesh::EntityId id, stk::mesh::Part& part)
  {
    stk::mesh::Entity entity = bulkData.get_entity(stk::topology::ELEM_RANK, id);
    if(bulkData.is_valid(entity) && bulkData.bucket(entity).owned())
    {
      bulkData.change_entity_parts(entity, stk::mesh::ConstPartVector{&part});
    }
  }
  stk::mesh::Part& create_part_with_id(stk::mesh::MetaData &metaData, int id, stk::topology topology)
  {
    stk::mesh::Part& part = metaData.declare_part_with_topology("block_"+std::to_string(id), topology);
    metaData.set_part_id(part, id);
    return part;
  }

  virtual const std::string get_filename() = 0;
};

class SkinAAWithModification : public SkinFileWithModification
{
protected:
  virtual const std::string get_filename() { return "AA.e"; }
  const SideTestUtil::TestCase AAExterior =   {"AA.e",   2, 10, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}};
  const SideTestUtil::TestCase AeAExterior =  {"AeA.e",   2, 10, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}};
  const SideTestUtil::TestCase AefAExterior = {"AefA.e",   2, 10, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}};
  const SideTestUtil::TestCase AAInterior =   {"AA.e",   2,  0, {}};
  const SideTestUtil::TestCase AeAInterior =  {"AeA.e",  3,  2, {{1, 5}, {shellId1, 0}, {shellId1, 1}, {2, 4}}};
  const SideTestUtil::TestCase AefAInterior = {"AefA.e", 3,  2, {{1, 5}, {shellId1, 0}, {shellId1, 1}, {shellId2, 0}, {shellId2, 1}, {2, 4}}};
};
TEST_F(SkinAAWithModification, TestSkinningWithNoMods)
{
  test_no_mods(AAExterior, AAInterior);
}
TEST_F(SkinAAWithModification, TestAddingOneShell)
{
  test_adding_one_shell(AeAExterior, AeAInterior);
}
TEST_F(SkinAAWithModification, TestAddingTwoShells)
{
  test_adding_two_shell(AefAExterior, AefAInterior);
}
TEST_F(SkinAAWithModification, TestAddTwoShellThenDeleteOne)
{
  test_adding_two_shells_then_delete_one(AeAExterior, AeAInterior);
}
TEST_F(SkinAAWithModification, TestAddTwoShellThenDeleteBoth)
{
  test_adding_two_shells_then_delete_both(AAExterior, AAInterior);
}
TEST_F(SkinAAWithModification, TestAddTwoShellThenDeleteBothInSameModCycle)
{
  test_adding_two_shells_then_delete_both_in_same_mod_cycle(AAExterior, AAInterior);
}
TEST_F(SkinAAWithModification, TestPartialCoincident)
{
  test_adding_partial_coincident_hex();
}
TEST_F(SkinAAWithModification, TestPartialCoincident2d)
{
  if(stk::parallel_machine_size(get_comm()) <= 2)
  {
    reset_mesh();
    set_spatial_dimension(2);
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    stk::mesh::Part& block1 = create_part_with_id(get_meta(), 1, stk::topology::QUAD_4_2D);
    stk::mesh::Part& block2 = create_part_with_id(get_meta(), 2, stk::topology::QUAD_4_2D);

    setup_for_skinning();

    create_2_quads_with_partial_coincident_quad_mesh(block1, block2);

    const SideTestUtil::TestCase exteriorCase = {"AA.e", 2, 9, {{1, 0}, {1, 1}, {1, 3},
                                                                {2, 1}, {2, 2}, {2, 3},
                                                                {3, 0}, {3, 1}, {3, 3}}};
    const SideTestUtil::TestCase interiorCase = {"AA.e", 2, 1, {{1, 2}, {3, 2}, {2, 0}}};
    test_skinning(exteriorCase, interiorCase);
  }
}

class SkinAWithModification : public SkinFileWithModification
{
protected:
  virtual const std::string get_filename() { return "A.e"; }
  const SideTestUtil::TestCase AExterior =   {"A.e",   1,  6, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}}};
  const SideTestUtil::TestCase AeExterior =  {"Ae.e",  2,  6, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {shellId1, 0}}};
  const SideTestUtil::TestCase AefExterior = {"Aef.e", 3,  6, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {shellId1, 0}, {shellId2, 0}}};
  const SideTestUtil::TestCase AInterior =   {"A.e",   2,  0, {}};
  const SideTestUtil::TestCase AeInterior =  {"Ae.e",  2,  1, {{1, 5}, {shellId1, 1}}};
  const SideTestUtil::TestCase AefInterior = {"Aef.e", 3,  1, {{1, 5}, {shellId1, 1}, {shellId2, 1}}};
};
TEST_F(SkinAWithModification, TestSkinningWithNoMods)
{
  test_no_mods(AExterior, AInterior);
}
TEST_F(SkinAWithModification, TestAddingOneShell)
{
  test_adding_one_shell(AeExterior, AeInterior);
}
TEST_F(SkinAWithModification, TestAddingTwoShells)
{
  test_adding_two_shell(AefExterior, AefInterior);
}
TEST_F(SkinAWithModification, TestAddTwoShellThenDeleteOne)
{
  test_adding_two_shells_then_delete_one(AeExterior, AeInterior);
}
TEST_F(SkinAWithModification, TestAddTwoShellThenDeleteBoth)
{
  test_adding_two_shells_then_delete_both(AExterior, AInterior);
}

}
