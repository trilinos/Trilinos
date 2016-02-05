/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include "mpi.h"                        // for MPI_COMM_WORLD
#include <stddef.h>                     // for size_t, nullptr
#include <string>                       // for string
#include <Ioss_IOFactory.h>             // for IOFactory
#include <Ioss_Region.h>                // for Region
#include <init/Ionit_Initializer.h>     // for Initializer
#include "stk_io/DatabasePurpose.hpp"
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include "../FaceCreationTestUtils.hpp"

namespace
{

class SkinWithModification : public stk::unit_test_util::MeshFixture
{
protected:
    SkinWithModification()
    {
        if(stk::parallel_machine_size(get_comm()) <= 2)
        {
            setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
            SideTestUtil::read_and_decompose_mesh("AA.e", get_bulk());
            boundaryPart = &get_meta().declare_part("boundary", get_meta().side_rank());
        }
    }
    void test_skinning(const SideTestUtil::TestCase& exteriorTestCase,
                       const SideTestUtil::TestCase& interiorTestCase)
    {
        test_exposed_boundary(exteriorTestCase);
        test_interior_block_boundary(interiorTestCase);
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
    void destroy_element(stk::mesh::Entity elem)
    {
        get_bulk().modification_begin();
        if(get_bulk().parallel_rank() == 0)
            get_bulk().destroy_entity(elem);
        get_bulk().modification_end();
    }
private:
    void test_exposed_boundary(const SideTestUtil::TestCase& testCase)
    {
        stk::mesh::create_exposed_boundary_sides(get_bulk(), get_things_to_skin(), *boundaryPart);
        SideTestUtil::expect_exposed_sides_connected_as_specified_in_test_case(get_bulk(), testCase, get_things_to_skin(), *boundaryPart);
        destroy_sides_in_boundary_part();
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
        SideTestUtil::create_interior_block_boundary_sides(get_bulk(), get_things_to_skin(), *boundaryPart);
        SideTestUtil::expect_interior_sides_connected_as_specified_in_test_case(get_bulk(), testCase, get_things_to_skin(), *boundaryPart);
        destroy_sides_in_boundary_part();
    }
    stk::mesh::Selector get_things_to_skin()
    {
        return get_meta().universal_part();
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
    const SideTestUtil::TestCase AAExterior =   {"AA.e",   2, 10, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}};
    const SideTestUtil::TestCase AeAExterior =  {"AA.e",   2, 10, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}};
    const SideTestUtil::TestCase AefAExterior = {"AA.e",   2, 10, {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}};
    const SideTestUtil::TestCase AAInterior =   {"AA.e",   2,  0, {}};
    const SideTestUtil::TestCase AeAInterior =  {"AeA.e",  3,  2, {{1, 5}, {shellId1, 0}, {shellId1, 1}, {2, 4}}};
    const SideTestUtil::TestCase AefAInterior = {"AefA.e", 3,  2, {{1, 5}, {shellId1, 0}, {shellId1, 1}, {shellId2, 0}, {shellId2, 1}, {2, 4}}};
private:
    stk::mesh::Part *boundaryPart;
};

TEST_F(SkinWithModification, BuildUpThenTearDownWithAura)
{
    if(stk::parallel_machine_size(get_comm()) <= 2)
    {
        test_skinning(AAExterior, AAInterior);

        stk::mesh::Entity shell13 = create_shell_with_id(shellId1);
        test_skinning(AeAExterior, AeAInterior);

        stk::mesh::Entity shell14 = create_shell_with_id(shellId2);
        test_skinning(AefAExterior, AefAInterior);

        destroy_element(shell14);
        test_skinning(AeAExterior, AeAInterior);

        destroy_element(shell13);
        test_skinning(AAExterior, AAInterior);
    }
}
}
