#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include "../../stk_mesh/stk_mesh/base/FEMHelpers.hpp"
#include "../../stk_mesh/stk_mesh/base/GetEntities.hpp"
#include <stk_mesh/base/DestroyElements.hpp>

namespace
{

void expect_valid(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities)
{
    for(stk::mesh::Entity entity : entities) {
        EXPECT_TRUE(bulk.is_valid(entity));
    }
}

void expect_invalid(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities)
{
    for(stk::mesh::Entity entity : entities) {
        EXPECT_FALSE(bulk.is_valid(entity));
    }
}

void expect_not_shared(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities)
{
    for(stk::mesh::Entity entity : entities) {
        EXPECT_TRUE(bulk.is_valid(entity));
        EXPECT_FALSE(bulk.bucket(entity).shared());
    }
}

void expect_shared(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &entities)
{
    for(stk::mesh::Entity entity : entities) {
        EXPECT_TRUE(bulk.is_valid(entity));
        EXPECT_TRUE(bulk.bucket(entity).shared());
    }
}

stk::mesh::EntityVector get_faces_for_entity(const stk::mesh::BulkData &bulk, const stk::mesh::Entity entity)
{
    unsigned numConnected = bulk.num_connectivity(entity, stk::topology::FACE_RANK);
    const stk::mesh::Entity* connectedEntities = bulk.begin(entity, stk::topology::FACE_RANK);
    stk::mesh::EntityVector entityFaces(connectedEntities, connectedEntities+numConnected);
    return entityFaces;
}

class HexMesh : public stk::unit_test_util::MeshFixture
{
protected:
    HexMesh()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        std::string meshDesc =
            "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
             0,2,HEX_8,2,9,10,3,6,11,12,7";
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
    }
};

TEST_F(HexMesh, DeleteOneElement)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
        stk::mesh::EntityVector orphanedNodes{
                    get_bulk().get_entity(stk::topology::NODE_RANK, 9),
                    get_bulk().get_entity(stk::topology::NODE_RANK, 10),
                    get_bulk().get_entity(stk::topology::NODE_RANK, 11),
                    get_bulk().get_entity(stk::topology::NODE_RANK, 12)
        };

        stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

        expect_valid(get_bulk(), orphanedNodes);
        expect_valid(get_bulk(), elementToDestroy);

        stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

        expect_invalid(get_bulk(), orphanedNodes);
        expect_invalid(get_bulk(), elementToDestroy);
    }
}

class TetMesh : public stk::unit_test_util::MeshFixture
{
protected:
    TetMesh()
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        std::string meshDesc =
               "0,1,TET_4,1,2,3,4\n\
                0,2,TET_4,2,5,3,4";

        if(get_bulk().parallel_size() == 2)
        {
            meshDesc =  "0,1,TET_4,1,2,3,4\n\
                         1,2,TET_4,2,5,3,4";
        }
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
    }
};

TEST_F(TetMesh, DeleteOneElement)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
        stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);

        stk::mesh::EntityVector orphanedNodes{
            get_bulk().get_entity(stk::topology::NODE_RANK, 5)
        };

        stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

        stk::mesh::EntityVector facesOfDestroyedElement = get_faces_for_entity(get_bulk(), elementToDestroy[0]);
        EXPECT_EQ(facesOfDestroyedElement.size(), 4u);
        for(stk::mesh::Entity face : facesOfDestroyedElement)
            EXPECT_TRUE(get_bulk().is_valid(face));

        expect_valid(get_bulk(), orphanedNodes);
        expect_valid(get_bulk(), elementToDestroy);

        stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

        expect_invalid(get_bulk(), orphanedNodes);
        expect_invalid(get_bulk(), elementToDestroy);

        unsigned numValid = 0;
        for(stk::mesh::Entity face : facesOfDestroyedElement)
            if(get_bulk().is_valid(face))
                numValid++;

        EXPECT_EQ(1u, numValid);
    }
}

TEST_F(TetMesh, DeleteElementOnProcBoundaryWithOwnedFace)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), {}, false);

        stk::mesh::EntityVector orphanedNodes{
            get_bulk().get_entity(stk::topology::NODE_RANK, 1)
        };

        stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1)};

        stk::mesh::EntityVector facesOfDestroyedElement;

        if(get_bulk().parallel_rank() == 0)
        {
            facesOfDestroyedElement = get_faces_for_entity(get_bulk(), elementToDestroy[0]);
            EXPECT_EQ(facesOfDestroyedElement.size(), 4u);

            for(stk::mesh::Entity face : facesOfDestroyedElement)
                EXPECT_TRUE(get_bulk().is_valid(face));

            expect_valid(get_bulk(), orphanedNodes);
            expect_valid(get_bulk(), elementToDestroy);
        }
        else if(get_bulk().parallel_rank() == 1)
        {
            stk::mesh::Entity element2 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2);
            EXPECT_TRUE(get_bulk().is_valid(element2));
        }
        stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

        expect_invalid(get_bulk(), orphanedNodes);
        expect_invalid(get_bulk(), elementToDestroy);

        if(get_bulk().parallel_rank() == 0)
        {
            for(stk::mesh::Entity face : facesOfDestroyedElement)
                EXPECT_FALSE(get_bulk().is_valid(face));
        }
        else if(get_bulk().parallel_rank() == 1)
        {
            stk::mesh::EntityVector faces = get_faces_for_entity(get_bulk(), get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2));
            expect_not_shared(get_bulk(), faces);
        }
    }
}

TEST_F(TetMesh, DeleteGhostedElement)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        stk::mesh::EntityVector orphanedNodes{
            get_bulk().get_entity(stk::topology::NODE_RANK, 1)
        };

        stk::mesh::EntityVector sharedNodes{
            get_bulk().get_entity(stk::topology::NODE_RANK, 2),
            get_bulk().get_entity(stk::topology::NODE_RANK, 3),
            get_bulk().get_entity(stk::topology::NODE_RANK, 4)
        };

        get_bulk().modification_begin();
        stk::mesh::Ghosting &ghosting = get_bulk().create_ghosting("ghost");
        stk::mesh::EntityProcVec ghostedElement;
        if(get_bulk().parallel_rank() == 0)
        {
            stk::mesh::Entity element1 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1);
            ghostedElement = {stk::mesh::EntityProc(element1, 1)};
        }
        get_bulk().change_ghosting(ghosting, ghostedElement);
        get_bulk().modification_end();

        stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 1)};
        expect_valid(get_bulk(), elementToDestroy);

        if(get_bulk().parallel_rank() == 0)
        {
            expect_valid(get_bulk(), orphanedNodes);
        }
        else if(get_bulk().parallel_rank() == 1)
        {
            stk::mesh::Entity element2 = get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2);
            EXPECT_TRUE(get_bulk().is_valid(element2));
        }
        stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

        expect_invalid(get_bulk(), orphanedNodes);
        expect_invalid(get_bulk(), elementToDestroy);

        if(get_bulk().parallel_rank() == 0)
        {
            expect_invalid(get_bulk(), sharedNodes);
        }
        else if(get_bulk().parallel_rank() == 1)
        {
            expect_not_shared(get_bulk(), sharedNodes);
        }
    }
}

class BeamMesh : public stk::unit_test_util::MeshFixture
{
protected:
    BeamMesh()
    {
        delete_meta();
        allocate_meta(2);
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        std::string meshDesc =
            "0,1,BEAM_2,1,2\n\
             0,2,BEAM_2,2,3";
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
    }
};

TEST_F(BeamMesh, DeleteOneElement)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    {
        stk::mesh::EntityVector orphanedNodes{
            get_bulk().get_entity(stk::topology::NODE_RANK, 3)
        };

        stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

        expect_valid(get_bulk(), orphanedNodes);
        expect_valid(get_bulk(), elementToDestroy);

        stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

        expect_invalid(get_bulk(), orphanedNodes);
        expect_invalid(get_bulk(), elementToDestroy);
    }
}

class QuadMesh : public stk::unit_test_util::MeshFixture
{
protected:
    QuadMesh()
    {
    }

    void setup_my_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        delete_meta();
        allocate_meta(2);
        setup_empty_mesh(auraOption);
        std::string meshDesc =
            "0,1,QUAD_4_2D,1,2,3,4\n\
             1,2,QUAD_4_2D,2,5,6,3";
        stk::unit_test_util::fill_mesh_using_text_mesh(meshDesc, get_bulk());
    }

    void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_my_mesh(auraOption);

        stk::mesh::EntityVector orphanedNodes{
            get_bulk().get_entity(stk::topology::NODE_RANK, 5),
                    get_bulk().get_entity(stk::topology::NODE_RANK, 6)
        };

        stk::mesh::EntityVector sharedNodes{
            get_bulk().get_entity(stk::topology::NODE_RANK, 2),
                    get_bulk().get_entity(stk::topology::NODE_RANK, 3)
        };

        stk::mesh::EntityVector elementToDestroy{get_bulk().get_entity(stk::topology::ELEMENT_RANK, 2)};

        expect_shared(get_bulk(), sharedNodes);

        if(get_bulk().parallel_rank() == 1)
        {
            expect_valid(get_bulk(), orphanedNodes);
            expect_valid(get_bulk(), elementToDestroy);
        }

        stk::mesh::destroy_elements(get_bulk(), elementToDestroy);

        expect_invalid(get_bulk(), orphanedNodes);
        expect_invalid(get_bulk(), elementToDestroy);

        if(get_bulk().parallel_rank() == 1)
        {
            expect_invalid(get_bulk(), sharedNodes);
        }
        else if(get_bulk().parallel_rank() == 0)
        {
            expect_not_shared(get_bulk(), sharedNodes);
        }
    }
};

TEST_F(QuadMesh, DeleteProcBoundaryElementWithoutAura)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
        run_test(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(QuadMesh, DeleteProcBoundaryElementWithAura)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
        run_test(stk::mesh::BulkData::AUTO_AURA);
}

}
