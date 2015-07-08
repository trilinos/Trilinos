#include <gtest/gtest.h>

#include <vector>
#include <algorithm>
#include <stdlib.h>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/ElemElemGraph.hpp>
#include <stk_mesh/base/ElemElemGraphImpl.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/getOption.h>

#include "UnitTestElementDeathUtils.hpp"

namespace
{



void test_active_part_membership(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& skin_faces_of_elem2, stk::mesh::Part& active)
{
    for(stk::topology::rank_t rank = stk::topology::BEGIN_RANK; rank < bulkData.mesh_meta_data().entity_rank_count(); rank++)
    {
        const stk::mesh::BucketVector &buckets = bulkData.get_buckets(rank, bulkData.mesh_meta_data().locally_owned_part());
        for(const stk::mesh::Bucket *bucket : buckets)
        {
            for(stk::mesh::Entity entity : *bucket)
            {
                bool skin_face_of_elem2 = false;
                if(rank == stk::topology::FACE_RANK)
                {
                    stk::mesh::EntityVector::iterator iter = std::find(skin_faces_of_elem2.begin(), skin_faces_of_elem2.end(), entity);
                    if(iter != skin_faces_of_elem2.end())
                    {
                        skin_face_of_elem2 = true;
                    }
                }

                if((rank == stk::topology::ELEM_RANK && bulkData.identifier(entity) == 2) || skin_face_of_elem2)
                {
                    EXPECT_FALSE(bucket->member(active))<< " is entity inactive: " << bulkData.entity_key(entity);
                }
                else
                {
                    EXPECT_TRUE(bucket->member(active)) << " is entity active" << bulkData.entity_key(entity);
                }
            }
        }
    }
}

void test_face_membership_for_death(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& internal_faces_of_elem2, stk::mesh::PartVector& boundary_mesh_parts)
{
    const stk::mesh::BucketVector &buckets = bulkData.get_buckets(stk::topology::FACE_RANK, bulkData.mesh_meta_data().locally_owned_part());
    for(const stk::mesh::Bucket *bucket : buckets)
    {
        for(stk::mesh::Entity entity : *bucket)
        {
            bool internal_face_of_elem2 = false;
            stk::mesh::EntityVector::iterator iter = std::find(internal_faces_of_elem2.begin(), internal_faces_of_elem2.end(), entity);
            if(iter!=internal_faces_of_elem2.end())
            {
                internal_face_of_elem2 = true;
            }

            if(internal_face_of_elem2)
            {
                EXPECT_TRUE(bucket->member_all(boundary_mesh_parts)) << " is entity in boundary parts: " << bulkData.entity_key(entity);
            }
            else
            {
                EXPECT_FALSE(bucket->member_all(boundary_mesh_parts)) << " is entity not in boundary parts" << bulkData.entity_key(entity);
            }
        }
    }
}

TEST(ElementDeath, keep_faces_after_element_death_after_calling_create_faces)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) <= 2)
    {
        unsigned spatialDim = 3;

        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::mesh::PartVector boundary_mesh_parts {&faces_part};
        stk::io::put_io_part_attribute(faces_part);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::mesh::Part& active = meta.declare_part("active"); // can't specify rank, because it gets checked against size of rank_names

        ASSERT_TRUE(active.primary_entity_rank() == stk::topology::INVALID_RANK);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        stk::mesh::create_faces(bulkData);

        ElementDeathUtils::put_mesh_into_part(bulkData, active);

        stk::mesh::ElemElemGraph graph(bulkData, active);

        size_t num_gold_edges = 6 / bulkData.parallel_size();
        ASSERT_EQ(num_gold_edges, graph.num_edges());

        stk::mesh::EntityVector deactivated_elems;

        stk::mesh::EntityId elem2Id = 2;
        stk::mesh::EntityId elem3Id = 3;

        stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem2Id);
        stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, elem3Id);

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(bulkData, entity_counts);
        ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 21);

        {
            stk::mesh::EntityVector skin_faces_of_elem2;
            stk::mesh::EntityVector internal_faces_of_elem2;
            if(bulkData.is_valid(elem2) && bulkData.bucket(elem2).owned())
            {
                deactivated_elems.push_back(elem2);
                unsigned num_faces = bulkData.num_faces(elem2);
                const stk::mesh::Entity* faces = bulkData.begin_faces(elem2);
                for(unsigned j=0;j<num_faces;++j)
                {
                    if(bulkData.num_elements(faces[j])==1)
                    {
                        skin_faces_of_elem2.push_back(faces[j]);
                    }
                    else
                    {
                        internal_faces_of_elem2.push_back(faces[j]);
                    }
                }
            }

            ElementDeathUtils::deactivate_elements(deactivated_elems, bulkData,  active);

            test_active_part_membership(bulkData, skin_faces_of_elem2, active);

            stk::mesh::perform_element_death(bulkData, graph, deactivated_elems, active, boundary_mesh_parts);

            stk::mesh::Entity face_between_elem2_and_elem3 = ElementDeathUtils::get_face_between_element_ids(graph, bulkData, elem2Id, elem3Id);

            ASSERT_TRUE(bulkData.is_valid(face_between_elem2_and_elem3));

            test_face_membership_for_death(bulkData, internal_faces_of_elem2, boundary_mesh_parts);
        }

        stk::mesh::comm_mesh_counts(bulkData, entity_counts);

        ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 21);

        // now kill 3

        {
            deactivated_elems.clear();

            stk::mesh::EntityVector skin_faces_of_elem3;
            stk::mesh::EntityVector internal_faces_of_elem3;
            if(bulkData.is_valid(elem3) && bulkData.bucket(elem3).owned())
            {
                deactivated_elems.push_back(elem3);
                unsigned num_faces = bulkData.num_faces(elem3);
                const stk::mesh::Entity* faces = bulkData.begin_faces(elem3);
                for(unsigned j=0;j<num_faces;++j)
                {
                    if(bulkData.num_elements(faces[j])==1)
                    {
                        skin_faces_of_elem3.push_back(faces[j]);
                    }
                    else
                    {
                        internal_faces_of_elem3.push_back(faces[j]);
                    }
                }
            }

            ElementDeathUtils::deactivate_elements(deactivated_elems, bulkData,  active);

            stk::mesh::Entity face_between_elem2_and_elem3 = ElementDeathUtils::get_face_between_element_ids(graph, bulkData, elem2Id, elem3Id);

            stk::mesh::EntityId face_id;

            EXPECT_TRUE(bulkData.is_valid(face_between_elem2_and_elem3));
            face_id = bulkData.identifier(face_between_elem2_and_elem3);
            ASSERT_FALSE(bulkData.bucket(face_between_elem2_and_elem3).member(active));

            stk::mesh::perform_element_death(bulkData, graph, deactivated_elems, active, boundary_mesh_parts);

            stk::mesh::Entity face_23 = bulkData.get_entity(stk::topology::FACE_RANK, face_id);
            EXPECT_TRUE(bulkData.is_valid(face_23));
        }

        stk::mesh::comm_mesh_counts(bulkData, entity_counts);
        ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 21);
    }
}

TEST(ElementDeath, keep_faces_after_element_death_without_calling_create_faces)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    if(stk::parallel_machine_size(comm) <= 2)
    {
        unsigned spatialDim = 3;

        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::Part& faces_part = meta.declare_part_with_topology("surface_5", stk::topology::QUAD_4);
        stk::mesh::PartVector boundary_mesh_parts {&faces_part};
        stk::io::put_io_part_attribute(faces_part);
        stk::mesh::BulkData bulkData(meta, comm);

        stk::mesh::Part& active = meta.declare_part("active"); // can't specify rank, because it gets checked against size of rank_names

        ASSERT_TRUE(active.primary_entity_rank() == stk::topology::INVALID_RANK);

        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", bulkData, comm);

        ElementDeathUtils::put_mesh_into_part(bulkData, active);

        stk::mesh::ElemElemGraph graph(bulkData, active);

        size_t num_gold_edges = 6 / bulkData.parallel_size();
        ASSERT_EQ(num_gold_edges, graph.num_edges());

        stk::mesh::EntityVector deactivated_elems;

        stk::mesh::EntityId elem2Id = 2;
        stk::mesh::EntityId elem3Id = 3;

        stk::mesh::Entity elem2 = bulkData.get_entity(stk::topology::ELEM_RANK, elem2Id);
        stk::mesh::Entity elem3 = bulkData.get_entity(stk::topology::ELEM_RANK, elem3Id);

        std::vector<size_t> entity_counts;
        stk::mesh::comm_mesh_counts(bulkData, entity_counts);
        ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 0);

        {
            stk::mesh::EntityVector skin_faces_of_elem2;
            stk::mesh::EntityVector internal_faces_of_elem2;
            if(bulkData.is_valid(elem2) && bulkData.bucket(elem2).owned())
            {
                deactivated_elems.push_back(elem2);
                unsigned num_faces = bulkData.num_faces(elem2);
                const stk::mesh::Entity* faces = bulkData.begin_faces(elem2);
                for(unsigned j=0;j<num_faces;++j)
                {
                    if(bulkData.num_elements(faces[j])==1)
                    {
                        skin_faces_of_elem2.push_back(faces[j]);
                    }
                    else
                    {
                        internal_faces_of_elem2.push_back(faces[j]);
                    }
                }
            }

            ElementDeathUtils::deactivate_elements(deactivated_elems, bulkData,  active);

            test_active_part_membership(bulkData, skin_faces_of_elem2, active);

            stk::mesh::perform_element_death(bulkData, graph, deactivated_elems, active, boundary_mesh_parts);

            stk::mesh::Entity face_between_elem2_and_elem3 = ElementDeathUtils::get_face_between_element_ids(graph, bulkData, elem2Id, elem3Id);

            ASSERT_TRUE(bulkData.is_valid(face_between_elem2_and_elem3));
        }

        stk::mesh::comm_mesh_counts(bulkData, entity_counts);

        ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 2);

        // now kill 3

        {
            deactivated_elems.clear();

            stk::mesh::EntityVector skin_faces_of_elem3;
            stk::mesh::EntityVector internal_faces_of_elem3;
            if(bulkData.is_valid(elem3) && bulkData.bucket(elem3).owned())
            {
                deactivated_elems.push_back(elem3);
                unsigned num_faces = bulkData.num_faces(elem3);
                const stk::mesh::Entity* faces = bulkData.begin_faces(elem3);
                for(unsigned j=0;j<num_faces;++j)
                {
                    if(bulkData.num_elements(faces[j])==1)
                    {
                        skin_faces_of_elem3.push_back(faces[j]);
                    }
                    else
                    {
                        internal_faces_of_elem3.push_back(faces[j]);
                    }
                }
            }

            ElementDeathUtils::deactivate_elements(deactivated_elems, bulkData,  active);

            stk::mesh::Entity face_between_elem2_and_elem3 = ElementDeathUtils::get_face_between_element_ids(graph, bulkData, elem2Id, elem3Id);

            stk::mesh::perform_element_death(bulkData, graph, deactivated_elems, active, boundary_mesh_parts);

            EXPECT_FALSE(bulkData.is_valid(face_between_elem2_and_elem3));
        }

        stk::mesh::comm_mesh_counts(bulkData, entity_counts);
        ASSERT_TRUE(entity_counts[stk::topology::FACE_RANK] == 2);
    }
}

} // end namespace


