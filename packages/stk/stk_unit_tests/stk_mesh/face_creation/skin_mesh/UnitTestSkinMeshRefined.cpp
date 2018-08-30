#include <gtest/gtest.h>

#include "UnitTestSkinMeshUseCaseUtils.hpp"
#include <stk_mesh/base/SkinBoundary.hpp>

namespace {

using namespace stk::mesh::impl;
using namespace stk::mesh;


TEST(ElementGraph, RefinedQuad)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;
    if (stk::parallel_machine_size(comm) <= 2)
    {
        const int spatialDim = 2;
        stk::mesh::MetaData meta(spatialDim);

        stk::mesh::Part &quad_part = meta.declare_part_with_topology("Quads", stk::topology::QUADRILATERAL_4_2D);
        stk::mesh::Part &skin = meta.declare_part_with_topology("Edges", stk::topology::LINE_2);
        stk::io::put_io_part_attribute(skin);
        stk::mesh::PartVector skin_parts = {&skin};
        stk::mesh::Part &active = meta.declare_part("active");

        stk::mesh::Field<double,stk::mesh::Cartesian> & node_coord = meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian>>(stk::topology::NODE_RANK, "coordinates");
        stk::mesh::put_field_on_mesh(node_coord, meta.universal_part(), 3, nullptr);
        stk::io::put_io_part_attribute(quad_part);
        meta.commit();

        stk::mesh::BulkData mesh(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);

        mesh.modification_begin();

        std::vector<std::vector<stk::mesh::EntityId>> elem_this_proc = {
                {1, 2, 3, 4, 5, 6},
                {7}
        };

        std::vector<stk::mesh::EntityId> ids = {1, 2, 3, 4, 5, 6, 7};
        std::vector<std::vector<stk::mesh::EntityId> > connectivity = {
                {1, 2, 3, 4},
                {1, 2, 6, 5},
                {2, 3, 7, 6},
                {8, 7, 3, 4},
                {1, 5, 8, 4},
                {5, 6, 7, 8},
                {10, 1, 4, 9}
        };

        bool running_2_procs = mesh.parallel_size() > 1;

        int my_proc_id = mesh.parallel_rank();
        std::vector<stk::mesh::EntityId> elems_this_proc;
        if(running_2_procs)
        {
            elems_this_proc = elem_this_proc[my_proc_id];
        }
        else
        {
            elems_this_proc = ids;
        }

        for(size_t i = 0; i < elems_this_proc.size(); ++i)
        {
            int index = static_cast<int>(elems_this_proc[i])-1;
            stk::mesh::Entity element = stk::mesh::declare_element(mesh, quad_part, elems_this_proc[i], connectivity[index]);
            if(ids[index]!=1)
            {
                mesh.change_entity_parts(element, stk::mesh::ConstPartVector{&active}, stk::mesh::ConstPartVector{});
            }
        }

        std::vector<stk::mesh::EntityId> nodes = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

        std::vector<std::vector<stk::mesh::EntityId>> nodes_this_proc = {
                {1, 2, 3, 4, 5, 6, 7, 8},
                {1, 4, 9, 10}
        };

        std::vector<std::vector<double> > node_coords = {
                {0, 0}, {3, 0}, {3, 3}, {0, 3}, {1, 1}, {2, 1}, {2, 2}, {1, 2}, {-3, 3}, {-3, 0}
        };

        std::vector<stk::mesh::EntityId> nodes_to_create;
        if(running_2_procs)
        {
            nodes_to_create = nodes_this_proc[my_proc_id];
        }
        else
        {
            nodes_to_create = nodes;
        }

        for(unsigned i = 0; i < nodes_to_create.size(); ++i)
        {
            stk::mesh::Entity const node = mesh.get_entity(stk::topology::NODE_RANK, nodes_to_create[i]);
            if(running_2_procs && (nodes_to_create[i]==1 || nodes_to_create[i]==4))
            {
                mesh.add_node_sharing(node,1-my_proc_id);
            }

            double * const coord = stk::mesh::field_data(node_coord, node);
            int index = static_cast<int>(nodes_to_create[i])-1;

            coord[0] = node_coords[index][0];
            coord[1] = node_coords[index][1];
        }

        mesh.modification_end();

        stk::mesh::Selector active_sel = active;
        stk::mesh::Selector air = !active;

        stk::mesh::create_exposed_block_boundary_sides(mesh, active_sel, {}, air);

        std::vector<size_t> mesh_counts;
        stk::mesh::comm_mesh_counts(mesh, mesh_counts);
        EXPECT_EQ(10u, mesh_counts[stk::topology::NODE_RANK]);
        EXPECT_EQ(7u, mesh_counts[stk::topology::ELEM_RANK]);
        EXPECT_EQ(6u, mesh_counts[meta.side_rank()]);

        stk::io::write_mesh("refined.g", mesh);

        std::vector<std::pair<stk::mesh::EntityId, int>> id_and_num_faces = {
                {1, 3},
                {2, 1},
                {3, 1},
                {4, 1},
                {5, 0},
                {6, 0},
                {7, 3}
        };

        for(size_t i = 0; i < id_and_num_faces.size(); ++i)
        {
            stk::mesh::EntityId id = id_and_num_faces[i].first;
            int gold_num_faces = id_and_num_faces[i].second;
            stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, id);
            if(mesh.is_valid(elem))
            {
                int num_faces = mesh.num_sides(elem);
                EXPECT_EQ(gold_num_faces, num_faces)<< "element " << id << " has topology " << mesh.bucket(elem).topology() << " with num faces " << num_faces << " not same as gold value " << gold_num_faces << std::endl;
            }
        }
    }
}

} // namespace
