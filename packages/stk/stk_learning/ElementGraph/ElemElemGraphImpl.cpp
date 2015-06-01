#include "ElemElemGraphImpl.hpp"

#include <vector>
#include <algorithm>

#include <stk_topology/topology.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace impl
{

void set_local_ids_and_fill_element_entities_and_topologies(stk::mesh::BulkData& bulkData, stk::mesh::EntityVector& local_id_to_element_entity, std::vector<stk::topology>& element_topologies)
{
    const stk::mesh::BucketVector& elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    size_t local_id = 0;
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        for(size_t j=0; j<bucket.size(); ++j)
        {
            local_id_to_element_entity[local_id] = bucket[j];
            element_topologies[local_id] = bucket.topology();
            bulkData.set_local_id(bucket[j], local_id);
            local_id++;
        }
    }
}

void fill_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph, SidesForElementGraph& via_sides)
{
    const stk::mesh::BucketVector& elemBuckets = bulkData.get_buckets(stk::topology::ELEM_RANK, bulkData.mesh_meta_data().locally_owned_part());
    for(size_t i=0; i<elemBuckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *elemBuckets[i];
        stk::topology topology = bucket.topology();
        int num_sides = topology.num_sides();
        std::vector<ElementSidePair> elem_side_pairs;
        stk::mesh::EntityVector side_nodes;
        stk::mesh::EntityVector connected_elements;
        for(size_t j=0; j<bucket.size(); ++j)
        {
            LocalId local_elem_id = bulkData.local_id(bucket[j]);
            const stk::mesh::Entity* elem_nodes = bucket.begin_nodes(j);
            elem_side_pairs.clear();
            for(int side_index=0; side_index<num_sides; ++side_index)
            {
                unsigned num_side_nodes = topology.side_topology(side_index).num_nodes();
                side_nodes.resize(num_side_nodes);
                topology.side_nodes(elem_nodes, side_index, side_nodes.begin());
                connected_elements.clear();
                stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, num_side_nodes, side_nodes.data(), connected_elements);
                for(size_t elem_index=0; elem_index<connected_elements.size(); ++elem_index)
                {
                    if (connected_elements[elem_index] != bucket[j])
                    {
                        elem_side_pairs.push_back(std::make_pair(bulkData.local_id(connected_elements[elem_index]),side_index));
                    }
                }
            }

            std::sort(elem_side_pairs.begin(), elem_side_pairs.end());
            std::vector<ElementSidePair>::iterator new_end = std::unique(elem_side_pairs.begin(), elem_side_pairs.end());
            elem_side_pairs.resize(new_end - elem_side_pairs.begin());
            for(size_t index=0; index<elem_side_pairs.size(); ++index)
            {
                elem_graph[local_elem_id].push_back(elem_side_pairs[index].first);
                via_sides[local_elem_id].push_back(elem_side_pairs[index].second);
            }
        }
    }
}

stk::mesh::EntityVector get_elements_to_communicate(const stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector elements_to_communicate;
    std::set<stk::mesh::Entity> element_set;
    const stk::mesh::BucketVector& shared_node_buckets = bulkData.get_buckets(stk::topology::NODE_RANK, bulkData.mesh_meta_data().globally_shared_part());
    for(size_t i=0; i<shared_node_buckets.size(); ++i)
    {
        const stk::mesh::Bucket& bucket = *shared_node_buckets[i];
        for(size_t node_index=0; node_index<bucket.size(); ++node_index)
        {
            stk::mesh::Entity node = bucket[node_index];
            const stk::mesh::Entity* elements = bulkData.begin_elements(node);
            unsigned num_elements = bulkData.num_elements(node);
            for(unsigned element_index=0; element_index<num_elements; ++element_index)
            {
                if (bulkData.bucket(elements[element_index]).owned())
                {
                    element_set.insert(elements[element_index]);
                }
            }
        }
    }
    elements_to_communicate.assign(element_set.begin(), element_set.end());
    return elements_to_communicate;
}

void pack_shared_side_nodes_of_elements(stk::CommSparse& comm, const stk::mesh::BulkData& bulkData, const stk::mesh::EntityVector& elements_to_communicate)
{
    for(size_t element_index=0; element_index<elements_to_communicate.size(); ++element_index)
    {
        stk::mesh::Entity elem = elements_to_communicate[element_index];
        stk::mesh::EntityId element_id = bulkData.identifier(elem);

        stk::topology topology = bulkData.bucket(elem).topology();
        const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(elem);
        unsigned num_sides = topology.num_sides();
        for(unsigned side_index=0; side_index<num_sides; ++side_index)
        {
            unsigned num_nodes_this_side = topology.side_topology(side_index).num_nodes();
            stk::mesh::EntityVector side_nodes(num_nodes_this_side);
            topology.side_nodes(elem_nodes, side_index, side_nodes.begin());

            std::vector<stk::mesh::EntityKey> side_node_entity_keys(num_nodes_this_side);
            for(size_t i=0; i<num_nodes_this_side; ++i)
            {
                side_node_entity_keys[i] = bulkData.entity_key(side_nodes[i]);
            }

            std::vector<int> sharing_procs;
            bulkData.shared_procs_intersection(side_node_entity_keys, sharing_procs);

            for(size_t proc_index=0; proc_index<sharing_procs.size(); ++proc_index)
            {
                comm.send_buffer(sharing_procs[proc_index]).pack<stk::mesh::EntityId>(element_id);
                comm.send_buffer(sharing_procs[proc_index]).pack<unsigned>(side_index);
                comm.send_buffer(sharing_procs[proc_index]).pack<unsigned>(num_nodes_this_side);
                for(size_t i=0; i<num_nodes_this_side; ++i)
                {
                    comm.send_buffer(sharing_procs[proc_index]).pack<stk::mesh::EntityKey>(side_node_entity_keys[i]);
                }
            }
        }
    }
}

//BeginDocExample2
void add_possibly_connected_elements_to_graph_using_side_nodes(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, const stk::mesh::EntityVector& side_nodes, ParallelGraphInfo& parallel_graph_info,
        LocalId other_element, int other_side, int other_proc)
{
    stk::mesh::EntityVector elements;
    unsigned num_side_nodes = side_nodes.size();
    stk::mesh::impl::find_locally_owned_elements_these_nodes_have_in_common(bulkData, num_side_nodes, side_nodes.data(), elements);
    for(size_t element_index=0; element_index<elements.size(); ++element_index)
    {
        stk::mesh::Entity elem = elements[element_index];
        stk::topology topology = bulkData.bucket(elem).topology();
        const stk::mesh::Entity* elem_nodes = bulkData.begin_nodes(elem);
        int num_sides = topology.num_sides();
        for(int side_index=0; side_index<num_sides; ++side_index)
        {
            unsigned num_nodes_this_side = topology.side_topology(side_index).num_nodes();
            if (num_nodes_this_side == num_side_nodes)
            {
                stk::mesh::EntityVector side_nodes_this_side(num_nodes_this_side);
                topology.side_nodes(elem_nodes, side_index, side_nodes_this_side.begin());

                std::pair<bool,unsigned> result = topology.side_topology(side_index).equivalent(side_nodes_this_side, side_nodes);
                if (result.first == true)
                {
                    LocalId local_elem_id = bulkData.local_id(elem);
                    elem_graph[local_elem_id].push_back(-1*other_element);
                    via_sides[local_elem_id].push_back(side_index);

                    parallel_graph_info.insert(std::make_pair(std::make_pair(local_elem_id, other_element), parallel_info(other_proc, other_side, result.second)));

                    break;
                }
            }
        }
    }
}
//EndDocExample2

void fill_parallel_graph(const stk::mesh::BulkData& bulkData, ElementGraph& elem_graph,
        SidesForElementGraph& via_sides, ParallelGraphInfo& parallel_graph_info)
{
    stk::mesh::EntityVector elements_to_communicate = get_elements_to_communicate(bulkData);

    stk::CommSparse comm(bulkData.parallel());

    pack_shared_side_nodes_of_elements(comm, bulkData, elements_to_communicate);
    comm.allocate_buffers();

    pack_shared_side_nodes_of_elements(comm, bulkData, elements_to_communicate);
    comm.communicate();

    for(int proc_id=0; proc_id<bulkData.parallel_size(); ++proc_id)
    {
        if (proc_id != bulkData.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                stk::mesh::EntityId element_id;
                comm.recv_buffer(proc_id).unpack<stk::mesh::EntityId>(element_id);
                unsigned side_index = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(side_index);
                unsigned num_side_nodes = 0;
                comm.recv_buffer(proc_id).unpack<unsigned>(num_side_nodes);
                stk::mesh::EntityVector side_nodes(num_side_nodes);
                for(unsigned i=0; i<num_side_nodes; ++i)
                {
                    stk::mesh::EntityKey key;
                    comm.recv_buffer(proc_id).unpack<stk::mesh::EntityKey>(key);
                    side_nodes[i] = bulkData.get_entity(key);
                }

                add_possibly_connected_elements_to_graph_using_side_nodes(bulkData, elem_graph, via_sides, side_nodes,
                        parallel_graph_info, element_id, side_index, proc_id);
            }
        }
    }
}

}
