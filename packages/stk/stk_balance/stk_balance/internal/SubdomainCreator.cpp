#include "SubdomainCreator.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>

namespace stk {
namespace balance {
namespace internal {

std::string SubdomainCreator::getSubdomainPartName(int subdomainId)
{
    std::ostringstream out;
    out << "subdomain_" << subdomainId;
    return out.str();
}

void SubdomainCreator::declare_all_subdomain_parts()
{
    for(int i=0;i<mNumSubdomains;++i)
    {
        std::string partNameForSubdomain = getSubdomainPartName(i);
        mMeta.declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK);
    }
}

stk::mesh::Part* SubdomainCreator::get_subdomain_part(size_t subdomain_num)
{
    std::string partNameForSubdomain = getSubdomainPartName(subdomain_num);
    return mMeta.get_part(partNameForSubdomain);
}

void SubdomainCreator::move_entities_into_subdomain_part(size_t i, const stk::mesh::EntityVector &entities)
{
    stk::mesh::PartVector partVector = get_parts_to_add_for_subdomain(i);
    for(size_t j = 0; j < entities.size(); j++)
        mBulk.change_entity_parts(entities[j], partVector);
}

stk::mesh::PartVector SubdomainCreator::get_parts_to_add_for_subdomain(size_t subdomain_num)
{
    stk::mesh::Part* part = get_subdomain_part(subdomain_num);
    return {part};
}

stk::io::EntitySharingInfo SubdomainCreator::get_node_sharing_info(unsigned subdomain)
{
    stk::mesh::EntityVector sharedNodes;
    std::vector<int> procsForSharedNodes;
    fill_shared_node_info_for_this_subdomain(subdomain, sharedNodes, procsForSharedNodes);

    stk::io::EntitySharingInfo nodeSharingInfo;
    for(size_t nodeIndex = 0; nodeIndex < sharedNodes.size(); nodeIndex++)
    {
        stk::mesh::EntityId nodeId = mBulk.identifier(sharedNodes[nodeIndex]);
        nodeSharingInfo.push_back({nodeId, procsForSharedNodes[nodeIndex]});
    }
    return nodeSharingInfo;
}

stk::mesh::EntityVector SubdomainCreator::get_nodes_shared_between_subdomains(int this_subdomain_index, int other_subdomain_index)
{
    stk::mesh::Selector selected_nodes = *get_subdomain_part(this_subdomain_index) &
                                         *get_subdomain_part(other_subdomain_index);
    stk::mesh::EntityVector nodes;
    stk::mesh::get_selected_entities(selected_nodes, mBulk.buckets(stk::topology::NODE_RANK), nodes);
    return nodes;
}

void SubdomainCreator::fill_shared_node_proc_info(stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes, int this_subdomain_num, int other_subdomain_num)
{
    stk::mesh::EntityVector nodes = get_nodes_shared_between_subdomains(this_subdomain_num, other_subdomain_num);
    shared_nodes.insert(shared_nodes.end(), nodes.begin(), nodes.end());
    procs_for_shared_nodes.resize(shared_nodes.size(), other_subdomain_num);
}

void SubdomainCreator::fill_shared_node_info_for_this_subdomain(const unsigned this_subdomain_num, stk::mesh::EntityVector& shared_nodes, std::vector<int>& procs_for_shared_nodes)
{
    for(int other_subdomain_num=0;other_subdomain_num<mNumSubdomains;++other_subdomain_num)
    {
        if(static_cast<int>(this_subdomain_num) != other_subdomain_num)
            fill_shared_node_proc_info(shared_nodes, procs_for_shared_nodes, this_subdomain_num, other_subdomain_num);
    }
}

void SubdomainCreator::create_subdomain_and_write(const std::string &filename,
                                                unsigned subdomain,
                                                int global_num_nodes,
                                                int global_num_elems,
                                                const stk::io::EntitySharingInfo &nodeSharingInfo,
                                                int numSteps, double timeStep)
{
    stk::mesh::MetaData newMeta;
    stk::mesh::BulkData newBulkData(newMeta, MPI_COMM_SELF);
    stk::tools::copy_mesh(mBulk, *mMeta.get_part(getSubdomainPartName(subdomain)), newBulkData);
    stk::io::write_file_for_subdomain(filename, subdomain, mNumSubdomains, global_num_nodes, global_num_elems, newBulkData, nodeSharingInfo, numSteps, timeStep);
}

}}}
