#include <stk_balance/internal/SubdomainCreator.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <vector>

namespace stk {
namespace balance {
namespace internal {

std::tuple<int,int> get_included_and_num_target_procs(stk::mesh::BulkData &bulk, stk::ParallelMachine comm)
{
    int includeMe = 0;
    {
        if(stk::mesh::count_selected_entities(bulk.mesh_meta_data().locally_owned_part(), bulk.buckets(stk::topology::ELEM_RANK))>0)
            includeMe = 1;
    }

    int numTarget = 0;
    stk::all_reduce_sum(comm, &includeMe, &numTarget, 1);

    return std::make_tuple(includeMe, numTarget);
}

int get_subdomain_index(int includeMe, stk::ParallelMachine comm)
{
    std::vector<int> myData {includeMe};
    std::vector<int> globalData(stk::parallel_machine_size(comm));
    stk::parallel_vector_concat(comm, myData, globalData);
    int mySubdomain = -1;
    if(includeMe)
    {
        mySubdomain = 0;
        for(int i = 0; i < stk::parallel_machine_rank(comm); ++i)
        {
            if(globalData[i] > 0)
                mySubdomain++;
        }
    }
    return mySubdomain;
}

void write_subdomain_files(stk::mesh::BulkData &bulk, int numTarget, int mySubdomain, const std::string& outputMesh)
{
    stk::balance::internal::SubdomainCreator subdomainCreator(bulk, numTarget);
    subdomainCreator.declare_all_subdomain_parts();
    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(bulk.mesh_meta_data().locally_owned_part(),
                                     bulk.buckets(stk::topology::ELEM_RANK),
                                     elements);
    bulk.modification_begin();
    if(mySubdomain >= 0)
        subdomainCreator.move_entities_into_subdomain_part(mySubdomain, elements);
    bulk.modification_end();

    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulk, counts);
    int global_num_nodes = counts[stk::topology::NODE_RANK];
    int global_num_elems = counts[stk::topology::ELEM_RANK];

    if(mySubdomain >= 0)
    {
        stk::io::EntitySharingInfo sharingInfo = subdomainCreator.get_node_sharing_info(mySubdomain);
        subdomainCreator.create_subdomain_and_write(outputMesh, mySubdomain, global_num_nodes, global_num_elems, sharingInfo);
    }
}

}}}


