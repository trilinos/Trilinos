#include "balanceMtoN.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/entityDataToField.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Comm.hpp>
#include "MxNutils.hpp"
#include "MtoNRebalancer.hpp"

namespace stk {
namespace balance {
namespace internal {

void execute_rebalanceMtoN(MtoNRebalancer& m2nRebalancer, const std::string& outputFilename, int numSteps, double timeStep)
{
    stk::mesh::BulkData& bulkData = m2nRebalancer.get_bulk();
    std::vector<unsigned> targetSubdomainsToProc = stk::balance::internal::assign_target_subdomains_roundrobin_to_procs(bulkData.parallel_size(), m2nRebalancer.get_num_target_subdomains());

    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(bulkData, counts);
    int global_num_nodes = counts[stk::topology::NODE_RANK];
    int global_num_elems = counts[stk::topology::ELEM_RANK];

    m2nRebalancer.move_subdomains_such_that_entire_subdomain_doesnt_span_proc_boundaries(targetSubdomainsToProc);

    for(unsigned subdomain = 0; subdomain < targetSubdomainsToProc.size(); subdomain++)
    {
        if(m2nRebalancer.does_this_proc_own_subdomain(targetSubdomainsToProc[subdomain]))
        {
            stk::io::EntitySharingInfo nodeSharingInfo = m2nRebalancer.get_node_sharing_info(subdomain);
            m2nRebalancer.create_subdomain_and_write(outputFilename, subdomain, global_num_nodes, global_num_elems, nodeSharingInfo, numSteps, timeStep);
        }
    }
}

bool rebalanceMtoN(stk::mesh::BulkData& bulkData, stk::mesh::Field<double> &targetDecompField, int num_target_procs, const std::string& outputFilename, int numSteps, double timeStep)
{
    stk::balance::GraphCreationSettings graphSettings;
    MtoNRebalancer m2nRebalancer(bulkData, targetDecompField, graphSettings, num_target_procs);
    execute_rebalanceMtoN(m2nRebalancer, outputFilename, numSteps, timeStep);

    return true;
}

bool rebalanceMtoN(stk::io::StkMeshIoBroker& ioBroker, stk::mesh::Field<double> &targetDecompField, int num_target_procs, const std::string& outputFilename, int numSteps, double timeStep)
{
    stk::balance::GraphCreationSettings graphSettings;
    MtoNRebalancer m2nRebalancer(ioBroker, targetDecompField, graphSettings, num_target_procs);
    execute_rebalanceMtoN(m2nRebalancer, outputFilename, numSteps, timeStep);

    return true;
}
}}}
