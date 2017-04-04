#include "ElemElemGraphImpl.hpp"
#include "ElemGraphCoincidentElems.hpp"
#include "ElemElemGraph.hpp"

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
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/util/SortAndUnique.hpp>

namespace stk
{
namespace mesh
{
namespace impl
{

bool are_side_nodes_degenerate(const stk::mesh::EntityVector &sideNodes)
{
    stk::mesh::EntityVector sortedNodes = sideNodes;
    stk::util::sort_and_unique(sortedNodes);

    return sortedNodes.size() != sideNodes.size();
}

struct TopologyChecker
{
    bool are_both_shells() const
    {
        return is_shell_or_beam2(localTopology) && is_shell_or_beam2(remoteTopology);
    }

    bool are_both_not_shells() const
    {
        return !is_shell_or_beam2(localTopology) && !is_shell_or_beam2(remoteTopology);
    }

    stk::topology localTopology;
    stk::topology remoteTopology;
};

bool is_side_node_permutation_positive(const stk::mesh::BulkData &bulkData, stk::mesh::Entity localElem, const stk::mesh::EntityVector& localElemSideNodes, unsigned sideIndex, const stk::mesh::EntityVector &otherElemSideNodes)
{
    EquivAndPositive result = stk::mesh::is_side_equivalent_and_positive(bulkData, localElem, sideIndex, otherElemSideNodes.data());
    return result.is_equiv && result.is_positive;
}

bool is_nondegenerate_coincident_connection(const stk::mesh::BulkData &bulkData,
                                            stk::mesh::Entity localElem,
                                            const stk::mesh::EntityVector& localElemSideNodes,
                                            unsigned sideIndex,
                                            stk::topology otherElemTopology,
                                            const stk::mesh::EntityVector &otherElemSideNodes)
{
    stk::topology localTopology = bulkData.bucket(localElem).topology();
    TopologyChecker topologyChecker {localTopology, otherElemTopology};
    if(topologyChecker.are_both_shells())
        return true;
    if(topologyChecker.are_both_not_shells())
        return is_side_node_permutation_positive(bulkData, localElem, localElemSideNodes, sideIndex, otherElemSideNodes);
    return false;
}

bool is_coincident_connection(const stk::mesh::BulkData &bulkData,
                              stk::mesh::Entity localElem,
                              const stk::mesh::EntityVector& localElemSideNodes,
                              unsigned sideIndex,
                              stk::topology otherElemTopology,
                              const stk::mesh::EntityVector &otherElemSideNodes)
{
    if(are_side_nodes_degenerate(otherElemSideNodes))
        return false;
    return is_nondegenerate_coincident_connection(bulkData, localElem, localElemSideNodes, sideIndex, otherElemTopology, otherElemSideNodes);
}

}}} // end namespaces stk mesh impl

