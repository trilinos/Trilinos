
#ifndef SKINMESHUTIL_HPP_
#define SKINMESHUTIL_HPP_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ParallelInfoForGraph.hpp>
#include <vector>

namespace stk {
namespace mesh {

class SkinMeshUtil {
public:
    SkinMeshUtil(ElemElemGraph& elemElemGraph,
                 const stk::mesh::PartVector& skinParts,
                 const stk::mesh::Selector& inputSkinSelector,
                 const stk::mesh::Selector* inputAirSelector = nullptr);

    std::vector<SideSetEntry> extract_skinned_sideset();

    std::vector<SideSetEntry> extract_interior_sideset();

    std::vector<SideSetEntry> extract_all_sides_sideset();

private:
    SkinMeshUtil();

    std::vector<int> get_exposed_sides(stk::mesh::impl::LocalId localId, int maxSidesThisElement);

    void add_exposed_sides_due_to_air_selector(impl::LocalId local_id, std::vector<int> &exposedSides);

    bool is_remote_element_air(const ParallelInfoForGraphEdges &parallelInfoForGraphEdges, const stk::mesh::GraphEdge &graphEdge);

    bool is_connected_element_air(const stk::mesh::GraphEdge &graphEdge);

    bool is_element_selected_and_can_have_side(const stk::mesh::Selector& selector, stk::mesh::Entity otherElement);

    std::vector<int> get_sides_exposed_on_other_procs(stk::mesh::impl::LocalId localId,
                                                      int numElemSides);

    std::vector<int> get_sides_for_skinning(const stk::mesh::Selector& skinSelector,
                                            const stk::mesh::Bucket& bucket,
                                            stk::mesh::Entity element,
                                            stk::mesh::impl::LocalId localId,
                                            const stk::mesh::Selector* airSelector = nullptr);

    void mark_local_connections(const stk::mesh::GraphEdge &graphEdge,
                                std::vector<bool> &isOnlyConnectedRemotely);

    void mark_remote_connections(const stk::mesh::GraphEdge &graphEdge,
                                 std::vector<bool> &isConnectedToRemoteElementInBodyToSkin);

    void mark_sides_exposed_on_other_procs(const stk::mesh::GraphEdge &graphEdge,
                                           std::vector<bool> &isConnectedToRemoteElementInBodyToSkin,
                                           std::vector<bool> &isOnlyConnectedRemotely);

    stk::mesh::ElemElemGraph& eeGraph;
    const stk::mesh::Selector& skinSelector;
    stk::mesh::impl::ParallelSelectedInfo remoteSkinSelector;
    const stk::mesh::Selector* airSelector;
    stk::mesh::impl::ParallelSelectedInfo remoteAirSelector;
};

}
}



#endif /* SKINMESHUTIL_HPP_ */
