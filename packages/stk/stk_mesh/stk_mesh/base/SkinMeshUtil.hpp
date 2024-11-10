// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef SKINMESHUTIL_HPP_
#define SKINMESHUTIL_HPP_

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ParallelInfoForGraph.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <vector>

namespace stk {
namespace mesh {

class BulkData;

class SkinMeshUtil {
public:
    SkinMeshUtil(ElemElemGraph& elemElemGraph,
                 const stk::mesh::Selector& inputSkinSelector,
                 const stk::mesh::Selector* inputAirSelector = nullptr);

    static std::vector<SideSetEntry> get_skinned_sideset(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector);
    static std::vector<SideSetEntry> get_skinned_sideset_excluding_region(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector, const stk::mesh::Selector& exclusionRegionSelector);
    static std::vector<SideSetEntry> get_interior_sideset(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector);
    static std::vector<SideSetEntry> get_all_sides_sideset(stk::mesh::BulkData & bulk, const stk::mesh::Selector& skinSelector, bool includeAuraElementSides = false);

private:
    SkinMeshUtil();

    std::vector<SideSetEntry> extract_skinned_sideset();

    std::vector<SideSetEntry> extract_interior_sideset();

    std::vector<SideSetEntry> extract_all_sides_sideset(bool includeAuraElementSides = false);

    void get_exposed_sides(stk::mesh::impl::LocalId localId, int maxSidesThisElement,
                           std::vector<int>& exposedSides);

    void add_exposed_sides_due_to_air_selector(impl::LocalId local_id, std::vector<int> &exposedSides);

    bool is_remote_element_air(const ParallelInfoForGraphEdges &parallelInfoForGraphEdges, const stk::mesh::GraphEdge &graphEdge);

    bool is_connected_element_air(const stk::mesh::GraphEdge &graphEdge);

    bool is_element_selected_and_can_have_side(const stk::mesh::Selector& selector, stk::mesh::Entity otherElement);

    void get_sides_exposed_on_other_procs(stk::mesh::impl::LocalId localId,
                                          int numElemSides,
                                          std::vector<int>& exposedSides);

    void get_sides_for_skinning(const stk::mesh::Bucket& bucket,
                                stk::mesh::Entity element,
                                stk::mesh::impl::LocalId localId,
                                std::vector<int>& exposedSides);

    void mark_local_connections(const stk::mesh::GraphEdge &graphEdge,
                                std::vector<int> &isOnlyConnectedRemotely);

    void mark_remote_connections(const stk::mesh::GraphEdge &graphEdge,
                                 std::vector<int> &isConnectedToRemoteElementInBodyToSkin);

    void mark_sides_exposed_on_other_procs(const stk::mesh::GraphEdge &graphEdge,
                                           std::vector<int> &isConnectedToRemoteElementInBodyToSkin,
                                           std::vector<int> &isOnlyConnectedRemotely);

    stk::mesh::ElemElemGraph& eeGraph;
    stk::mesh::Selector skinSelector;
    stk::mesh::impl::ParallelSelectedInfo remoteSkinSelector;
    const bool useAirSelector = false;
    stk::mesh::Selector airSelector;
    stk::mesh::impl::ParallelSelectedInfo remoteAirSelector;
    std::vector<int> m_isConnectedToRemoteElem;
    std::vector<int> m_isOnlyConnectedRemotely;
};

}
}



#endif /* SKINMESHUTIL_HPP_ */
