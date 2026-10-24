/*
 * Akri_SmallClusterCleanup.cpp
 *
 *  Created on: Apr 14, 2026
 *      Author: drnoble
 */

#include <Akri_SmallClusterCleanup.hpp>

#include <Akri_CDMesh_Utils.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_SideAttachedElements.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

namespace krino {

void cleanup_small_clusters_of_phase(stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const NamedPhase & fromPhase, const NamedPhase & toPhase, const size_t smallClusterSize)
{
  const stk::mesh::Selector fromSelector = phaseSupport.get_selector_of_blocks_with_phase(fromPhase.tag());
  const std::vector<stk::mesh::Entity> elemsInSmallCluster =
     (smallClusterSize > 0) ?
     find_selected_owned_elements_that_are_in_a_side_attached_group_smaller_than_size(mesh, fromSelector, smallClusterSize) :
     find_owned_elements_that_are_not_in_the_largest_group_of_selected_side_attached_elements(mesh, fromSelector);
  const std::map<int,int> phaseConversionMap = phaseSupport.build_phase_conversion_map(fromPhase.tag(), toPhase.tag());
  krinolog << "Converting " << elemsInSmallCluster.size() << " elements in small clusters of " << fromPhase << " to " << toPhase << "\n";
  if (stk::is_true_on_any_proc(mesh.parallel(), !elemsInSmallCluster.empty()))
    batch_convert_elements_and_their_sides(mesh, phaseConversionMap, elemsInSmallCluster);
}

void cleanup_small_clusters_of_phases(stk::mesh::BulkData & mesh, const Phase_Support & phaseSupport, const std::vector<std::tuple<NamedPhase,NamedPhase,unsigned>> & fromPhaseToPhaseAndClusterSize)
{
  for (const auto & [fromPhase, toPhase, smallClusterSize] : fromPhaseToPhaseAndClusterSize)
    cleanup_small_clusters_of_phase(mesh, phaseSupport, fromPhase, toPhase, smallClusterSize);
}

}


