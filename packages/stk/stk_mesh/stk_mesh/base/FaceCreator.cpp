#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/base/SideSetEntry.hpp>
#include <stk_mesh/base/FaceCreator.hpp>

namespace stk {
namespace mesh {

FaceCreator::FaceCreator(stk::mesh::BulkData& bulkData, stk::mesh::ElemElemGraph& elemElemGraph)
: m_bulkData(bulkData), m_eeGraph(elemElemGraph)
{

}

FaceCreator::~FaceCreator()
{

}

void FaceCreator::fill_side_ordinals(size_t element_side_index, const SideSet &skinnedSideSet, std::vector<int>& ordinals)
{
    size_t j=element_side_index;
    do
    {
        ordinals.push_back(skinnedSideSet[j].side);
        j++;
    }
    while(j<skinnedSideSet.size() && skinnedSideSet[j].element==skinnedSideSet[element_side_index].element);
}

std::vector<int> FaceCreator::get_side_ordinals_of_element(size_t element_side_index, const SideSet &skinnedSideSet)
{
    std::vector<int> ordinals;
    const int max_num_sides_per_any_element = 20;
    ordinals.reserve(max_num_sides_per_any_element);
    fill_side_ordinals(element_side_index, skinnedSideSet, ordinals);
    return ordinals;
}

size_t FaceCreator::create_face_entities_per_element(size_t element_side_index, const SideSet &skinnedSideSet, const stk::mesh::PartVector& skinParts, std::vector<stk::mesh::sharing_info> &sharedModified)
{
    std::vector<int> ordinals = get_side_ordinals_of_element(element_side_index, skinnedSideSet);
    stk::mesh::impl::LocalId localId = m_eeGraph.get_local_element_id(skinnedSideSet[element_side_index].element);
    m_eeGraph.create_side_entities(ordinals, localId, skinParts, sharedModified);
    return ordinals.size();
}

void FaceCreator::create_side_entities_given_sideset(const SideSet &skinnedSideSet, const stk::mesh::PartVector& skinParts, bool doLocalModCycle)
{
  std::vector<stk::mesh::sharing_info> sharedModified;
  if (doLocalModCycle) {
    m_bulkData.modification_begin();
  }
  for(size_t i=0;i<skinnedSideSet.size();) {
    i += create_face_entities_per_element(i, skinnedSideSet, skinParts, sharedModified);
  }

  std::sort(sharedModified.begin(), sharedModified.end(), SharingInfoLess());
  if(!sharedModified.empty()) {
    for(size_t i=1; i<sharedModified.size(); i++) {
      if(sharedModified[i].m_entity == sharedModified[i-1].m_entity) {
        sharedModified[i].m_owner = sharedModified[i-1].m_owner;
      }
    }
  }

  if (doLocalModCycle) {
    m_bulkData.make_mesh_parallel_consistent_after_skinning(sharedModified);
  }
}

}
}
