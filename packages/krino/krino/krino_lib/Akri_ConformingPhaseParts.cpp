#include <Akri_ConformingPhaseParts.hpp>

#include <Akri_Phase_Support.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_io/IossBridge.hpp>

namespace krino {

void determine_original_undecomposed_part_changes_for_current_parts(
    const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentParts,
    const stk::mesh::EntityRank entityRank,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  const auto & allDecomposedBlocksSelector = phaseSupport.get_all_decomposed_blocks_selector();
  for(auto * partPtr : currentParts)
  {
    if( partPtr->primary_entity_rank() == entityRank && (stk::io::is_part_io_part(*partPtr) || allDecomposedBlocksSelector(partPtr)) )
    {
      stk::mesh::Part * originalPartPtr = const_cast<stk::mesh::Part *>(phaseSupport.find_original_part(*partPtr));
      if (nullptr != originalPartPtr && originalPartPtr != partPtr)
      {
        addParts.push_back(originalPartPtr);
        removeParts.push_back(partPtr);

        for(auto * supersetPartPtr : partPtr->supersets())
          if (!stk::mesh::is_auto_declared_part(*supersetPartPtr))
            removeParts.push_back(supersetPartPtr);
      }
    }
  }
}

void determine_original_undecomposed_part_changes_for_entities(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Bucket & bucket,
    const Phase_Support & phaseSupport,
    stk::mesh::Part & childPart,
    stk::mesh::Part & parentPart,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts)
{
  addParts.clear();
  removeParts.clear();
  const stk::mesh::PartVector & currentParts = bucket.supersets();
  determine_original_undecomposed_part_changes_for_current_parts(phaseSupport, currentParts, bucket.entity_rank(), addParts, removeParts);

  if (!bucket.member(activePart)) addParts.push_back(&activePart);
  if (bucket.member(childPart)) removeParts.push_back(&childPart);
  if (bucket.member(parentPart)) removeParts.push_back(&parentPart);
}

}
