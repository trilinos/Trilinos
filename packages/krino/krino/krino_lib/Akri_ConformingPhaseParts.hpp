#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CONFORMINGPHASEPARTS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CONFORMINGPHASEPARTS_HPP_

#include <stk_mesh/base/Types.hpp>

namespace krino {

class Phase_Support;
class PhaseTag;

void determine_original_undecomposed_part_changes_for_current_parts(
    const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentParts,
    const stk::mesh::EntityRank entityRank,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts);

void determine_original_undecomposed_part_changes_for_entities(
    const stk::mesh::BulkData & mesh,
    const stk::mesh::Bucket & bucket,
    const Phase_Support & phaseSupport,
    stk::mesh::Part & childPart,
    stk::mesh::Part & parentPart,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts);

void append_conforming_part_changes_for_current_parts(const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentParts,
    const stk::mesh::EntityRank entityRank,
    const PhaseTag & phase,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts);

void append_nonconforming_part_changes_for_current_parts(const Phase_Support & phaseSupport,
    const stk::mesh::PartVector & currentParts,
    const stk::mesh::EntityRank entityRank,
    stk::mesh::Part & childPart,
    stk::mesh::Part & parentPart,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts);

void determine_part_changes_for_side(
    const Phase_Support & phaseSupport,
    const stk::mesh::EntityRank sideRank,
    const stk::mesh::PartVector & currentSideParts,
    const bool doesSideHaveActiveElements,
    stk::mesh::Part & blockBoundarySide,
    stk::mesh::Part & childPart,
    stk::mesh::Part & parentPart,
    stk::mesh::Part & activePart,
    stk::mesh::PartVector & addParts,
    stk::mesh::PartVector & removeParts);

void determine_child_conforming_parts(const stk::mesh::MetaData & meta,
    const Phase_Support & phaseSupport,
    const stk::topology topology,
    const stk::mesh::PartVector & parentParts,
    const stk::mesh::PartVector & attributeParts,
    stk::mesh::Part & childPart,
    stk::mesh::Part & activePart,
    const PhaseTag & phase,
    stk::mesh::PartVector & childParts);
}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CONFORMINGPHASEPARTS_HPP_ */
