#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CONFORMINGPHASEPARTS_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CONFORMINGPHASEPARTS_HPP_

#include <stk_mesh/base/Types.hpp>

namespace krino {

class Phase_Support;

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

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CONFORMINGPHASEPARTS_HPP_ */
