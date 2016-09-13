#ifndef STKIOUTILS_HPP_
#define STKIOUTILS_HPP_

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"

namespace stk {
namespace io {

inline
size_t get_entities(stk::mesh::Part &part,
                    stk::mesh::EntityRank type,
                    const stk::mesh::BulkData &bulk,
                    stk::mesh::EntityVector &entities,
                    bool include_shared,
                    const stk::mesh::Selector *subset_selector)
{
    stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);

    stk::mesh::Selector own_share = meta.locally_owned_part();
    if(include_shared)
        own_share |= meta.globally_shared_part();

    stk::mesh::Selector selector = part & own_share;
    if(subset_selector)
        selector &= *subset_selector;

    get_selected_entities(selector, bulk.buckets(type), entities);
    return entities.size();
}

inline
stk::mesh::EntityRank part_primary_entity_rank(const stk::mesh::Part &part)
{
  if (stk::mesh::MetaData::get(part).universal_part() == part) {
    return stk::topology::NODE_RANK;
  }
  else {
    return part.primary_entity_rank();
  }
}

}}


#endif /* STKIOUTILS_HPP_ */
