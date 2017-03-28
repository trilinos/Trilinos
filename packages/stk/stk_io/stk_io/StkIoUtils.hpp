#ifndef STKIOUTILS_HPP_
#define STKIOUTILS_HPP_

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"

namespace stk { namespace io { class StkMeshIoBroker; } };
namespace stk { namespace io { class MetaData; } };
namespace stk { namespace io { class BulkData; } };

namespace stk {
namespace io {

size_t get_entities(stk::mesh::Part &part,
                    stk::mesh::EntityRank type,
                    const stk::mesh::BulkData &bulk,
                    stk::mesh::EntityVector &entities,
                    bool include_shared,
                    const stk::mesh::Selector *subset_selector);

stk::mesh::EntityRank part_primary_entity_rank(const stk::mesh::Part &part);

typedef std::map<std::string, std::vector<std::string>> IossBlockMembership;

IossBlockMembership get_block_memberships(stk::io::StkMeshIoBroker& stkIo);

void fill_block_parts_given_names(const std::vector<std::string>& side_block_names,
                                              stk::mesh::MetaData& meta,
                                              std::vector<const stk::mesh::Part*>& blocks);

void create_bulkdata_sidesets(stk::mesh::BulkData& bulkData);

void clear_bulkdata_sidesets(stk::mesh::BulkData& bulkData);

bool isSidesetSupported(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &sides, const stk::mesh::impl::ParallelPartInfo &parallelPartInfo);

stk::mesh::FieldVector get_transient_fields(stk::mesh::MetaData &meta);
stk::mesh::FieldVector get_transient_fields(stk::mesh::MetaData &meta, const stk::mesh::EntityRank rank);

}}


#endif /* STKIOUTILS_HPP_ */
