#ifndef STK_STK_MESH_STK_MESH_BASE_SIDESETUTIL_HPP_
#define STK_STK_MESH_STK_MESH_BASE_SIDESETUTIL_HPP_

#include <stddef.h>                                       // for size_t
#include <map>                                            // for map, etc
#include <string>                                         // for string, etc
#include <vector>                                         // for vector
#include "stk_mesh/base/Types.hpp"                        // for EntityRank, etc
#include "stk_mesh/base/ExodusTranslator.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class SideSet; } }

namespace stk {
namespace mesh {

void reconstruct_sideset(stk::mesh::BulkData& bulkData, const stk::mesh::Part& surfacePart);

void fill_sideset(const stk::mesh::Part& sidesetPart, stk::mesh::BulkData& bulkData, const stk::mesh::Selector& elementSelector);

void create_bulkdata_sidesets(stk::mesh::BulkData& bulkData);

bool should_reconstruct_sideset(const stk::mesh::BulkData& bulkData, const stk::mesh::Part& surfacePart);

bool does_not_contain_internal_sideset(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &sides, const stk::mesh::impl::ParallelPartInfo &parallelPartInfo);

const stk::mesh::Part& get_sideset_parent(const stk::mesh::Part& sidesetPart);

std::pair<bool,bool> is_positive_sideset_polarity(const stk::mesh::BulkData &bulk,
                                                  const stk::mesh::Part& sideSetPart,
                                                  stk::mesh::Entity face,
                                                  const stk::mesh::Part* activePart = nullptr,
                                                  const stk::mesh::SideSet* inputSidesetPtr = nullptr);

std::pair<bool,bool> is_positive_sideset_face_polarity(const stk::mesh::BulkData &bulk, stk::mesh::Entity face,
                                                       const stk::mesh::Part* activePart = nullptr);

std::vector<const stk::mesh::Part*> get_sideset_io_parts(const stk::mesh::BulkData& bulkData, stk::mesh::Entity face);

void toggle_sideset_updaters(stk::mesh::BulkData& bulk, bool flag);

}
}

#endif /* STK_STK_MESH_STK_MESH_BASE_SIDESETUTIL_HPP_ */
