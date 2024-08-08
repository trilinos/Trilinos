#ifndef STK_STK_MESH_STK_MESH_BASEIMPL_SIDESETUTILIMPL_HPP_
#define STK_STK_MESH_STK_MESH_BASEIMPL_SIDESETUTILIMPL_HPP_

#include <stddef.h>                                       // for size_t
#include <map>                                            // for map, etc
#include <string>                                         // for string, etc
#include <vector>                                         // for vector
#include "stk_mesh/base/Types.hpp"                        // for EntityRank, etc
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class SideSet; } }
namespace stk { namespace mesh { class ElemElemGraph; } }
namespace stk { namespace mesh { struct GraphEdge; } }

namespace stk {
namespace mesh {
namespace impl {

struct PolarityDataCache
{
public:
  PolarityDataCache(const BulkData &bulk_, const Selector& activeSelector_, const ParallelPartInfo& parallelPartInfo_);

  const BulkData &bulk;
  const MetaData& meta;
  const Selector &activeSelector;
  const ParallelPartInfo& parallelPartInfo;
  const EntityRank sideRank;

  std::vector<bool> isConnectedElementActive;
  std::vector<bool> isConnectedElementInSideset;
  std::vector<bool> isConnectedElementPermutationPositive;

private:
  PolarityDataCache() = delete;
};

std::pair<bool,bool> internal_is_positive_sideset_polarity(const Part& sideSetPart,
                                                           const Entity face,
                                                           PolarityDataCache& cache,
                                                           const SideSet* inputSidesetPtr = nullptr);

std::pair<bool,bool> internal_is_positive_sideset_face_polarity(const Entity face, PolarityDataCache& cache);

}
}
}

#endif /* STK_STK_MESH_STK_MESH_BASEIMPL_SIDESETUTILIMPL_HPP_ */
