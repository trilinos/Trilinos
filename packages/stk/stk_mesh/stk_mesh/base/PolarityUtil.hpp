#ifndef STK_STK_MESH_STK_MESH_BASE_POLARITYUTIL_HPP_
#define STK_STK_MESH_STK_MESH_BASE_POLARITYUTIL_HPP_

#include <stddef.h>                                       // for size_t
#include <map>                                            // for map, etc
#include <string>                                         // for string, etc
#include <vector>                                         // for vector
#include <memory>
#include "stk_mesh/base/Types.hpp"                        // for EntityRank, etc
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
#include "stk_mesh/baseImpl/SideSetUtilImpl.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class SideSet; } }
namespace stk { namespace mesh { class ElemElemGraph; } }
namespace stk { namespace mesh { struct GraphEdge; } }

namespace stk {
namespace mesh {

class PolarityUtil
{
public:
  PolarityUtil(const BulkData &bulk, const Selector& activeSelector);

  std::pair<bool,bool> is_positive_sideset_polarity(const Part& sideSetPart,
                                                    const Entity face,
                                                    const SideSet* inputSidesetPtr = nullptr);

  std::pair<bool,bool> is_positive_sideset_face_polarity(const Entity face);

protected:
  const BulkData &m_bulk;
  Selector m_activeSelector;
  const impl::ParallelPartInfo& m_parallelPartInfo;
  impl::PolarityDataCache m_polarityDataCache;
};

}
}



#endif /* STK_STK_MESH_STK_MESH_BASE_POLARITYUTIL_HPP_ */
