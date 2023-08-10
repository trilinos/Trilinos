#include "stk_mesh/base/BulkData.hpp"
#include <stk_mesh/base/PolarityUtil.hpp>
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/SidesetUpdater.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_mesh/baseImpl/elementGraph/GraphTypes.hpp"
#include <stk_util/environment/RuntimeWarning.hpp>

namespace stk {
namespace mesh {

namespace {

const impl::ParallelPartInfo& get_parallel_part_info(const stk::mesh::BulkData& bulk)
{
  std::vector<std::shared_ptr<SidesetUpdater>> updaters = bulk.get_observer_type<SidesetUpdater>();
  STK_ThrowRequireMsg(updaters.size() > 0, "ERROR: Bulk data does not have a SidesetUpdater enabled");
  return updaters[0]->get_parallel_part_info();
}

}

PolarityUtil::PolarityUtil(const BulkData &bulk, const Selector& activeSelector)
  : m_bulk(bulk)
  , m_activeSelector(activeSelector)
  , m_parallelPartInfo(get_parallel_part_info(bulk))
  , m_polarityDataCache(m_bulk, m_activeSelector, m_parallelPartInfo)
{
  STK_ThrowRequireMsg(bulk.has_face_adjacent_element_graph(), "ERROR: Bulk data does not have face adjacency graph enabled");
}


std::pair<bool,bool> PolarityUtil::is_positive_sideset_polarity(const Part& sideSetPart,
                                                                const Entity face,
                                                                const SideSet* inputSidesetPtr)
{
  return impl::internal_is_positive_sideset_polarity(sideSetPart, face, m_polarityDataCache, inputSidesetPtr);
}

std::pair<bool,bool> PolarityUtil::is_positive_sideset_face_polarity(const Entity face)
{
  return impl::internal_is_positive_sideset_face_polarity(face, m_polarityDataCache);
}

}
}
