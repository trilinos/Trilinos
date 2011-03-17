#include <stk_mesh/base/Trace.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <set>

namespace stk {
namespace mesh {

std::vector<Watch*>& watch_vector()
{
  static std::vector<Watch*> watch_vector;
  return watch_vector;
}

} // namespace mesh
} // namespace stk
