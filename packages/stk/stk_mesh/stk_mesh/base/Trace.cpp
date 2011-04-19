#include <stk_mesh/base/Trace.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <set>

namespace stk {
namespace mesh {

//We used to store a static vector 'watch_vector', but its
//contents (pointers) didn't get deleted which resulted in
//memory leaks.
//Instead, we now wrap the watch_vector in this 'WatchVectorHolder'
//struct, which has a destructor which cleans up the contents
//of watch_vector. This eliminates the memory leaks.
struct WatchVectorHolder {
 WatchVectorHolder() : watch_vector() {}
 ~WatchVectorHolder()
  {
    for(std::vector<Watch*>::iterator it=watch_vector.begin();
        it!=watch_vector.end(); ++it) {
      delete *it;
    }
  }

 std::vector<Watch*> watch_vector;
};

std::vector<Watch*>& watch_vector()
{
  static WatchVectorHolder watch_vector_holder;
  return watch_vector_holder.watch_vector;
}

} // namespace mesh
} // namespace stk
