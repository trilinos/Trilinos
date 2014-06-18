#ifndef STK_SIERRA_MESH_ARRAYMESH_ARRAYUTILS_HPP
#define STK_SIERRA_MESH_ARRAYMESH_ARRAYUTILS_HPP

#include <vector>
#include <algorithm>

namespace sierra {
namespace mesh {

/** insert item into vec, if item is not already present, keeping vec sorted.
*/
inline
void insert_sorted(std::vector<int>& vec, int item)
{
  std::vector<int>::iterator it = std::lower_bound(vec.begin(), vec.end(), item);
  if (it == vec.end() || *it != item) {
    vec.insert(it, item);
  }
}

}//namespace mesh
}//namespace sierra

#endif

