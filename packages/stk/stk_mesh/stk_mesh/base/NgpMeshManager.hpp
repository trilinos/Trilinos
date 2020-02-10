#ifndef STK_MESH_NGPMESHMANAGER_HPP
#define STK_MESH_NGPMESHMANAGER_HPP

#include "stk_mesh/base/Ngp.hpp"

namespace stk {
namespace mesh {

class NgpMeshManager
{
public:
  NgpMeshManager() = default;
  virtual ~NgpMeshManager() = default;

  virtual stk::mesh::NgpMesh & get_mesh() = 0;

  virtual void update_mesh() = 0;
};

}
}

#endif
