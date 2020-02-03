#ifndef HOSTMESHMANAGER_HPP
#define HOSTMESHMANAGER_HPP

#include "stk_mesh/base/NgpMeshManager.hpp"
#include "stk_mesh/base/NgpMesh.hpp"

namespace stk {
namespace mesh {

class BulkData;

class HostMeshManager : public NgpMeshManager
{
public:
  HostMeshManager(stk::mesh::BulkData & bulk);
  ~HostMeshManager() override;

  stk::mesh::NgpMesh & get_mesh() override;

  void update_mesh() override;

private:
  stk::mesh::HostMesh m_hostMesh;
};

}
}

#endif
