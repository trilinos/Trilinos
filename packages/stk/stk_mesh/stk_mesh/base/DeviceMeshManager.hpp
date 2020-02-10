#ifndef DEVICEMESHMANAGER_HPP
#define DEVICEMESHMANAGER_HPP

#include "stk_mesh/base/NgpMeshManager.hpp"
#include "stk_mesh/base/NgpMesh.hpp"

namespace stk {
namespace mesh {

class BulkData;

class DeviceMeshManager : public NgpMeshManager
{
public:
  DeviceMeshManager(stk::mesh::BulkData & bulk);
  ~DeviceMeshManager() override;

  stk::mesh::NgpMesh & get_mesh() override;

  void update_mesh() override;

private:
  stk::mesh::DeviceMesh m_deviceMesh;
};

}
}

#endif
