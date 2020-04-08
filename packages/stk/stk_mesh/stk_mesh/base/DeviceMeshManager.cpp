#include "DeviceMeshManager.hpp"
#include "BulkData.hpp"
#include "NgpMesh.hpp"
#include "stk_util/stk_kokkos_macros.h"

namespace stk {
namespace mesh {

DeviceMeshManager::DeviceMeshManager(const stk::mesh::BulkData & bulk)
  : NgpMeshManager(),
    m_deviceMesh(bulk)
{
}

DeviceMeshManager::~DeviceMeshManager()
{
}

stk::mesh::NgpMesh &
DeviceMeshManager::get_mesh()
{
#ifndef STK_USE_DEVICE_MESH
  ThrowErrorMsg("DeviceMeshManager does not have host mesh");
  static HostMesh hostMesh;
  return hostMesh;  // Never used; Avoids build problems on the CPU
#else
  return m_deviceMesh;
#endif
}

void DeviceMeshManager::update_mesh()
{
  m_deviceMesh.update_mesh();
}

}
}
