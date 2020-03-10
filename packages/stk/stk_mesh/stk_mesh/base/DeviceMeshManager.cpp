#include "DeviceMeshManager.hpp"
#include "BulkData.hpp"
#include "stk_mesh/base/NgpMesh.hpp"

namespace stk {
namespace mesh {

DeviceMeshManager::DeviceMeshManager(stk::mesh::BulkData & bulk)
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
#ifndef KOKKOS_ENABLE_CUDA
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
