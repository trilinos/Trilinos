#include "HostMeshManager.hpp"
#include "BulkData.hpp"
#include "NgpMesh.hpp"
#include "stk_util/stk_kokkos_macros.h"

namespace stk {
namespace mesh {

HostMeshManager::HostMeshManager(const stk::mesh::BulkData & bulk)
  : NgpMeshManager(),
    m_hostMesh(bulk)
{
}

HostMeshManager::~HostMeshManager()
{
}

stk::mesh::NgpMesh &
HostMeshManager::get_mesh()
{
#ifdef STK_USE_DEVICE_MESH
  ThrowErrorMsg("HostMeshManager does not have device mesh");
  static DeviceMesh deviceMesh;
  return deviceMesh;  // Never used; Avoids build problems on the GPU
#else
  return m_hostMesh;
#endif
}

void HostMeshManager::update_mesh()
{
  m_hostMesh.update_mesh();
}

}
}
