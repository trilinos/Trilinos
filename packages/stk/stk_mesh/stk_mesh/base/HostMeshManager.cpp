#include "HostMeshManager.hpp"
#include "BulkData.hpp"
#include "stk_mesh/base/NgpMesh.hpp"

namespace stk {
namespace mesh {

HostMeshManager::HostMeshManager(stk::mesh::BulkData & bulk)
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
#ifdef KOKKOS_ENABLE_CUDA
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
