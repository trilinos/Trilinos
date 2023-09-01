#ifndef STK_MIDDLE_MESH_MESH_SCATTER_SPEC_H
#define STK_MIDDLE_MESH_MESH_SCATTER_SPEC_H

#include "mesh.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "utils.hpp"


namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshScatterSpec {
public:
  explicit MeshScatterSpec(MPI_Comm comm, std::shared_ptr<mesh::Mesh> mesh) :
    m_comm(comm),
    m_mesh(mesh)
  {}

  MPI_Comm get_comm() const { return m_comm; }

  std::shared_ptr<mesh::Mesh> get_mesh() const { return m_mesh; }

  void add_destination(mesh::MeshEntityPtr entity, int destRank);

  void get_destinations(mesh::MeshEntityPtr entity, std::vector<int>& destRanks);

private:
  using MapKey = std::pair<int, mesh::MeshEntityType>;
  using MapValue = std::vector<int>;

  MPI_Comm m_comm;
  std::shared_ptr<mesh::Mesh> m_mesh;
  std::map<MapKey, MapValue> m_entityProcMap;
};


}
}
}
}

#endif
