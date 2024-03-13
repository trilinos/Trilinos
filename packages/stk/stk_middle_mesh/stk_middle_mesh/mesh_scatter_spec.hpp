#ifndef STK_MIDDLE_MESH_MESH_SCATTER_SPEC_H
#define STK_MIDDLE_MESH_MESH_SCATTER_SPEC_H

#include "mesh.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "utils.hpp"
#include "variable_size_field.hpp"


namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshScatterSpec {
public:
  explicit MeshScatterSpec(MPI_Comm comm, std::shared_ptr<mesh::Mesh> mesh) :
    m_comm(comm),
    m_mesh(mesh),
    m_destRanks(mesh ? create_variable_size_field<int>(mesh, FieldShape(0, 0, 1)) : nullptr)
  {}

  MPI_Comm get_comm() const { return m_comm; }

  std::shared_ptr<mesh::Mesh> get_mesh() const { return m_mesh; }

  void add_destination(mesh::MeshEntityPtr entity, int destRank);

  void get_destinations(mesh::MeshEntityPtr entity, std::vector<int>& destRanks);

private:
  MPI_Comm m_comm;
  std::shared_ptr<mesh::Mesh> m_mesh;
  VariableSizeFieldPtr<int> m_destRanks;
};


}
}
}
}

#endif
