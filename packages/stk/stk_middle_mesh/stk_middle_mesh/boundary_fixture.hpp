#ifndef BOUNDARY_FIXTURE_H
#define BOUNDARY_FIXTURE_H

#include "field.hpp"
#include "mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// a callable class that returns true if a mesh vertex is *not* on
// the boundary
class BoundaryFixture
{
  public:
    BoundaryFixture(std::shared_ptr<Mesh> mesh) { set_field(mesh); }

    bool operator()(MeshEntityPtr e) { return (*m_field)(e, 0, 0); }

  private:
    using Bool = int_least8_t; // Field<bool> doesn't work

    void set_field(std::shared_ptr<Mesh> mesh)
    {
      m_field     = create_field<Bool>(mesh, FieldShape(1, 0, 0), 1, true);
      auto& field = *m_field;
      for (auto& edge : mesh->get_edges())
        if (edge && edge->count_up() == 1 && edge->count_remote_shared_entities() == 0)
        {
          field(edge->get_down(0), 0, 0) = false;
          field(edge->get_down(1), 0, 0) = false;
        }
    }

    std::shared_ptr<Field<Bool>> m_field;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
