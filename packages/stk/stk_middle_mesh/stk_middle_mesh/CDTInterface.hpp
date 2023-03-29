#ifndef CDT_INTERFACE_H
#define CDT_INTERFACE_H

#include "abstract_cdt_interface.hpp"
#include "mesh.hpp"
#include "mesh_io.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class CDTInterface : public AbstractCDTInterface
{
  public:
    explicit CDTInterface(std::shared_ptr<mesh::Mesh> mesh)
      : m_mesh(mesh)
    {}

    void triangulate(const utils::impl::Projection& proj) override;

  private:
    std::shared_ptr<mesh::Mesh> m_mesh;
    const bool m_output = false;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
