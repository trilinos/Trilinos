#ifndef FACE_CONTAINER
#define FACE_CONTAINER

#include <stk_middle_mesh/element_operations_2d.hpp>
#include <stk_middle_mesh/mesh.hpp>

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

class FaceContainer
{
  public:
    explicit FaceContainer(const double eps = 1e-13)
      : m_eps(eps)
    {}

    void set_projection(utils::impl::Projection* proj) { m_elemOps.set_projection(proj); }

    bool contains(mesh::MeshEntityPtr face, const utils::Point& pt)
    {
      utils::Point ptXi = m_elemOps.compute_xi_coords(face, pt);
      if (face->get_type() == mesh::MeshEntityType::Triangle)
      {
        double xi3 = 1 - ptXi.get_x() - ptXi.get_y();
        bool val   = ptXi.get_x() >= -m_eps && ptXi.get_x() <= 1 + m_eps && ptXi.get_y() >= -m_eps &&
                   ptXi.get_y() <= 1 + m_eps && xi3 >= -m_eps && xi3 <= 1 + m_eps;
        return val;
      } else // quad
      {
        return ptXi.get_x() >= -m_eps && ptXi.get_y() >= -m_eps && ptXi.get_x() <= 1 + m_eps &&
               ptXi.get_y() <= 1 + m_eps;
      }
    }

    bool contains(mesh::MeshEntityPtr face, mesh::MeshEntityPtr v1) { return contains(face, v1->get_point_orig(0)); }

    double get_eps() const { return m_eps; }

    void set_eps(const double val) { m_eps = val; }

  private:
    mesh::impl::ElementOperations2D m_elemOps;
    double m_eps;
};

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
