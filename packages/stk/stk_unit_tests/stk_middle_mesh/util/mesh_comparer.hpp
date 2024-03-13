#include "stk_middle_mesh/mesh.hpp"

namespace stk {
namespace middle_mesh {
namespace impl {

// Verifies that two meshes have identical mesh data structures.  It uses centroids
// to locate corresponding elements and then verifies the data structures are
// identical (including entity orientations)
class MeshComparer
{
  public:
    MeshComparer(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2, double tol);

    void compare();

  private:
    void compare_entity_counts();

    void compare_centroids();

    mesh::MeshEntityPtr find_closest_mesh2_entity(int dim, const utils::Point& centroid);

    void compare_downward(mesh::MeshEntityPtr mesh1Entity, mesh::MeshEntityPtr mesh2Entity);

    void compare_upward(mesh::MeshEntityPtr mesh1Entity, mesh::MeshEntityPtr mesh2Entity);

    mesh::MeshEntityPtr find_upward_entity(mesh::MeshEntityPtr entity, const utils::Point& centroid);

    double compute_dist(const utils::Point& pt1, const utils::Point& pt2);

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    double m_tol;
};

void compare_meshes(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2, double tol);

} // namespace impl
} // namespace middle_mesh
} // namespace stk
