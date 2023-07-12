#ifndef MESH_SNAPPER_H
#define MESH_SNAPPER_H

#include "adjacency_search.hpp"
#include "mesh.hpp"
#include "mesh_snapper_opts.hpp"
#include "predicates/intersection.hpp"
#include <set>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshSnapper
{
  public:
    // TODO: change the constructor to take the meshes, have snap() take zero arguments
    explicit MeshSnapper(const MeshSnapperOpts& opts)
      : m_opts(opts)
    {}

    void snap(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2);

  private:
    template <typename T>
    using SetType = std::set<T, MeshEntityCompare>;

    void snap_el(MeshEntityPtr el1, const std::vector<MeshEntityPtr> mesh2Els);

    MeshSnapperOpts m_opts;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_classifier;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
