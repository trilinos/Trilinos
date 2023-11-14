#ifndef ACTIVE_VERT_DATA_H
#define ACTIVE_VERT_DATA_H

#include <array>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "mat2x2.hpp"
#include "matrix.hpp"
#include "mesh.hpp"
#include "plane_projection.hpp"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

using stk::middle_mesh::mesh::MeshEntityPtr;

class ActiveVertData
{
  public:
    explicit ActiveVertData(std::shared_ptr<mesh::Mesh> mesh, MeshEntityPtr vert) :
      m_mesh(mesh)
    {
      set_sizes(vert);
      set_triangles(vert);
      set_current_indices(vert);
    }

    // gives number of (unique) vertices
    int get_num_verts() const { return m_allPoints.size(); }

    // gives number of unique verts that are present (but not necessarily owned) on this process
    int get_num_local_verts() const { return m_localVertsUnique.size(); }

    // gives number of elements
    int get_num_elements() const { return m_triVertIndices.extent0(); };


    // returns the verts of the ith element
    std::array<stk::middle_mesh::utils::Point, 3> get_element_verts(const int tri)
    {
      assert(tri >= 0 && tri < get_num_elements());
      for (int i=0; i < 3; ++i)
      {
        int idx = m_triVertIndices(tri, i);
        if (idx < int(m_localVertsUnique.size()))
          m_allPoints[idx] = m_localVertsUnique[idx]->get_point_orig(0);
      }

      return {m_allPoints[m_triVertIndices(tri, 0)], 
              m_allPoints[m_triVertIndices(tri, 1)],
              m_allPoints[m_triVertIndices(tri, 2)]};
    }

    std::array<int, 3> get_element_vert_ids(int tri) const 
    { 
      assert(tri >= 0 && tri < get_num_elements());
      return { m_triVertIndices(tri, 0), m_triVertIndices(tri, 1), m_triVertIndices(tri, 2)};
    }

    std::array<stk::middle_mesh::utils::Point, 3> get_element_verts_orig(const int tri)
    {
      assert(tri >= 0 && tri < get_num_elements());
      return {m_allPointsOrig[m_triVertIndices(tri, 0)], m_allPointsOrig[m_triVertIndices(tri, 1)], m_allPointsOrig[m_triVertIndices(tri, 2)]};
    }

    // returns the vert that this patch is centered on
    MeshEntityPtr get_current_vert() const { return m_localVertsUnique[m_triVertIndices(0, m_currentEntityIdx[0])]; }

    stk::middle_mesh::utils::Point get_current_point_orig() const
    {
      return m_allPointsOrig[m_triVertIndices(0, m_currentEntityIdx[0])];
    }

    int get_current_vert_idx(const int tri)
    {
      assert(tri >= 0 && tri < get_num_elements());
      return m_currentEntityIdx[tri];
    }

    std::vector<utils::Point>& get_unique_verts()
    {
      for (size_t i=0; i < m_localVertsUnique.size(); ++i)
        m_allPoints[i] = m_localVertsUnique[i]->get_point_orig(0);
      return m_allPoints;
    }

    const std::vector<MeshEntityPtr>& get_local_verts() const { return m_localVertsUnique; }

    std::vector<MeshEntityPtr>& get_local_verts() { return m_localVertsUnique; }

    // returns the coordinates of the unique verts when the patch
    // was first created
    const std::vector<stk::middle_mesh::utils::Point>& get_points_orig() const { return m_allPointsOrig; }

    mesh::RemoteSharedEntity get_vert_owner(int vertIdx)
    {
      if (vertIdx < int(m_localVertsUnique.size()))
        return mesh::get_owner_remote(m_mesh, m_localVertsUnique[vertIdx]);
      else
        return m_remoteVertOwners[vertIdx - m_localVertsUnique.size()];
    }

    void add_local_vert(mesh::MeshEntityPtr vert);

    void add_remote_vert(const mesh::RemoteSharedEntity& owner);

    void add_remote_element(const std::array<int, 3>& triVertIndices);

    void finish_initialization()
    {
      m_allPointsOrig = m_allPoints;
    }

  private:
    using LocalInt = uint_least16_t; // int size used for indexing
                                     // into local arrays

    void set_sizes(MeshEntityPtr vert);

    void set_triangles(MeshEntityPtr vert);

    // gets the local index of vert in each triangle and writes it to
    // m_current_entity
    void set_current_indices(MeshEntityPtr vert);

    int count_num_triangles(MeshEntityPtr vert, MeshEntityPtr el);

    // gets the local index of vert within el's verts
    int get_local_idx(MeshEntityPtr el, MeshEntityPtr vert);

    // returns the vertices of the triangles in tri_verts, and returns
    // the number of triangles.  If there is only 1 triangle, only
    // the first 3 entries of tri_verts are populated
    int get_triangle_verts(MeshEntityPtr el, MeshEntityPtr vert, std::array<MeshEntityPtr, 12>& triVerts);

    std::shared_ptr<mesh::Mesh> m_mesh;
    std::vector<MeshEntityPtr> m_localVertsUnique;
    std::vector<utils::Point> m_allPoints;
    std::vector<utils::Point> m_allPointsOrig;
    std::vector<mesh::RemoteSharedEntity> m_remoteVertOwners;
    utils::impl::Matrix<LocalInt> m_triVertIndices; // nelem x 3 array, indexes into
                                              // m_verts_unique to give the
                                              // vertices of the triangle
    std::vector<LocalInt> m_currentEntityIdx; // gives index of current
                                              // entity in each element

    const static std::array<int, 3> M_VERTMAP_QUAD_TO_TRI1;
    const static std::array<int, 3> M_VERTMAP_QUAD_TO_TRI2;
    const static std::array<int, 3> M_VERTMAP_QUAD_TO_TRI3;
    const static std::array<int, 3> M_VERTMAP_QUAD_TO_TRI4;
    const static std::array<int, 4> M_EXCLUDED_TRI;
};

} // namespace impl

} // namespace opt
} // namespace middle_mesh
} // namespace stk
#endif
