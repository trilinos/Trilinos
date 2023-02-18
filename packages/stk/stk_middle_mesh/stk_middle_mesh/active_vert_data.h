#ifndef ACTIVE_VERT_DATA_H
#define ACTIVE_VERT_DATA_H

#include <array>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "mat2x2.h"
#include "matrix.h"
#include "mesh.h"
#include "plane_projection.h"

namespace stk {
namespace middle_mesh {
namespace opt {
namespace impl {

using stk::middle_mesh::mesh::MeshEntityPtr;

class ActiveVertData
{
  public:
    explicit ActiveVertData(MeshEntityPtr vert)
    {
      set_sizes(vert);
      get_triangles(vert);
      get_current_indices(vert);
      // computeW();
    }

    // gives number of (unique) vertices
    int get_num_verts() const { return m_vertsUnique.size(); }

    // gives number of elements
    int get_num_elements() const { return m_triVerts.extent0(); };

    // returns the verts of the ith element
    std::array<MeshEntityPtr, 3> get_element_verts(const int i)
    {
      assert(i >= 0 && i < get_num_elements());
      return {m_vertsUnique[m_triVerts(i, 0)], m_vertsUnique[m_triVerts(i, 1)], m_vertsUnique[m_triVerts(i, 2)]};
    }

    std::array<stk::middle_mesh::utils::Point, 3> get_element_verts_orig(const int i)
    {
      assert(i >= 0 && i < get_num_elements());
      return {m_pointsOrig[m_triVerts(i, 0)], m_pointsOrig[m_triVerts(i, 1)], m_pointsOrig[m_triVerts(i, 2)]};
    }

    // returns the vert that this patch is centered on
    MeshEntityPtr get_current_vert() const { return m_vertsUnique[m_triVerts(0, m_currentEntity[0])]; }

    stk::middle_mesh::utils::Point get_current_point_orig() const
    {
      return m_pointsOrig[m_triVerts(0, m_currentEntity[0])];
    }

    int get_current_vert_idx(const int i)
    {
      assert(i >= 0 && i < get_num_elements());
      return m_currentEntity[i];
    }
    /*
        utils::impl::Mat2x2<double> getW(const int el, utils::impl::PlaneProjection val)
        {
          assert(el >= 0 && el < get_num_elements());
          return m_w_mats(el, static_cast<int>(val));
        }
    */
    const std::vector<MeshEntityPtr>& get_unique_verts() const { return m_vertsUnique; }

    // returns the coordinates of the unique verts when the patch
    // was first created
    const std::vector<stk::middle_mesh::utils::Point>& get_points_orig() const { return m_pointsOrig; }

  private:
    using LocalInt = uint_least16_t; // int size used for indexing
                                     // into local arrays

    void set_sizes(MeshEntityPtr vert);

    void get_triangles(MeshEntityPtr vert);

    void compute_w();

    // gets the local index of vert in each triangle and writes it to
    // m_current_entity
    void get_current_indices(MeshEntityPtr vert);

    int count_num_triangles(MeshEntityPtr vert, MeshEntityPtr el);

    // gets the local index of vert within el's verts
    int get_local_idx(MeshEntityPtr el, MeshEntityPtr vert);

    // returns the vertices of the triangles in tri_verts, and returns
    // the number of triangles.  If there is only 1 triangle, only
    // the first 3 entries of tri_verts are populated
    int get_triangle_verts(MeshEntityPtr el, MeshEntityPtr vert, std::array<MeshEntityPtr, 12>& triVerts);

    std::vector<MeshEntityPtr> m_vertsUnique;
    std::vector<utils::Point> m_pointsOrig;
    utils::impl::Matrix<LocalInt> m_triVerts; // nelem x 3 array, indexes into
                                              // m_verts_unique to give the
                                              // vertices of the triangle
    std::vector<LocalInt> m_currentEntity;    // gives index of current
                                              // entity in each element
    //    utils::impl::Matrix< utils::impl::Mat2x2<double> > m_w_mats;  // nelem x 6 containing the W
    //                                        // matrix for all 6 possible
    //                                        // projections
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
