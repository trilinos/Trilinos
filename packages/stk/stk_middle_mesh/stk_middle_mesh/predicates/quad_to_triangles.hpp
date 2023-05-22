#ifndef PREDICATES_QUAD_TO_TRIANGLES
#define PREDICATES_QUAD_TO_TRIANGLES

#include "intersection_common.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "triangle_coord_utils.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

class QuadToTriangles
{
  public:
    QuadToTriangles()
    {
      m_mesh   = mesh::make_empty_mesh();
      verts[0] = m_mesh->create_vertex(0, 0);
      verts[1] = m_mesh->create_vertex(0, 0);
      verts[2] = m_mesh->create_vertex(0, 0);
      verts[3] = m_mesh->create_vertex(0, 0);
      el1      = m_mesh->create_triangle_from_verts(verts[VERTMAP_TRI1_TO_QUAD[0]], verts[VERTMAP_TRI1_TO_QUAD[1]],
                                                    verts[VERTMAP_TRI1_TO_QUAD[2]]);

      el2 = m_mesh->create_triangle_from_verts(verts[VERTMAP_TRI2_TO_QUAD[0]], verts[VERTMAP_TRI2_TO_QUAD[1]],
                                               verts[VERTMAP_TRI2_TO_QUAD[2]]);
    }

    void set_triangles(mesh::MeshEntityPtr el)
    {
      assert(el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad);
      mesh::MeshEntityPtr elVerts[mesh::MAX_DOWN];
      get_downward(el, 0, elVerts);

      for (int i = 0; i < 4; ++i)
        verts[i]->set_point_orig(0, elVerts[i]->get_point_orig(0));
    }

    // given the PointRecords for the same point on m_el1 and m_el2, figures
    // out the PointRecord for the point on the quad
    // If a point is classified on a vertex or edge of one of the triangles,
    // that triangle will be r1. Ties a settled in favor of
    // m_el1
    PointRecord get_quad_record(mesh::MeshEntityPtr quad, PointRecordForTriangle& r1, PointRecordForTriangle& r2);

    PointRecord create_record(mesh::MeshEntityPtr quad, int edgeId, double edgeXi);

    PointRecord create_record(mesh::MeshEntityPtr quad, int vertId);

    // enforces consistency conditions between the records
    void enforce_record_consistency(impl::PointRecordForTriangle& r1, PointRecordForTriangle& r2);

    utils::Point compute_xyz_coords(const PointRecord& quadRecord, bool allowExterior=false);

    utils::Point get_quad_xi_coords(const PointRecord& quadRecord, bool allowExterior=false);

    PointRecord classify_onto(const PointRecord& record, mesh::MeshEntityPtr el);


    // Data describing the two triangles the quad is broken up into
    mesh::MeshEntityPtr el1;
    mesh::MeshEntityPtr el2;
    mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
    // TODO: make constexpr
    // TODO: use std::array
    const static int VERTMAP_TRI1_TO_QUAD[3];
    const static int VERTMAP_TRI2_TO_QUAD[3];
    const static int EDGEMAP_TRI1_TO_QUAD[3];
    const static int EDGEMAP_TRI2_TO_QUAD[3];
    // for each quad vert, gives the triangle id and the triangle vert id
    const static std::array<std::pair<int, int>, 4> VERTMAP_QUAD_TO_TRI;
    // for each quad edge, gives the triangle id and triangle edge id
    const static std::array<std::pair<int, int>, 4> EDGEMAP_QUAD_TO_TRI;
    const static int INTERIOR_EDGE; // id of the edge of both triangles
                                    // that is not an edge of the quad

  private:
    // returns the vertex index on the triangle that is on the given quad
    // vertex id.  Returns -1 if not found
    int find_common_vertex(const int* vertmap, int quadId);

    std::shared_ptr<mesh::Mesh> m_mesh;
    TriangleCoordUtils m_triangleUtils;
};

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
