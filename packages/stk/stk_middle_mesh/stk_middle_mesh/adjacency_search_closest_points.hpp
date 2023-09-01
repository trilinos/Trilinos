#ifndef ADJACENCY_SEARCH_CLOSEST_POINTS_H
#define ADJACENCY_SEARCH_CLOSEST_POINTS_H

#include "field.hpp"
#include "mesh.hpp"

#include <set>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

// TODO: finish this
struct DistanceCompare
{
    DistanceCompare(const utils::Point& pt)
      : m_pt(pt)
    {}

    bool operator()(const MeshEntityPtr& lhs, const MeshEntityPtr& rhs) const
    {
      utils::Point ptLhs = lhs->get_point_orig(0);
      utils::Point ptRhs = rhs->get_point_orig(0);

      double distLhs = dot(ptLhs - m_pt, ptLhs - m_pt);
      double distRhs = dot(ptRhs - m_pt, ptRhs - m_pt);

      return distLhs < distRhs;
    }

  private:
    utils::Point m_pt;
};

// Does an adjacency-based search to find the mesh vertices closest (in xyz space)
// do the given vertex.  Optionally allows for returning a filtered set of vertices
// via the vertex_status field.  Only vertices for which the field returns true will
// be returned
class AdjacencySearchClosestPoints
{
    using DistanceSortedSet = std::multiset<MeshEntityPtr, DistanceCompare>;
    using SetType           = std::set<MeshEntityPtr, MeshEntityCompare>;

  public:
    using Bool = int_least8_t;

    AdjacencySearchClosestPoints(std::shared_ptr<Mesh> mesh, mesh::FieldPtr<Bool> vertexStatus = nullptr)
      : m_mesh(mesh)
      , m_vertexStatus(vertexStatus)
    {}

    void search(MeshEntityPtr startVertex, int numVertices, std::vector<MeshEntityPtr>& closestVerts);

  private:
    void add_one_layer_surrounding_vertices(MeshEntityPtr vert, SetType& seenVerts, DistanceSortedSet& candidateVerts);

    void add_two_layers_surrounding_vertices(MeshEntityPtr vert, SetType& seenVerts, DistanceSortedSet& candidateVerts);
    void add_surrounding_vertices(MeshEntityPtr verts, SetType& seenVerts, DistanceSortedSet& candidateVerts,
                                  std::vector<MeshEntityPtr>& surroundingVerts);

    bool is_vertex_included(MeshEntityPtr vert);

    std::shared_ptr<Mesh> m_mesh;
    mesh::FieldPtr<Bool> m_vertexStatus;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
