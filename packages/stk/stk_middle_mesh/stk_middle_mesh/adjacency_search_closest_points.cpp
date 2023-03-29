#include "adjacency_search_closest_points.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void AdjacencySearchClosestPoints::search(MeshEntityPtr startVertex, int numVertices,
                                          std::vector<MeshEntityPtr>& closestVerts)
{
  closestVerts.clear();
  DistanceSortedSet candidateVerts(DistanceCompare(startVertex->get_point_orig(0)));
  SetType seenVerts;
  seenVerts.insert(startVertex);
  add_one_layer_surrounding_vertices(startVertex, seenVerts, candidateVerts);

  while (!(candidateVerts.size() == 0) && closestVerts.size() < size_t(numVertices))
  {
    auto it                   = candidateVerts.begin();
    MeshEntityPtr currentVert = *it;
    candidateVerts.erase(it);

    if (is_vertex_included(currentVert))
      closestVerts.push_back(currentVert);

    add_two_layers_surrounding_vertices(currentVert, seenVerts, candidateVerts);
  }
}

void AdjacencySearchClosestPoints::add_one_layer_surrounding_vertices(MeshEntityPtr vert, SetType& seenVerts,
                                                                      DistanceSortedSet& candidateVerts)
{
  std::vector<MeshEntityPtr> surroundingVerts;
  add_surrounding_vertices(vert, seenVerts, candidateVerts, surroundingVerts);
}

void AdjacencySearchClosestPoints::add_two_layers_surrounding_vertices(MeshEntityPtr vert, SetType& seenVerts,
                                                                       DistanceSortedSet& candidateVerts)
{
  std::vector<MeshEntityPtr> surroundingVerts, surroundingVertsUnused;
  add_surrounding_vertices(vert, seenVerts, candidateVerts, surroundingVerts);
  for (auto& surroundingVert : surroundingVerts)
    add_surrounding_vertices(surroundingVert, seenVerts, candidateVerts, surroundingVertsUnused);
}

void AdjacencySearchClosestPoints::add_surrounding_vertices(MeshEntityPtr vert, SetType& seenVerts,
                                                            DistanceSortedSet& candidateVerts,
                                                            std::vector<MeshEntityPtr>& surroundingVerts)
{
  get_bridge_adjacent(vert, 2, 0, surroundingVerts);
  for (auto& adjacentVert : surroundingVerts)
    if (seenVerts.count(adjacentVert) == 0)
    {
      candidateVerts.insert(adjacentVert);
      seenVerts.insert(adjacentVert);
    }
}

bool AdjacencySearchClosestPoints::is_vertex_included(MeshEntityPtr vert)
{
  return !m_vertexStatus || (*m_vertexStatus)(vert, 0, 0);
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
