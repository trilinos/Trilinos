#include "adjacency_search.hpp"
#include "mesh_entity.hpp"
#include "plane_projection.hpp"
#include <cassert>
#include <limits>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

mesh::MeshEntityPtr AdjacencySearch::get_next(std::vector<mesh::MeshEntityPtr>& els2)
{
  if (m_mesh1Queue.size() == 0)
    return nullptr;

  std::pair<mesh::MeshEntityPtr, mesh::MeshEntityPtr> p = m_mesh1Queue.front();
  m_mesh1Queue.pop();
  get_next(p.first, p.second, els2);
  assert(els2.size() > 0);

  return p.first;
}

std::pair<mesh::MeshEntityPtr, mesh::MeshEntityPtr> AdjacencySearch::get_first_element()
{
  // find vertex with minimum number of upward edges (which implies minimum
  // number of elements)
  mesh::MeshEntityPtr v1 = nullptr;
  int nV1                = std::numeric_limits<int>::max();
  for (auto& v : m_mesh1->get_vertices())
    if (v && v->count_up() < nV1)
    {
      v1  = v;
      nV1 = v->count_up();
    }

  // find the vertex on mesh2 closest to v1;
  utils::Point p1        = v1->get_point_orig(0);
  mesh::MeshEntityPtr v2 = nullptr;
  double dist            = std::numeric_limits<double>::max();
  for (auto& v : m_mesh2->get_vertices())
  {
    double distV = compute_d2(p1, v->get_point_orig(0));
    if (distV < dist)
    {
      v2   = v;
      dist = distV;
    }
  }

  // (arbitrarily) pick the first upward element of each vertex
  std::vector<mesh::MeshEntityPtr> entities;
  [[maybe_unused]] int nUp = get_upward(v1, 2, entities);
  assert(nUp > 0);
  mesh::MeshEntityPtr el1 = entities[0];

  assert(v2);
  nUp = get_upward(v2, 2, entities);
  assert(nUp > 0);
  mesh::MeshEntityPtr el2 = entities[0];

  return {el1, el2};
}

void AdjacencySearch::get_next(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2, std::vector<mesh::MeshEntityPtr>& els2)
{
  if (m_output)
  {
    std::cout << "\nEntered getNext" << std::endl;
    std::cout << "el1 = " << el1 << std::endl;
    std::cout << "el2 = " << el2 << std::endl;
  }
  // do search to find el2 that contains el1's centroid

  // find el2 that contains el1s centroid

  // if (any el2 vertex within el1)
  //   get all vertices within el1
  //   get elements from vertices
  //   return
  // else if (any el2 edges intersect el1's edges)
  //   do edge-adjacent element search
  //   add any elements whose centroids are contained within el1 or whose
  //   edges intersect el1
  //   return
  // else // el2 contains el1
  //   add el2 to list
  //   return

  // get element that contains the centroid of el1
  el2 = get_mesh2_element(el1, el2);

  if (m_output)
  {
    std::cout << "new el2 = " << el2 << std::endl;
  }

  mesh::MeshEntityPtr v, edge2;
  if ((v = m_preds.any_vertices_contained(el1, el2)))
  {
    if (m_output)
      std::cout << "found contained vertices" << std::endl;
    SetType<mesh::MeshEntityPtr> verts;
    get_contained_vertices(el1, v, verts);
    get_elements_from_vertices(el1, verts, els2);
  } else if ((edge2 = m_preds.any_edges_intersect(el1, el2)))
  {
    if (m_output)
      std::cout << "getting intersecting edge" << std::endl;
    do_edge_search(el1, edge2, els2);
  } else // el2 contains el1
  {
    if (m_output)
      std::cout << "el2 contains el1" << std::endl;

    els2.clear();
    els2.push_back(el2);
  }

  update_mesh1_queue(el1, els2);
}

mesh::MeshEntityPtr AdjacencySearch::get_mesh2_element(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr el2)
{
  // TODO: this search adds elements in all directions, not just those towards
  //       el1's centroid.  This means the queue will continually expand away
  //       from the target element until the target element is found
  //       Unfortunately, we can't guarantee any kind of montonicity in the
  //       approach to the target element, so its hard to filter out elements
  utils::Point centroid1 = compute_centroid_3d(el1);
  std::vector<mesh::MeshEntityPtr> entities;
  std::queue<mesh::MeshEntityPtr> que;
  SetType<mesh::MeshEntityPtr> seenEntities;

  que.push(el2);
  seenEntities.insert(el2);
  // mesh::MeshEntityPtr verts2[mesh::MaxDown];
  while (que.size() > 0)
  {
    mesh::MeshEntityPtr elI = que.front();
    que.pop();
    bool isContained = m_preds.is_point_contained(elI, centroid1);
    if (isContained)
    {
      return elI;
    }

    // else keep searching
    get_bridge_adjacent(elI, 1, 2, entities);
    for (auto& e : entities)
      if (seenEntities.count(e) == 0)
      {
        que.push(e);
        seenEntities.insert(e);
      }
  }

  throw std::runtime_error("could not find an element on mesh2 that contains "
                           "the centroid of el1");
}

void AdjacencySearch::get_contained_vertices(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr v,
                                             SetType<mesh::MeshEntityPtr>& verts)
{
  std::queue<mesh::MeshEntityPtr> que;
  verts.clear();

  que.push(v);
  verts.insert(v);
  //TODO: the same exterior vert can be processed more than once.  This is slightly inefficient
  while (que.size() > 0)
  {
    auto vI = que.front();
    que.pop();
    for (int j = 0; j < vI->count_up(); ++j)
    {
      auto edgeJ = vI->get_up(j);
      auto elJ   = edgeJ->get_up(0);
      auto vJ    = get_other_vert(vI, edgeJ);

      if (verts.count(vJ) == 0 && m_preds.is_point_contained(el1, elJ, vJ->get_point_orig(0)))
      {
        que.push(vJ);
        verts.insert(vJ);
      }     
    }
  }

  if (m_output)
  {
    std::cout << "contained vertices: " << std::endl;
    for (auto& vert : verts)
      std::cout << "  " << vert << std::endl;
  }
}

void AdjacencySearch::get_elements_from_vertices(mesh::MeshEntityPtr el1, const SetType<mesh::MeshEntityPtr>& verts,
                                                 std::vector<mesh::MeshEntityPtr>& els)
{
  if (m_output)
    std::cout << "\nEntered get_elementsFromVertices" << std::endl;

  SetType<mesh::MeshEntityPtr> elsUnique;
  for (auto& v : verts)
  {
    if (m_output)
      std::cout << "v = " << v << std::endl;
    get_upward(v, 2, els);

    elsUnique.insert(els.begin(), els.end());
  }

  if (m_output)
    std::cout << "number of contained verts = " << verts.size() << std::endl;

  // search for edges that pass through el1 but don't have any vertices inside
  //TODO: this checks every edge, including those between 2 interior points, which cant possible
  //      intersect an edge of el1
  std::queue<mesh::MeshEntityPtr> que;
  for (auto& el : elsUnique)
    for (int i = 0; i < el->count_down(); ++i)
      que.push(el->get_down(i));

  while (que.size() > 0)
  {
    auto edge = que.front();
    que.pop();

    for (int i = 0; i < edge->count_up(); ++i)
    {
      if (elsUnique.count(edge->get_up(i)) == 0)
      {
        if (m_preds.any_edge_intersect(el1, edge))
        {
          if (m_output)
            std::cout << "found edge with intersection: " << edge << std::endl;
          auto el = edge->get_up(i);
          elsUnique.insert(el);
          for (int j = 0; j < el->count_down(); ++j)
            if (el->get_down(j) != edge)
              que.push(el->get_down(j));
        }
      }
    }
  }

  if (m_output)
    std::cout << "finished edge search" << std::endl;

  els.assign(elsUnique.begin(), elsUnique.end());
}

void AdjacencySearch::do_edge_search(mesh::MeshEntityPtr el1, mesh::MeshEntityPtr edge2,
                                     std::vector<mesh::MeshEntityPtr>& els)
{
  std::queue<mesh::MeshEntityPtr> que;
  SetType<mesh::MeshEntityPtr> seenEls;
  std::vector<mesh::MeshEntityPtr> entities;

  for (int i = 0; i < edge2->count_up(); ++i)
  {
    mesh::MeshEntityPtr elI = edge2->get_up(i);
    que.push(elI);
    seenEls.insert(elI);
  }

  els.clear();
  while (que.size() > 0)
  {
    mesh::MeshEntityPtr elI = que.front();
    que.pop();

    if (m_output)
      std::cout << "\nconsidering element " << elI << std::endl;

    // check if the centroid of el_i is contained in el1, or if
    // any of el_is edges intersect el1
    if (m_preds.is_point_contained(el1, elI, compute_centroid_3d(elI)) || m_preds.any_edges_intersect(el1, elI))
    {
      if (m_output)
        std::cout << "adding element" << std::endl;
      els.push_back(elI);

      // do an edge adjacent search for more elements
      // Could also do a vertex adjacent search, but doing edge based seems
      // less likely to add unneeded elements to the que (which would result
      // in more calls to anyEdgesIntersect, which could be expensive)
      get_bridge_adjacent(elI, 1, 2, entities);
      for (auto& e : entities)
        if (seenEls.count(e) == 0)
        {
          que.push(e);
          seenEls.insert(e);
        }
    } else
    {
      if (m_output)
        std::cout << "not adding element" << std::endl;
    }
  }
}

void AdjacencySearch::update_mesh1_queue(mesh::MeshEntityPtr el1, const std::vector<mesh::MeshEntityPtr>& els)
{
  std::vector<mesh::MeshEntityPtr> entities;
  get_bridge_adjacent(el1, 0, 2, entities);
  for (auto& e : entities)
    if (!(*m_field)(e, 0, 0))
    {
      mesh::MeshEntityPtr el2 = get_nearest_element(e, els);
      (*m_field)(e, 0, 0)     = true;
      m_mesh1Queue.push(std::make_pair(e, el2));
    }
}

mesh::MeshEntityPtr AdjacencySearch::get_nearest_element(mesh::MeshEntityPtr el1,
                                                         const std::vector<mesh::MeshEntityPtr>& els2)
{
  utils::Point p1         = compute_centroid_3d(el1);
  mesh::MeshEntityPtr el2 = nullptr;
  double dist             = std::numeric_limits<double>::max();

  for (auto& e : els2)
  {
    double distE = compute_d2(p1, compute_centroid_3d(e));
    if (distE < dist)
    {
      el2  = e;
      dist = distE;
    }
  }

  return el2;
}

double compute_d2(const utils::Point& p1, const utils::Point& p2)
{
  double dx = p2.get_x() - p1.get_x();
  double dy = p2.get_y() - p1.get_y();
  double dz = p2.get_z() - p1.get_z();

  return dx * dx + dy * dy + dz * dz;
}

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
