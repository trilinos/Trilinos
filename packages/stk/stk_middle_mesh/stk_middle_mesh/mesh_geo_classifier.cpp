#include "mesh_geo_classifier.hpp"
#include <algorithm>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void MeshGeoClassifier::classify(std::shared_ptr<Mesh> mesh)
{
  set_all_interior(mesh);
  set_boundary_edge_dim(mesh);
  auto verts = set_vertex_verts(mesh);
  set_boundary_edge_ids(mesh, verts);
  set_boundary_vert_ids(mesh);
}

// classify all mesh entities on geometric face 0
void MeshGeoClassifier::set_all_interior(std::shared_ptr<Mesh> mesh)
{
  GeoClassification vertG(2, 0);
  for (int dim = 0; dim <= 2; ++dim)
    for (auto& e : mesh->get_mesh_entities(dim))
      if (e)
        mesh->set_geo_class(e, vertG);
}

// sets the dimension (but not id) of boundary edges
void MeshGeoClassifier::set_boundary_edge_dim(std::shared_ptr<Mesh> mesh)
{
  GeoClassification edgeG(1, std::numeric_limits<int>::max());
  for (auto& edge : mesh->get_edges())
    if (edge && edge->count_up() == 1)
      mesh->set_geo_class(edge, edgeG);
}

// classify vertices on geometric vertices
// returns vector of mesh verticies classified on geometric vertices
std::vector<MeshEntityPtr> MeshGeoClassifier::set_vertex_verts(std::shared_ptr<Mesh> mesh)
{
  // do these in a deterministic order so that classifying several meshes
  // on the same geometry will give the same result, regardless of the
  // order the verts were created
  std::vector<MeshEntityPtr> verts;
  for (auto& v : mesh->get_vertices())
    if (v && v->count_up() == 2)
      verts.push_back(v);

  // sort
  auto cmp = [](const MeshEntityPtr& v1, const MeshEntityPtr v2) {
    auto pt1 = v1->get_point_orig(0);
    auto pt2 = v2->get_point_orig(0);
    if (pt1.x < pt2.x)
      return true;
    else if (pt1.x > pt2.x)
      return false;
    else if (pt1.y < pt2.y)
      return true;
    else if (pt1.y > pt2.y)
      return false;
    else if (pt1.z < pt2.z)
      return true;
    else if (pt1.z > pt2.z)
      return false;
    else
      return false;
  };

  std::sort(verts.begin(), verts.end(), cmp);

  for (unsigned int i = 0; i < verts.size(); ++i)
    mesh->set_geo_class(verts[i], GeoClassification(0, i));

  return verts;
}

void MeshGeoClassifier::set_boundary_edge_ids(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts)
{
  int iter        = 0;
  MeshEntityPtr v = get_start_vertex(mesh, verts);
  GeoClassification(1, iter);

  do
  {
    MeshEntityPtr edge = get_next_edge(mesh, v, iter);
    mesh->set_geo_class(edge, GeoClassification(1, iter));

    if (edge->get_down(0) == v)
      v = edge->get_down(1);
    else
      v = edge->get_down(0);

    // found a geometric vertex, which divides geometric edges
    if (mesh->get_geo_class(v).dim == 0)
      iter++;

  } while (v != verts[0]);
}

MeshEntityPtr MeshGeoClassifier::get_start_vertex(std::shared_ptr<Mesh> mesh, std::vector<MeshEntityPtr>& verts)
{
  if (verts.size() > 0)
    return verts[0];

  // else pick an a vertex on a boundary edge at random
  for (auto& edge : mesh->get_edges())
    if (edge && mesh->get_geo_class(edge).dim == 1)
      return edge->get_down(0);

  throw std::runtime_error("could not find starting entity");
}

// returns a mesh edge that is:
//   1. adjacent to v
//   2. is classified on a geometric edge
//   3. has geometric edge id greater than iter
//   4. is closest parallel with the vector <1, 0, 0>
MeshEntityPtr MeshGeoClassifier::get_next_edge(std::shared_ptr<Mesh> mesh, MeshEntityPtr v, const int iter)
{
  utils::Point dirVec(1, 0, 0);
  utils::Point ptI = v->get_point_orig(0);
  double minDev    = std::numeric_limits<double>::max();
  int minId        = -1;
  for (int i = 0; i < v->count_up(); ++i)
  {
    auto edge = v->get_up(i);
    auto g    = mesh->get_geo_class(edge);
    if (g.dim != 1 || g.id <= iter)
      continue;

    utils::Point otherPt;
    if (edge->get_down(0) == v)
      otherPt = edge->get_down(1)->get_point_orig(0);
    else
    {
      assert(edge->get_down(1) == v);
      otherPt = edge->get_down(0)->get_point_orig(0);
    }

    auto bVec   = otherPt - ptI;
    double devI = std::abs(dot(dirVec, bVec) / std::sqrt(dot(bVec, bVec)) - 1);
    if (devI < minDev)
    {
      minDev = devI;
      minId  = i;
    }
  }

  assert(minId != -1);
  return v->get_up(minId);
}

// classifies verts on geometric edges
void MeshGeoClassifier::set_boundary_vert_ids(std::shared_ptr<Mesh> mesh)
{
  for (auto& edge : mesh->get_edges())
    if (edge)
    {
      auto g = mesh->get_geo_class(edge);
      if (g.dim == 1)
      {
        auto v1 = edge->get_down(0);
        auto v2 = edge->get_down(1);

        if (mesh->get_geo_class(v1).dim != 0)
          mesh->set_geo_class(v1, g);

        if (mesh->get_geo_class(v2).dim != 0)
          mesh->set_geo_class(v2, g);
      }
    }
}

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
