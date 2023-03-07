#include "mesh_snapper.hpp"

#include <limits>

#include "mesh_io.hpp" //TODO: DEBUGGING
#include <iostream>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void MeshSnapper::snap(std::shared_ptr<Mesh> mesh1, std::shared_ptr<Mesh> mesh2)
{
  m_classifier = std::make_shared<predicates::impl::PointClassifierNormalWrapper>(
      mesh2, middle_mesh::impl::PointClassifierNormalWrapperTolerances(m_opts.tol, m_opts.tol));
  AdjacencySearch search(mesh1, mesh2);
  std::vector<MeshEntityPtr> mesh2Els;
  MeshEntityPtr el1;

  while ((el1 = search.get_next(mesh2Els)))
  {
    snap_el(el1, mesh2Els);
  }
}

void MeshSnapper::snap_el(MeshEntityPtr el1, const std::vector<MeshEntityPtr> mesh2Els)
{
  // get unique vertices
  SetType<MeshEntityPtr> verts2;
  MeshEntityPtr vertsI[MAX_DOWN];
  for (auto& el2 : mesh2Els)
  {
    int nverts = get_downward(el2, 0, vertsI);
    for (int i = 0; i < nverts; ++i)
      verts2.insert(vertsI[i]);
  }

  MeshEntityPtr verts1[MAX_DOWN], edges1[MAX_DOWN];
  get_downward(el1, 0, verts1);
  get_downward(el1, 1, edges1);
  for (auto& v : verts2)
  {
    MeshEntityPtr el2               = v->get_up(0)->get_up(0);
    predicates::impl::PointRecord r = m_classifier->classify(el1, el2, v->get_point_orig(0));
    if (r.type == predicates::impl::PointClassification::Vert)
    {
      v->set_point_orig(0, verts1[r.id]->get_point_orig(0));
    } else if (r.type == predicates::impl::PointClassification::Edge)
    {
      auto edge   = edges1[r.id];
      auto edgeXi = m_classifier->get_edge_xi(r);
      v->set_point_orig(0, compute_edge_coords_orig(edge, edgeXi));
    }
  }
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
